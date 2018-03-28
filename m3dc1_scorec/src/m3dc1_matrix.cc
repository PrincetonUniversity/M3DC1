/******************************************************************************
  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifdef M3DC1_PETSC
#include "m3dc1_matrix.h"
#include "m3dc1_matrix_allocate.h"
#include "m3dc1_mesh.h"
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <petscsys.h>
#include <cassert>
#include <iostream>
#include <vector>
#ifdef PETSC_USE_COMPLEX
#include <complex>
using std::complex;
#endif
void las_init(int * argc, char ** argv[], MPI_Comm cm)
{
  PETSC_COMM_WORLD = cm;
  PetscInitialize(argc,argv,PETSC_NULL,PETSC_NULL);
}
void get_num_blocks(m3dc1_matrix * mat, int * ent_dim, int * eid, int * num_blks)
{
  // assuming only verts have nodes
  int vrt_dim = 0;
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(ent_dim,eid,&vrt_dim,&num_adj_vrts);
  int blks_per_nd = mat->get_field()->get_num_value();
  *num_blks = blks_per_nd * num_adj_vrts;
}
void get_block_ids(m3dc1_matrix * mat, int * ent_dim, int * eid, int * blk_ids)
{
  int is_par = mat->is_parallel();
  int vrt_dim = 0;
  m3dc1_field * fld = mat->get_field();
  // assuming only vrts have nodes
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(ent_dim,eid,&vrt_dim,&num_adj_vrts);
  int * adj_vrts = new int[num_adj_vrts];
  m3dc1_ent_getadj(ent_dim,eid,&vrt_dim,&adj_vrts[0],&num_adj_vrts,&num_adj_vrts);
  int dofs_per_blk = fld->get_dof_per_value();
  int blks_per_nd = fld->get_num_value();
  int dofs_per_nd = dofs_per_blk * blks_per_nd;
  int * lcl_dof_ids = new int[dofs_per_nd];
  int * gbl_dof_ids = new int[dofs_per_nd];
  int * dof_ids[] = {&lcl_dof_ids[0],&gbl_dof_ids[0]};
  int dof_cnt = 0;
  int nds_per_ent = num_adj_vrts;
  for(int ent_nd = 0; ent_nd < nds_per_ent; ++ent_nd)
  {
    get_ent_localdofid(fld,adj_vrts[ent_nd],&lcl_dof_ids[0],&dof_cnt);
    get_ent_globaldofid(fld,adj_vrts[ent_nd],&gbl_dof_ids[0],&dof_cnt);
    for(int nd_blk = 0; nd_blk < blks_per_nd; ++nd_blk)
    {
      int & dof_id = dof_ids[is_par][nd_blk*dofs_per_nd];
      int blk_id = dof_id / dofs_per_blk;
      assert(dof_id % dofs_per_blk == 0);
      blk_ids[ent_nd * blks_per_nd + nd_blk] = blk_id;
    }
  }
  delete [] gbl_dof_ids;
  delete [] lcl_dof_ids;
  delete [] adj_vrts;
}
// todo : account for complex type
void insert_element_blocks(m3dc1_matrix * mat, int * ent_dim, int * eid, double * vals)
{
  int num_ent_blks = 0;
  get_num_blocks(mat,ent_dim,eid,&num_ent_blks);
  int * blk_ids = new int[num_ent_blks];
  // DBG(memset());
  get_block_ids(mat,ent_dim,eid,&blk_ids[0]);
  // assuming only verts hold nodes
  int num_adj_vrts = 0;
  int vrt_dim = 0;
  m3dc1_ent_getnumadj(ent_dim,eid,&vrt_dim,&num_adj_vrts);
  int blks_per_nd = mat->get_field()->get_num_value();
  int blks_per_elt = blks_per_nd * num_adj_vrts;
  mat->add_blocks(blks_per_elt,blk_ids,blks_per_elt,blk_ids,vals);
  delete [] blk_ids;
}
void insert_node_blocks(m3dc1_matrix * mat, int * ent_dim, int * eid, int * nd1, int * nd2, double * vals)
{
  int vrt_dim = 0;
  int blks_per_ent = 0;
  get_num_blocks(mat,ent_dim,eid,&blks_per_ent);
  int * blk_ids = new int[blks_per_ent];
  get_block_ids(mat,ent_dim,eid,&blk_ids[0]);
  //DBG(memset());
  //assuming only verts have nodes and that all nodes have the same number of blocks(values..)
  int blks_per_nd = mat->get_field()->get_num_value();
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(ent_dim,eid,&vrt_dim,&num_adj_vrts);
  assert(num_adj_vrts * blks_per_nd == blks_per_ent);
  mat->add_blocks(blks_per_nd,&blk_ids[*nd1 * blks_per_nd],blks_per_nd,&blk_ids[*nd2 * blks_per_nd],vals);
  delete [] blk_ids;
}
// ***********************************
//              HELPER
// ***********************************
void describeMatrix(Mat A)
{
  MatInfo info;
  MatGetInfo(A,MAT_LOCAL,&info);
  MatType tp;
  MatGetType(A,&tp);
  int gbl_rws = 0;
  int gbl_cls = 0;
  MatGetSize(A,&gbl_rws,&gbl_cls);
  int lcl_rws = 0;
  int lcl_cls = 0;
  MatGetLocalSize(A,&lcl_rws,&lcl_cls);
  int fst_rw = 0;
  int lst_rw = 0;
  MatGetOwnershipRange(A,&fst_rw,&lst_rw);
  int bs = sqrt(info.block_size);
  std::cout << "Matrix created[" << PCU_Comm_Self() << "]: " << std::endl
            << "\t Mat type: " << tp << std::endl
            << "\t Block size: " << bs << std::endl
            << "\t Global sizes: (" << gbl_rws << ", " << gbl_cls << ")" << std::endl
            << "\t Local sizes: (" << lcl_rws << ", " << lcl_cls << ")" << std::endl
            << "\t Ownership range: [" << fst_rw << ", " << lst_rw << ")" << std::endl
            << "\t Global sizes (blocks): (" << gbl_rws / bs << ", " << gbl_cls / bs << ")" << std::endl
            << "\t Local sizes (blocks): (" << lcl_rws / bs << ", " << lcl_cls / bs << ")" << std::endl
            << "\t Ownership range (blocks): [" << fst_rw / bs << ", " << lst_rw / bs << ")" << std::endl
            << "\t NNZ allocated: " << info.nz_allocated << std::endl;
}
void printMemStat()
{
  PetscLogDouble mem, mem_max;
  PetscMemoryGetCurrentUsage(&mem);
  PetscMemoryGetMaximumUsage(&mem_max);
  std::cout << "\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank " << PCU_Comm_Self() << " current " << mem/1e6 << std::endl;
}
void field2Vec(MPI_Comm cm, m3dc1_field * fld, Vec V, int st)
{
  bool lcl = cm == MPI_COMM_SELF;
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // exernal var
  apf::Mesh2 * msh = m3dc1_msh->get_mesh();
  int num_own_nds = m3dc1_msh->get_own_count(0); // assuming only verts have dofs
  int num_own_dof = 0;
  FieldID fid = fld->get_id();
  m3dc1_field_getnumowndof(&fid, &num_own_dof);
  int dof_per_nd = num_own_dof / num_own_nds;
  VecCreateMPI(cm,num_own_dof,PETSC_DETERMINE,&V);
  int num_lcl_nds = m3dc1_msh->get_local_count(0);
  int * dof_ids = new int[dof_per_nd];
  memset(&dof_ids[0],0,sizeof(int)*dof_per_nd);
  // assumes st = 0 means real, st=1 means complex
  int sz = dof_per_nd * (st+1);
  double * dof_data = new double[sz];
  //DBG(memset(&dof_data[0],0,sizeof(double)*sz));
#ifdef PETSC_USE_COMPLEX
  std::vector<PetscComplex> cplx_data(dof_per_nd);
#endif
  int vrt_tp = 0;
  for(int nd = 0; nd < num_lcl_nds; ++nd)
  {
    apf::MeshEntity * ent = apf::getMdsEntity(msh,vrt_tp,nd);
    if(!msh->isOwned(ent) && !lcl)
      continue;
    int num_dof = 0;
    m3dc1_ent_getdofdata(&vrt_tp,&nd,&fid,&num_dof,&dof_data[0]);
    m3dc1_ent_getglobaldofid(&vrt_tp,&nd,&fid,&dof_ids[0],&num_dof);
#ifdef PETSC_USE_COMPLEX
    for(int ii = 0; ii < dof_per_nd; ++ii)
      cplx_data[ii] = dof_data[ii*2] + dof_data[ii*2+1] * PETSC_i;
    VecSetValues(V,dof_per_nd,dof_ids,&cplx_data[0],INSERT_VALUES);
#else
    VecSetValues(V,dof_per_nd,dof_ids,&dof_data[0],INSERT_VALUES);
#endif
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  delete [] dof_ids;
  delete [] dof_data;
}
void vec2Field(MPI_Comm cm, m3dc1_field * fld, Vec V, int st)
{
  bool lcl = cm == MPI_COMM_SELF;
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // external variable
  apf::Mesh2 * msh = m3dc1_msh->get_mesh();
  FieldID fid = fld->get_id();
  int num_own_nds = m3dc1_msh->get_own_count(0);
  int num_own_dof = 0;
  m3dc1_field_getnumowndof(&fid,&num_own_dof);
  int vrt_tp = 0;
  int dof_per_nd = num_own_dof / num_own_nds;
  int num_lcl_nd = m3dc1_msh->get_local_count(0);
  int * dof_ids = new int[dof_per_nd];
  //DBG(memset(&dof_ids[0],0,sizeof(int)*dof_per_nd));
  int sz = dof_per_nd*(st+1);
  double * dof_data = new double[sz];
#ifdef PETSC_USE_COMPLEX
  std::vector<PetscComplex> cplx_data(dof_per_nd);
#endif
  //DBG(memset(&dof_data[0],0,sizeof(double)*sz);
  for(int nd = 0; nd < num_lcl_nd; ++nd)
  {
    apf::MeshEntity * ent = apf::getMdsEntity(msh,vrt_tp,nd);
    if(!msh->isOwned(ent) && !lcl)
      continue;
    int dof_cnt = 0;
    m3dc1_ent_getglobaldofid(&vrt_tp,&nd,&fid,&dof_ids[0],&dof_cnt);
#ifdef PETSC_USE_COMPLEX
    VecGetValues(V,dof_per_ent,&dof_ids[0],&cplx_data[0]);
    for(int ii = 0; ii < dof_per_ent; ++ii)
    {
      dof_data[2*ii] = cplx_data[ii].real();
      dof_data[2*ii+1] = cplx_data[ii].imag();
    }
#else
    VecGetValues(V,dof_per_nd,&dof_ids[0],&dof_data[0]);
#endif
    m3dc1_ent_setdofdata(&vrt_tp,&nd,&fid,&dof_per_nd,&dof_data[0]);
  }
  // could get rid of this with some work
  m3dc1_field_sync(&fid);
  delete [] dof_ids;
  delete [] dof_data;
}
m3dc1_solver* m3dc1_solver::_instance=NULL;
m3dc1_solver* m3dc1_solver::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_solver();
  return _instance;
}
m3dc1_matrix::m3dc1_matrix(int i, int s, m3dc1_mesh * msh, m3dc1_field * f, MPI_Comm c)
  : cm(c)
  , A()
  , x()
  , b()
  , ksp()
  , id(i)
  , scalar_type(s)
  , fixed(false)
  , fld(f)
  , is_par(c != MPI_COMM_SELF)
  , blk_sz(0)
{
  const char * mat_tps[][2] = { {MATSEQAIJ, MATSEQBAIJ}, {MATMPIAIJ, MATMPIBAIJ} };
  bool is_par = cm != MPI_COMM_SELF;
  size_t num_ent[2] = { msh->get_mesh()->count(0), (size_t)apf::countOwned(msh->get_mesh(),0) }; // assumes that only verts hold dofs
  int blk_sz = fld->get_dof_per_value();
  int dof_per_ent = fld->get_num_value() * blk_sz;
  int num_lcl_dof = num_ent[is_par] * dof_per_ent;
  MatCreate(cm,&A);
  const char * mat_tp = mat_tps[is_par][blk_sz > 1];
  MatSetType(A,mat_tp);
  MatSetSizes(A,num_lcl_dof,num_lcl_dof,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetBlockSize(A,blk_sz);
  // call preallocate
  allocateMatrix(A,msh,fld);
  //MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  if(!(blk_sz-1)) // only supported for AIJ not BAIJ
    MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
  MatCreateVecs(A,&x,&b);
#ifdef DEBUG
  if(!PCU_Comm_Self())
    describeMatrix(A);
#endif
  KSPCreate(cm,&ksp);
  KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,0.000001, 0.00000001, PETSC_DEFAULT, 1000);
  if(msh->get_mesh()->getDimension() == 2)
  {
    KSPSetType(ksp,KSPPREONLY);
    PC pc;
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
  }
  KSPSetFromOptions(ksp);
}
m3dc1_matrix::~m3dc1_matrix()
{
  MatDestroy(&A);
}
void m3dc1_matrix::add_blocks(int blk_rw_cnt, int * blk_rws, int blk_col_cnt, int * blk_cols, double * vals)
{
  MatSetValuesBlocked(A,blk_rw_cnt,blk_rws,blk_col_cnt,blk_cols,vals,ADD_VALUES);
}
void m3dc1_matrix::fix()
{
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  fixed = true;
}
void m3dc1_matrix::get_values(std::vector<int>& rows, std::vector<int>& n_columns, std::vector<int>& columns, std::vector<double>& values)
{
#ifdef PETSC_USE_COMPLEX
   if (!PCU_Comm_Self())
     std::cout << "[M3DC1 ERROR] " << __func__ << ": not supported for complex\n";
   return M3DC1_FAILURE;
#else
  PetscErrorCode ierr;
  PetscInt rstart, rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  ierr = MatGetOwnershipRange(A, &rstart, &rend);
  for (PetscInt row=rstart; row<rend; ++row)
  {
    MatGetRow(A, row, &ncols, &cols, &vals);
    rows.push_back(row);
    n_columns.push_back(ncols);
    for (int i=0; i<ncols; ++i)
    {
      columns.push_back(cols[i]);
      values.push_back(vals[i]);
    }
    MatRestoreRow(A, row, &ncols, &cols, &vals);
  }
#endif
}
void m3dc1_matrix::add_values(int rsize, int * rows, int csize, int * cols, double * vals)
{
  MatSetValues(A, rsize, rows, csize, cols, vals, ADD_VALUES);
}
void m3dc1_matrix::set_values(int rsize, int * rows, int csize, int * cols, double * vals)
{
  MatSetValues(A, rsize, rows, csize, cols, vals, INSERT_VALUES);
}
void m3dc1_matrix::write(const char * fn)
{
  PetscViewer view;
  MPI_Comm cm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)A,&cm);
  PetscViewerASCIIOpen(cm, fn, &view);
  PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB);
  MatView(A, view);
  PetscViewerDestroy(&view);
}
void m3dc1_matrix::multiply(m3dc1_field * in, m3dc1_field * out)
{
  MPI_Comm cm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)A,&cm);
  field2Vec(cm,in,x,get_scalar_type());
  MatMult(A, x, b);
  vec2Field(cm,out,b,get_scalar_type());
}
void m3dc1_matrix::zero()
{
  MatZeroEntries(A);
  fixed = false;
};
void m3dc1_matrix::solve(m3dc1_field * rhs)
{
  MPI_Comm cm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)A, &cm);
  field2Vec(cm,rhs,b,get_scalar_type()); // copy field values into vector
  KSPSolve(ksp, b, x);
  int itr = -1;
  KSPGetIterationNumber(ksp, &itr);
  if(!PCU_Comm_Self())
    std::cout << "\t-- # solver iterations " << itr << std::endl;
  vec2Field(cm,rhs,x,get_scalar_type());
}
int m3dc1_matrix::solve_iteration_count()
{
  int itr = -1;
  KSPGetIterationNumber(ksp,&itr);
  return itr;
}
#endif //#ifndef M3DC1_MESHGEN
