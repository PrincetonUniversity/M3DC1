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

// function defined in m3dc1
// for internal test, define the two functions as blank
// add two lines in your main.cc to get rid of undefined symbol
// extern "C" int setPETScMat(int matrixid, Mat * A) {};
// extern "C" int setPETScKSP(int matrixid, KSP * ksp, Mat * A){};
extern "C" int setPETScMat(int matrixid, Mat * A);
#ifndef PETSCMASTER
// petsc 3.8.2 and older
extern "C" int setPETScKSP(int matrixid, KSP * ksp, Mat * A);
#else
extern "C" int setPETScKSP(int matrixid, KSP * ksp, Mat * A, Vec *b, Vec *x);
#endif
extern "C" int get_iter_num_(int * iter_num_p) {*iter_num_p=-1; return *iter_num_p;}

inline void double2petsc(const int n, double * in, PetscScalar *& out)
{
#ifdef PETSC_USE_COMPLEX
  for(int ii = 0; ii < n; ++ii)
    out[ii] = in[ii*2] + in[ii*2+1] * PETSC_i;
#else
  (void)n;
  out = in;
#endif
}
inline void petsc2double(const int n, PetscScalar * in, double *& out)
{
#ifdef PETSC_USE_COMPLEX
  for(int ii = 0; ii < n; ++ii)
  {
    out[ii*2] = in[ii].real();
    out[ii*2+1] = in[ii].imag();
  }
#else
  (void)n;
  out = in;
#endif
}
void las_init(int * argc, char ** argv[], MPI_Comm cm)
{
  PETSC_COMM_WORLD = cm;
  PetscInitialize(argc,argv,PETSC_NULL,PETSC_NULL);
}
void get_num_blocks(m3dc1_field * fld, int ent_dim, int eid, int * num_blks)
{
  // assuming only verts have nodes
  int vrt_dim = 0;
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(&ent_dim,&eid,&vrt_dim,&num_adj_vrts);
  int blks_per_nd = fld->getBlocksPerNode();
  *num_blks = blks_per_nd * num_adj_vrts;
}
void get_block_ids(m3dc1_field * fld,
                   bool lcl,
                   int ent_dim,
                   int eid,
                   int * blk_ids)
{
  // assuming only vrts have nodes
  int vrt_dim = 0;
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(&ent_dim,&eid,&vrt_dim,&num_adj_vrts);
  if(ent_dim == 0)
    num_adj_vrts = 1;
  int * adj_vrts = new int[num_adj_vrts];
  m3dc1_ent_getadj(&ent_dim,&eid,&vrt_dim,&adj_vrts[0],&num_adj_vrts,&num_adj_vrts);
  if(ent_dim == 0)
    adj_vrts[0] = eid;
  int dofs_per_blk = fld->getDofsPerBlock();
  int blks_per_nd = fld->getBlocksPerNode();
  int dofs_per_nd = dofs_per_blk * blks_per_nd;
  int * dof_ids = new int[dofs_per_nd];
  int dof_cnt = 0;
  int nds_per_ent = num_adj_vrts;
  for(int ent_nd = 0; ent_nd < nds_per_ent; ++ent_nd)
  {
    get_ent_xxxdofid(lcl,fld,adj_vrts[ent_nd],&dof_ids[0],&dof_cnt);
    for(int nd_blk = 0; nd_blk < blks_per_nd; ++nd_blk)
    {
      int & dof_id = dof_ids[nd_blk*dofs_per_blk];
      int blk_id = dof_id / dofs_per_blk;
      assert(dof_id % dofs_per_blk == 0);
      blk_ids[ent_nd * blks_per_nd + nd_blk] = blk_id;
    }
  }
  delete [] dof_ids;
  delete [] adj_vrts;
}
// todo : account for complex type
void insert_element_blocks(m3dc1_matrix * mat, int ent_dim, int eid, double * vals)
{
  bool lcl = !mat->is_parallel();
  m3dc1_field * fld = mat->get_field();
  int num_ent_blks = 0;
  get_num_blocks(fld,ent_dim,eid,&num_ent_blks);
  int * blk_ids = new int[num_ent_blks]();
  get_block_ids(fld,lcl,ent_dim,eid,&blk_ids[0]);
  // assuming only verts hold nodes
  int num_adj_vrts = 0;
  int vrt_dim = 0;
  m3dc1_ent_getnumadj(&ent_dim,&eid,&vrt_dim,&num_adj_vrts);
  int blks_per_nd = fld->getBlocksPerNode();
  int blks_per_elt = blks_per_nd * num_adj_vrts;
  mat->add_blocks(blks_per_elt,blk_ids,blks_per_elt,blk_ids,vals);
  delete [] blk_ids;
}
void insert_node_blocks(m3dc1_matrix * mat,
                        int ent_dim,
                        int eid,
                        int nd1,
                        int nd2,
                        double * vals)
{
  bool lcl = !mat->is_parallel();
  int vrt_dim = 0;
  int blks_per_ent = 0;
  m3dc1_field * fld = mat->get_field();
  get_num_blocks(fld,ent_dim,eid,&blks_per_ent);
  int * blk_ids = new int[blks_per_ent];
  get_block_ids(fld,lcl,ent_dim,eid,&blk_ids[0]);
  //assuming only verts have nodes and that all nodes have the same number of blocks (values..)
  int blks_per_nd = fld->getBlocksPerNode();
  int num_adj_vrts = 0;
  m3dc1_ent_getnumadj(&ent_dim,&eid,&vrt_dim,&num_adj_vrts);
  assert(num_adj_vrts * blks_per_nd == blks_per_ent);
  mat->add_blocks(blks_per_nd,&blk_ids[nd1 * blks_per_nd],blks_per_nd,&blk_ids[nd2 * blks_per_nd],vals);
  delete [] blk_ids;
}
void describeMatrix(Mat A)
{
  MPI_Comm cm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)A,&cm);
  MatInfo info;
  MatGetInfo(A,MAT_GLOBAL_SUM,&info);
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
  int rnk = -1;
  MPI_Comm_rank(cm,&rnk);
  if(rnk == 0)
  {
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
}
void printMemStat()
{
  PetscLogDouble mem, mem_max;
  PetscMemoryGetCurrentUsage(&mem);
  PetscMemoryGetMaximumUsage(&mem_max);
  std::cout << "\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank " << PCU_Comm_Self() << " current " << mem/1e6 << std::endl;
}
void write_vec(Vec v, const char * fn)
{
  PetscViewer view;
  MPI_Comm cm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)v,&cm);
  PetscViewerASCIIOpen(cm, fn, &view);
  PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB);
  VecView(v, view);
  PetscViewerDestroy(&view);
}
void field2Vec(m3dc1_field * fld, Vec V)
{
  MPI_Comm vcm;
  PetscObjectGetComm((PetscObject)V,&vcm);
  int sz = 0;
  MPI_Comm_size(vcm,&sz);
  bool lcl = sz == 1;
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // exernal var
  int num_lcl_vtx = m3dc1_msh->get_local_count(0);
  FieldID fid = fld->getId();
  int blk_per_vtx = fld->getBlocksPerNode();
  int dof_per_blk = fld->getDofsPerBlock();
  int dof_per_vtx = dof_per_blk * blk_per_vtx;
  int dof_sz = dof_per_vtx;
#ifdef PETSC_COMPLEX
  dof_sz *= 2;
#endif
  double * dof_data = new double[dof_sz]();
  PetscScalar * dof_data_cmplx = new PetscScalar[dof_sz]();
  PetscScalar * data = dof_data_cmplx;
  int vtx_dim = 0;
  // the number of blocks per vtx is const for the same field
  int * blk_ids = new int[blk_per_vtx]();
  for(int vtx = 0; vtx < num_lcl_vtx; ++vtx)
  {
    int num_dof = 0;
    m3dc1_ent_getdofdata(&vtx_dim,&vtx,&fid,&num_dof,&dof_data[0]);
    double2petsc(dof_per_vtx,dof_data,data);
    get_block_ids(fld,lcl,0,vtx,&blk_ids[0]);
    VecSetValuesBlocked(V,blk_per_vtx,blk_ids,&data[0],INSERT_VALUES);
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  delete [] blk_ids;
  delete [] dof_data_cmplx;
  delete [] dof_data;
}
void vec2Field(Vec V, m3dc1_field * fld)
{
  MPI_Comm vcm;
  PetscObjectGetComm((PetscObject)V,&vcm);
  int sz = 0;
  MPI_Comm_size(vcm,&sz);
  bool lcl = sz == 1;
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // external variable
  int num_lcl_vtxs = m3dc1_msh->get_local_count(0);
  FieldID fid = fld->getId();
  int blk_per_vtx = fld->getBlocksPerNode();
  int dof_per_blk = fld->getDofsPerBlock();
  int dof_per_vtx = dof_per_blk * blk_per_vtx;
  int dof_sz = dof_per_vtx;
#ifdef PETSC_COMPLEX
  dof_sz *= 2;
#endif
  double * dof_data = new double[dof_sz]();
  double * data = dof_data;
  PetscScalar * dof_data_cmplx = new PetscScalar[dof_sz]();
  int * dof_ids = new int[dof_per_vtx]();
  int vtx_dim = 0;
  int lw = -1;
  int hg = -1;
  VecGetOwnershipRange(V,&lw,&hg);
  for(int vtx = 0; vtx < num_lcl_vtxs; ++vtx)
  {
    int num_dof = 0;
    get_ent_xxxdofid(lcl,fld,vtx,dof_ids,&num_dof);
    if(dof_ids[0] >= lw && dof_ids[dof_per_vtx - 1] < hg)
    {
      VecGetValues(V,dof_per_vtx,dof_ids,&dof_data_cmplx[0]);
      petsc2double(dof_per_vtx,dof_data_cmplx,data);
      m3dc1_ent_setdofdata(&vtx_dim,&vtx,&fid,&dof_per_vtx,&data[0]);
    }
  }
  //MPI_Comm ocm = PCU_Get_Comm();
  //PCU_Switch_Comm(vcm);
  //m3dc1_field_sync(&fid);
  //PCU_Switch_Comm(ocm);
  delete [] dof_ids;
  delete [] dof_data_cmplx;
  delete [] dof_data;
}

int calcRowOwner(Mat A, int lcl_rws, int rw)
{
  MPI_Comm cm;
  PetscObjectGetComm((PetscObject)A,&cm);
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  int rnk = 0;
  MPI_Comm_rank(cm,&rnk);
  std::vector<int> rws(sz,0);
  MPI_Allgather(&lcl_rws,1,MPI_INTEGER,
                &rws[0],1,MPI_INTEGER,
                cm);
  int rw_hd = 0;
  int peer = 0;
  while(rw > rw_hd)
    rw_hd += rws[peer++];
  return --peer;
}
int getRowOwner(Mat A, int rw)
{
  MPI_Comm cm;
  PetscObjectGetComm((PetscObject)A,&cm);
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  const int * rngs = NULL;
  MatGetOwnershipRanges(A,&rngs);
  int rnk = 0;
  while(rw >= rngs[rnk])
  { ++rnk; }
  --rnk;
  return rnk;
}
bool ownRow(Mat A, int rw)
{
  int fst_rw = 0;
  int lst_rw = 0;
  MatGetOwnershipRange(A,&fst_rw,&lst_rw);
  return (rw <= fst_rw && rw < lst_rw);
}
bool willOwnRow(Mat A, int lcl_rws, int rw)
{
  int ownr = calcRowOwner(A,lcl_rws,rw);
  MPI_Comm cm;
  PetscObjectGetComm((PetscObject)A,&cm);
  int rnk = -1;
  MPI_Comm_rank(cm,&rnk);
  return ownr == rnk;
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
  , ownership(NULL)
  , low_row(-1)
  , hgh_row(-1)
  , non_lcl_nnz_cnt(0)
  , lcl_nnz_cnt(0)
{
  const char * mat_tps[][2] = { {MATSEQAIJ, MATSEQBAIJ}, {MATMPIAIJ, MATMPIBAIJ} };
  bool is_par = cm != MPI_COMM_SELF;
  size_t num_ent[2] = { msh->get_mesh()->count(0),
                        static_cast<size_t>(apf::countOwned(msh->get_mesh(),0,msh->get_ownership())) }; // assumes that only verts hold dofs
  blk_sz = fld->getDofsPerBlock();
  int dof_per_ent = fld->getBlocksPerNode() * blk_sz;
  int num_lcl_dof = num_ent[is_par] * dof_per_ent;
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  ownership = new int[sz];
  memset(&ownership[0],0,sz*sizeof(int));
  int num_lcl_blks = num_lcl_dof / blk_sz;
  MPI_Allgather(&num_lcl_blks,1,MPI_INTEGER,
                &ownership[0],1,MPI_INTEGER,
                cm);
  MatCreate(cm,&A);
  const char * mat_tp = mat_tps[is_par][blk_sz > 1];
  MatSetType(A,mat_tp);
  MatSetSizes(A,num_lcl_dof,num_lcl_dof,PETSC_DETERMINE,PETSC_DETERMINE);
  int num_gbl_dof = PETSC_DECIDE;
  PetscSplitOwnershipBlock(cm,blk_sz,&num_lcl_dof,&num_gbl_dof);
  MatSetBlockSize(A,blk_sz);
  allocateMatrix(A,this,msh,fld);
  MatSetUp(A);
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  if(!(blk_sz-1)) // only supported for AIJ not BAIJ
    MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
  MatCreateVecs(A,&x,&b);
#ifdef DEBUG
  describeMatrix(A);
#endif
  KSPCreate(cm,&ksp);
  KSPSetTolerances(ksp,0.000001, 0.00000001, PETSC_DEFAULT, 1000);
  if(msh->get_mesh()->getDimension() == 2)
  {
    KSPSetType(ksp,KSPPREONLY);
    PC pc;
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    //PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
    PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST);
  }
  KSPSetFromOptions(ksp);
  MatGetOwnershipRange(A,&low_row,&hgh_row);
  low_row /= blk_sz;
  hgh_row /= blk_sz;
}
m3dc1_matrix::~m3dc1_matrix()
{
  MatDestroy(&A);
}
void m3dc1_matrix::multiply(m3dc1_field * in, m3dc1_field * out)
{
  field2Vec(in,x);
  double nrm = 0.0;
  VecNorm(x,NORM_2,&nrm);
  PetscViewer x_vwr;
  PetscViewerASCIIOpen(cm,"mult_x.m",&x_vwr);
  VecView(x,x_vwr);
  PetscViewerDestroy(&x_vwr);
  if(!PCU_Comm_Self())
    std::cout << "Pre-multiply multiplicand norm is: " << nrm << std::endl;
  MatMult(A, x, b);
  VecScatter ctx;
  Vec lb;
  VecScatterCreateToAll(b,&ctx,&lb);
  VecScatterBegin(ctx,b,lb,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,b,lb,INSERT_VALUES,SCATTER_FORWARD);
  VecNorm(b,NORM_2,&nrm);
  if(!PCU_Comm_Self())
    std::cout << "Post-multiply result norm is: " << nrm << std::endl;
  PetscViewer b_vwr;
  PetscViewerASCIIOpen(cm,"mult_b.m",&b_vwr);
  VecView(b,b_vwr);
  PetscViewerDestroy(&b_vwr);
  vec2Field(b,out);
}
void m3dc1_matrix::solve(m3dc1_field * lhs)
{
  field2Vec(fld,b); // copy field values into vector
  PetscViewer rhs_vwr;
  PetscViewerASCIIOpen(cm,"solve_rhs.m",&rhs_vwr);
  VecView(b,rhs_vwr);
  PetscViewerDestroy(&rhs_vwr);
  double nrm = 0.0;
  VecNorm(b,NORM_2,&nrm);
  if(!PCU_Comm_Self())
    std::cout << "Pre-solve rhs norm: " << nrm << std::endl;
  KSPSetOperators(ksp,A,A);
  KSPSolve(ksp, b, x);
  int itr = -1;
  KSPGetIterationNumber(ksp, &itr);
  if(!PCU_Comm_Self())
    std::cout << "\t-- # solver iterations " << itr << std::endl;
  VecScatter ctx;
  Vec lx;
  VecScatterCreateToAll(x,&ctx,&lx);
  VecScatterBegin(ctx,x,lx,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,lx,INSERT_VALUES,SCATTER_FORWARD);
  VecNorm(x,NORM_2,&nrm);
  if(!PCU_Comm_Self())
    std::cout << "Post-solve lhs norm: " << nrm << std::endl;
  PetscViewer lhs_vwr;
  PetscViewerASCIIOpen(cm,"solve_lhs.m",&lhs_vwr);
  VecView(x,lhs_vwr);
  PetscViewerDestroy(&lhs_vwr);
  vec2Field(lx,lhs);
}
bool m3dc1_matrix::willOwn(int rw)
{
  int ownr = whoOwns(rw);
  int slf = -1;
  MPI_Comm_rank(cm,&slf);
  return ownr == slf;
}
int m3dc1_matrix::whoOwns(int rw)
{
  int rw_hd = 0;
  int peer = 0;
  while(rw_hd <= rw)
    rw_hd += ownership[peer++];
  return --peer;
}
int m3dc1_matrix::calcFirstRow()
{
  int slf = -1;
  MPI_Comm_rank(cm,&slf);
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  int frst_rw = 0;
  for(int rnk = 0; rnk < slf; ++rnk)
    frst_rw += ownership[rnk];
  return frst_rw;
}
int m3dc1_matrix::calcLastRowP1()
{
  int slf = -1;
  MPI_Comm_rank(cm,&slf);
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  int lst_rw = 0;
  for(int rnk = 0; rnk <= slf; ++rnk)
    lst_rw += ownership[rnk];
  return lst_rw;
}
void m3dc1_matrix::add_blocks(int blk_rw_cnt,
                              int * blk_rws,
                              int blk_col_cnt,
                              int * blk_cols,
                              double * vals)
{
  int sz = blk_sz * blk_sz * blk_rw_cnt * blk_col_cnt;
  PetscScalar * cplx_data = new PetscScalar[sz]();
  PetscScalar * data = cplx_data;
  double2petsc(sz,vals,data);
  MatSetValuesBlocked(A,blk_rw_cnt,blk_rws,blk_col_cnt,blk_cols,data,ADD_VALUES);
  for(int rw = 0; rw < blk_rw_cnt; ++rw)
  {
    if(blk_rws[rw] < low_row || blk_rws[rw] >= hgh_row)
      non_lcl_nnz_cnt += blk_sz * blk_sz * blk_col_cnt;
    else
      lcl_nnz_cnt += blk_sz * blk_sz * blk_col_cnt;
  }
  fixed = false;
  delete [] cplx_data;
}
void m3dc1_matrix::fix()
{
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  fixed = true;
}
void m3dc1_matrix::get_values(std::vector<int>& rows,
                              std::vector<int>& n_columns,
                              std::vector<int>& columns,
                              std::vector<double>& values)
{
#ifdef PETSC_USE_COMPLEX
  (void)rows;
  (void)n_columns;
  (void)columns;
  (void)values;
   if (!PCU_Comm_Self())
     std::cout << "[M3DC1 ERROR] " << __func__ << ": not supported for complex\n";
#else
  PetscInt rstart, rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  MatGetOwnershipRange(A, &rstart, &rend);
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
void m3dc1_matrix::add_values(int rw_cnt, int * rows, int cl_cnt, int * cols, double * vals)
{
  int sz = rw_cnt * cl_cnt;
  PetscScalar * cplx_data = new PetscScalar[sz]();
  PetscScalar * data = cplx_data;
  double2petsc(sz,vals,data);
  MatSetValues(A, rw_cnt, rows, cl_cnt, cols, data, ADD_VALUES);
  fixed = false;
  delete [] cplx_data;
}
void m3dc1_matrix::set_values(int rw_cnt, int * rows, int cl_cnt, int * cols, double * vals)
{
  int sz = rw_cnt * cl_cnt;
  PetscScalar * cplx_data = new PetscScalar[sz]();
  PetscScalar * data = cplx_data;
  double2petsc(sz,vals,data);
  MatSetValues(A, rw_cnt, rows, cl_cnt, cols, data, INSERT_VALUES);
  fixed = false;
  delete [] cplx_data;
}
void m3dc1_matrix::write(const char * fn)
{
  PetscViewer view;
  PetscViewerASCIIOpen(cm, fn, &view);
  PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB);
  MatView(A, view);
  PetscViewerDestroy(&view);
}
void m3dc1_matrix::printAssemblyInfo()
{
  int d_stsh = 0;
  int d_mllc = 0;
  int b_stsh = 0;
  int b_mllc = 0;
  MatStashGetInfo(A,&d_stsh,&d_mllc,&b_stsh,&b_mllc);
  int stsh = d_stsh;
  int mllc = d_mllc;;
  if(blk_sz != 1)
  {
    stsh = b_stsh;
    mllc = b_mllc;
  }
  int sz = 0;
  MPI_Comm_size(cm,&sz);
  int * stshs = new int[sz];
  int * mllcs = new int[sz];
  MPI_Gather(&stsh,1,MPI_INTEGER,
             &stshs[0],1,MPI_INTEGER,
             0,cm);
  MPI_Gather(&mllc,1,MPI_INTEGER,
             &mllcs[0],1,MPI_INTEGER,
             0,cm);
  if(!PCU_Comm_Self())
  {
    std::cout << (blk_sz == 1 ? "double" : "block") << " stashes :\n";
    for(int rnk = 0; rnk < sz; ++rnk)
      std::cout << stshs[rnk] << " ";
    std::cout << "\n"
              << (blk_sz == 1 ? "double" : "block") << " mallocs :\n";
    for(int rnk = 0; rnk < sz; ++rnk)
      std::cout << mllcs[rnk] << " ";
    std::cout << std::endl;
  }
}
void m3dc1_matrix::printNNZStats()
{
  int rnk = -1;
  int sz = 0;
  MPI_Comm_rank(cm,&rnk);
  MPI_Comm_size(cm,&sz);
  MatInfo inf;
  MatInfoType tp = (sz == 1 ? MAT_LOCAL : MAT_GLOBAL_SUM);
  MatGetInfo(A,tp,&inf);
  if(rnk == 0)
  {
    std::cout << "Matrix nonzero info:\n"
              << "  Allocated: " << inf.nz_allocated << "\n"
              << "  Used     : " << inf.nz_used      << "\n"
              << "  Unneeded : " << inf.nz_unneeded  << "\n"
              << "  Mallocs  : " << inf.mallocs      << "\n"
              << "  Memory   : " << inf.memory       << std::endl;
  }
}
void m3dc1_matrix::zero()
{
  MatZeroEntries(A);
  fixed = false;
}
void m3dc1_matrix::zero_rows(int rw_cnt, int * rows)
{
  MatZeroRows(A,rw_cnt,rows,0.0,PETSC_NULL,PETSC_NULL);
  fixed = false;
}
int m3dc1_matrix::solver_iteration_count()
{
  int itr = -1;
  KSPGetIterationNumber(ksp,&itr);
  return itr;
}
#endif //#ifndef M3DC1_MESHGEN
