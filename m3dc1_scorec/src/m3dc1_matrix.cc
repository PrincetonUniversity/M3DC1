/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifdef M3DC1_PETSC
#include "m3dc1_matrix.h"
#include "apf.h"
#include "apfMesh.h"
#include "apfMDS.h"
#include <vector>
#include "PCU.h"
#include "m3dc1_mesh.h"
#include <assert.h>
#include <iostream>

#include "m3dc1_matrix_allocate.h"

#ifdef PETSC_USE_COMPLEX
#include "petscsys.h" // for PetscComplex
#include <complex>
using std::complex;
#endif

using std::vector;

// todo : account for complex type
void mat_insert_element_block(m3dc1_matrix * mat, m3dc1_mesh * msh, apf::MeshEntity * ent, double * vals)
{
  int fid = mat->get_fieldOrdering();
  m3dc1_field * fld = (*msh->field_container)[fid]; // todo: get rid of direct access to member variables
  // todo: assuming that only vertices hold nodes, should loop over dim and check in the field shape has nodes
  // build an element with the correct field (shape) and count the effecting nodes, use the apf functions to get the dofs
  //  far more straightforward
  apf::Element * elt = apf::createElement(fld->get_field(),ent);
  apf::Downward dn;
  msh->mesh->getDownward(ent,0,dn);
  int nds_per_elt = apf::countNodes(elt);
  int blks_per_nd = fld->get_num_value();
  int dofs_per_blk = fld->get_dof_per_value();
  int dofs_per_nd = dofs_per_blk * blks_per_nd;
  int blks_per_elt = blks_per_nd * nds_per_elt;
  // blk ids depend on whether the matrix is local or not
  int is_par = mat->is_parallel();
  int * mem = new int[blks_per_elt+2*dofs_per_nd];
  // DBG(memset(&mem[0],0,blks_per_elt+2*dofs_per_nd);
  int * blk_ids = &mem[0];
  int * lcl_dof_ids = &mem[blks_per_elt];
  // DBG(memset(&gbl_dof_ids[0],0,dofs_per_nd*sizeof(int));
  int * gbl_dof_ids = &mem[blks_per_elt+dofs_per_nd];
  int * dof_ids[] = {&lcl_dof_ids[0],&gbl_dof_ids[0]};
  // DBG(memset(&lcl_dof_ids[0],0,dofs_per_nd*sizeof(int));
  // assumes that the number of adjacent verts = number of elt nodes
  int dof_cnt = 0;
  for(int elt_nd = 0; elt_nd < nds_per_elt; ++elt_nd)
  {
    apf::MeshEntity * vrt = dn[elt_nd];
    //int ids[] = {get_ent_localid(msh->mesh,vrt), get_ent_globalid(msh->mesh,vrt)};
    get_ent_localdofid(fld,get_ent_localid(msh->mesh,vrt),dof_ids[0],&dof_cnt);
    get_ent_globaldofid(fld,get_ent_globalid(msh->mesh,vrt),dof_ids[1], &dof_cnt);
    for(int nd_blk = 0; nd_blk < blks_per_nd; ++nd_blk)
      blk_ids[elt_nd*blks_per_nd + nd_blk] = dof_ids[is_par][nd_blk*dofs_per_nd] / dofs_per_blk;
  }
  mat->add_blocks(blks_per_elt,blk_ids,blks_per_elt,blk_ids,vals);
  //dof_ids = NULL;
  //lcl_dof_ids = NULL;
  //gbl_dof_ids = NULL;
  //blk_ids = NULL;
  delete [] mem;
  apf::destroyElement(elt);
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
  std::cout<<"\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank "<<PCU_Comm_Self()<<" current "<<mem/1e6<<std::endl;
}

int copyField2PetscVec(FieldID field_id, Vec& petscVec, int scalar_type)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance(); // external var
  apf::Mesh2 * m = msh->mesh;
  int num_own_ent= msh->num_own_ent[0]; //assumes only verts have dofs
  int num_own_dof=0;
  int vertex_type=0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  int ierr = VecCreateMPI(MPI_COMM_WORLD, num_own_dof, PETSC_DECIDE, &petscVec);
  CHKERRQ(ierr);
  VecAssemblyBegin(petscVec);

  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0];

  double dof_data[FIXSIZEBUFF];
  assert(sizeof(dof_data)>=dofPerEnt*2*sizeof(double));
  int nodeCounter=0;

  apf::MeshEntity* ent;
  for (int inode=0; inode<num_vtx; inode++)
  {
    ent = getMdsEntity(m,vertex_type,inode);
    if (!is_ent_original(m,ent)) continue;
      nodeCounter++;
    int num_dof;
    m3dc1_ent_getdofdata (&vertex_type, &inode, &field_id, &num_dof, dof_data);
    assert(num_dof*(1+scalar_type)<=sizeof(dof_data)/sizeof(double));
    int start_global_dof_id, end_global_dof_id_plus_one;
    // FIXME: for blocked DOF's
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx=0;
    for (int i=0; i<dofPerEnt; i++)
    {
      PetscScalar value;
      if (scalar_type == M3DC1_REAL) value = dof_data[startIdx++];
      else
      {
#ifdef PETSC_USE_COMPLEX
        value = complex<double>(dof_data[startIdx*2],dof_data[startIdx*2+1]);
#else
        if (!PCU_Comm_Self())
          std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
        startIdx++;
      }
      ierr = VecSetValue(petscVec, start_global_dof_id+i, value, INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }
  assert(nodeCounter==num_own_ent);
  ierr=VecAssemblyEnd(petscVec);
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

void field2Vec(m3dc1_field * fld, Vec V, int st)
{
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // exernal var
  apf::Mesh2 * msh = m3dc1_msh->mesh;
  int num_own_nds = m3dc1_msh->num_own_ent[0]; // assuming only verts have dofs
  int num_own_dof = 0;
  FieldID fid = fld->get_id();
  m3dc1_field_getnumowndof(&fid, &num_own_dof);
  int dof_per_nd = num_own_dof / num_own_nds;
  VecCreateMPI(MPI_COMM_WORLD,num_own_dof,PETSC_DETERMINE,&V);
  int num_lcl_nds = m3dc1_msh->num_local_ent[0];
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
    if(!is_ent_original(msh,ent))
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

int copyPetscVec2Field(Vec& petscVec, FieldID field_id, int scalar_type)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0],num_own_dof=0, vertex_type=0;
  // FIXME: do not use API
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;

  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  std::vector<PetscInt> ix(dofPerEnt);
  std::vector<PetscScalar> values(dofPerEnt);
  std::vector<double> dof_data(dofPerEnt*(1+scalar_type));
  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0];

  int ierr;

  apf::MeshEntity* ent;
  for (int inode=0; inode<num_vtx; inode++)
  {
    ent = getMdsEntity(m, vertex_type,inode);
    if (!is_ent_original(m, ent)) continue;
    int start_global_dof_id, end_global_dof_id_plus_one;
    // FIXME: for blocked DOF's
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;

    for (int i=0; i<dofPerEnt; i++)
      ix.at(i)=startIdx+i;
    ierr=VecGetValues(petscVec, dofPerEnt, &ix[0], &values[0]); CHKERRQ(ierr);
    startIdx=0;
    for (int i=0; i<dofPerEnt; i++)
    {
      if (scalar_type == M3DC1_REAL)
      {
#ifdef PETSC_USE_COMPLEX
        dof_data.at(startIdx++)= values.at(i).real();
#else
        dof_data.at(startIdx++)= values.at(i);
#endif
      }
      else
      {
#ifdef PETSC_USE_COMPLEX
        dof_data.at(2*startIdx)=values.at(i).real();
        dof_data.at(2*startIdx+1)=values.at(i).imag();
        startIdx++;
#else
        if (!PCU_Comm_Self())
          std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, &field_id, &dofPerEnt, &dof_data[0]);
  }
  m3dc1_field_sync(&field_id);
  return M3DC1_SUCCESS;
}

void vec2Field(m3dc1_field * fld, Vec V, int st)
{
  m3dc1_mesh * m3dc1_msh = m3dc1_mesh::instance(); // external variable
  apf::Mesh2 * msh = m3dc1_msh->mesh;
  FieldID fid = fld->get_id();
  int num_own_nds = m3dc1_msh->num_own_ent[0];
  int num_own_dof = 0;
  m3dc1_field_getnumowndof(&fid,&num_own_dof);
  int vrt_tp = 0;
  int dof_per_nd = num_own_dof / num_own_nds;
  int num_lcl_nd = m3dc1_msh->num_local_ent[0];
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
    if(!is_ent_original(msh,ent))
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

// ***********************************
//              M3DC1_SOLVER
// ***********************************

m3dc1_solver* m3dc1_solver::_instance=NULL;
m3dc1_solver* m3dc1_solver::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_solver();
  return _instance;
}

m3dc1_solver::~m3dc1_solver()
{
  if (matrix_container!=NULL)
    matrix_container->clear();
  matrix_container=NULL;
  delete _instance;
}

void m3dc1_solver::add_matrix(int matrix_id, m3dc1_matrix* matrix)
{
  assert(matrix_container->find(matrix_id)==matrix_container->end());
  matrix_container->insert(std::map<int, m3dc1_matrix*>::value_type(matrix_id, matrix));
}

m3dc1_matrix* m3dc1_solver::get_matrix(int matrix_id)
{
  std::map<int, m3dc1_matrix*>::iterator mit = matrix_container->find(matrix_id);
  if (mit == matrix_container->end())
    return (m3dc1_matrix*)NULL;
  return mit->second;
}

// ***********************************
//              M3DC1_MATRIX
// ***********************************

m3dc1_matrix::m3dc1_matrix(int i, int s, FieldID f)
  : id(i)
  , scalar_type(s)
  , fieldOrdering(f)
{
  mat_status = M3DC1_NOT_FIXED;
}

m3dc1_matrix::~m3dc1_matrix()
{
  MatDestroy(&A);
}

inline void m3dc1_matrix::add_blocks(int blk_rw_cnt, int * blk_rws, int blk_col_cnt, int * blk_cols, double * vals)
{
  MatSetValuesBlocked(A,blk_rw_cnt,blk_rws,blk_col_cnt,blk_cols,vals,ADD_VALUES);
}

int m3dc1_matrix::get_values(vector<int>& rows, vector<int>& n_columns, vector<int>& columns, vector<double>& values)
{
  if (mat_status != M3DC1_FIXED)
    return M3DC1_FAILURE;
#ifdef PETSC_USE_COMPLEX
   if (!PCU_Comm_Self())
     std::cout<<"[M3DC1 ERROR] "<<__func__<<": not supported for complex\n";
   return M3DC1_FAILURE;
#else
  PetscErrorCode ierr;
  PetscInt rstart, rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  ierr = MatGetOwnershipRange(A, &rstart, &rend);
  CHKERRQ(ierr);
  for (PetscInt row=rstart; row<rend; ++row)
  {
    ierr = MatGetRow(A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
    rows.push_back(row);
    n_columns.push_back(ncols);
    for (int i=0; i<ncols; ++i)
    {
      columns.push_back(cols[i]);
      values.push_back(vals[i]);
    }
    ierr = MatRestoreRow(A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
  }
  assert(rows.size()==rend-rstart);
#endif
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::set_value(int row, int col, int operation, double real_val, double imag_val) //insertion/addition with global numbering
{
  if (mat_status == M3DC1_FIXED)
    return M3DC1_FAILURE;
  PetscErrorCode ierr;

  if (scalar_type==M3DC1_REAL) // real
  {
    if (operation)
      ierr = MatSetValue(A, row, col, real_val, ADD_VALUES);
    else
      ierr = MatSetValue(A, row, col, real_val, INSERT_VALUES);
  }
  else // complex
  {
#ifdef PETSC_USE_COMPLEX
    PetscScalar value = complex<double>(real_val,imag_val);
    if (operation)
      ierr = MatSetValue(A, row, col, value, ADD_VALUES);
    else
      ierr = MatSetValue(A, row, col, value, INSERT_VALUES);
#else
    if (!PCU_Comm_Self())
      std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
      abort();
#endif
  }
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::add_values(int rsize, int * rows, int csize, int * columns, double* values)
{
  if (mat_status == M3DC1_FIXED)
    return M3DC1_FAILURE;
  PetscErrorCode ierr;
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  vector<PetscScalar> petscValues(rsize*csize);
  for (int i=0; i<rsize; i++)
  {
    //if (id==22)
      //std::cout<<std::endl<<"id "<<id<<" row "<<rows[i]<<std::endl;
    for (int j=0; j<csize; j++)
    {
      //if (id==22)
        //std::cout<<" colum "<<columns[j]<<" "<<values[i*csize+j]<<" ";
      if (scalar_type==M3DC1_REAL) petscValues.at(i*csize+j)=values[i*csize+j];
      else
      {
#ifdef PETSC_USE_COMPLEX
        petscValues.at(i*csize+j)=complex<double>(values[2*i*csize+2*j], values[2*i*csize+2*j+1]);
#else
        if (!PCU_Comm_Self())
        std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
  }

  ierr = MatSetValues(A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
#else
  ierr = MatSetValues(A, rsize, rows, csize, columns, (PetscScalar*)values, ADD_VALUES);
#endif
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int matrix_solve::setUpRemoteAStruct()
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int dofPerVar = 6, vertex_type=0;
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];

  std::vector<int> nnz_remote(num_values*num_vtx);
  int brgType = m->getDimension();

  apf::MeshEntity* ent;
  for (int inode=0; inode<num_vtx; inode++)
  {
    ent = getMdsEntity(m, vertex_type, inode);
    int owner=get_ent_ownpartid(m, ent);
    if (owner!=PCU_Comm_Self())
    {
      apf::Adjacent elements;
      getBridgeAdjacent(m, ent, brgType, 0, elements);
      int num_elem=0;
      for (int i=0; i<elements.getSize(); ++i)
      {
        if (!m->isGhost(elements[i]))
          ++num_elem;
      }

      remoteNodeRow[owner][inode]=num_elem+1;
      remoteNodeRowSize[owner]+=num_elem+1;
      for (int i=0; i<num_values; i++)
        nnz_remote[inode*num_values+i]=(num_elem+1)*num_values;
    }
    else
    {
      apf::Copies remotes;
      m->getRemotes(ent,remotes);
      APF_ITERATE(apf::Copies, remotes, it)
        remotePidOwned.insert(it->first);
    }
  }
  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF,&remoteA);
  CHKERRQ(ierr);
  ierr = MatSetType(remoteA, MATSEQBAIJ);CHKERRQ(ierr);
  ierr = MatSetBlockSize(remoteA, dofPerVar); CHKERRQ(ierr);
  ierr = MatSetSizes(remoteA, total_num_dof*num_vtx, total_num_dof*num_vtx, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  MatSeqBAIJSetPreallocation(remoteA, dofPerVar, 0, &nnz_remote[0]);
  ierr = MatSetUp (remoteA);CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::write (const char* file_name)
{
  PetscErrorCode ierr;
  PetscViewer lab;
  if (get_type()==0)
  {
    char name_buff[256];
    sprintf(name_buff, "%s-%d.m",file_name,PCU_Comm_Self());
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name_buff, &lab); CHKERRQ(ierr);
  }
  else
  {
    ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, file_name, &lab); CHKERRQ(ierr);
  }
  ierr = PetscViewerPushFormat(lab, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
  ierr = MatView(A, lab); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&lab); CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::printInfo()
{
  MatInfo info;
  MatGetInfo(A, MAT_LOCAL,&info);
  std::cout<<"Matrix "<<id<<" info "<<std::endl;
  std::cout<<"\t nz_allocated,nz_used,nz_unneeded "<<info.nz_allocated<<" "<<info.nz_used<<" "<<info.nz_unneeded<<std::endl;
  std::cout<<"\t memory mallocs "<<info.memory<<" "<<info.mallocs<<std::endl;
  PetscInt nstash, reallocs, bnstash, breallocs;
  MatStashGetInfo(A,&nstash,&reallocs,&bnstash,&breallocs);
  std::cout<<"\t nstash, reallocs, bnstash, breallocs "<<nstash<<" "<<reallocs<<" "<<bnstash<<" "<<breallocs<<std::endl;
  return M3DC1_SUCCESS;
}
// ***********************************
//              MATRIX_MULTIPLY
// ***********************************
matrix_mult::matrix_mult(int i , int s, FieldID fld)
  : m3dc1_matrix(i,s,fld)
{
  is_par = 0;
  m3dc1_mesh * msh = m3dc1_mesh::instance(); //external variable, pass in or use field
  int num_ent = msh->num_local_ent[0]; // assumes that only verts hold dofs
  m3dc1_field * mf = (*msh->field_container)[fld];
  int blk_sz = mf->get_dof_per_value();
  int dof_per_ent = mf->get_num_value() * blk_sz;
  int num_lcl_dof = num_ent * dof_per_ent;
  MatCreate(PETSC_COMM_SELF,&A);
  const char * mat_tp = (blk_sz == 1 ? MATSEQAIJ : MATSEQBAIJ);
  MatSetType(A,mat_tp);
  MatSetSizes(A,num_lcl_dof,num_lcl_dof,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetBlockSize(A,blk_sz);
  // call preallocate
  allocateMatrix(A,msh,mf);
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  if(!(blk_sz-1)) // only supported for AIJ not BAIJ
    MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
#ifdef DEBUG
  if(!PCU_Comm_Self())
    describeMatrix(A);
#endif
}

int matrix_mult::assemble()
{
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  set_status(M3DC1_FIXED);
  return M3DC1_SUCCESS;
}

int matrix_mult::multiply(FieldID in_field, FieldID out_field)
{
  if (!localMat)
  {
    Vec b, c;
    copyField2PetscVec(in_field, b, get_scalar_type());
    int ierr = VecDuplicate(b, &c);CHKERRQ(ierr);
    MatMult(A, b, c);
    copyPetscVec2Field(c, out_field, get_scalar_type());
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
    return 0;
  }
  else
  {
    Vec b, c;
    m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[in_field];
    int num_dof = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();

#ifdef DEBUG
    m3dc1_field * mf2 = (*(m3dc1_mesh::instance()->field_container))[out_field];
    int num_dof2 = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
    assert(num_dof==num_dof2);
#endif
    int bs;
    int ierr;
    MatGetBlockSize(A, &bs);
    PetscScalar * array[2];
    // FIXME: for blocked DOF's
    int vid=0;
    m3dc1_field_getdataptr(&in_field, (double**)array);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double * array_org = (double*)array[0];
      array[0] = new PetscScalar[num_dof];
      for (int i=0; i<num_dof; i++)
      {
        array[0][i]=array_org[i];
      }
    }
#endif
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[0],&b); CHKERRQ(ierr);
    // FIXME: for blocked DOF's
    m3dc1_field_getdataptr(&out_field, (double**)array+1);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double * array_org = (double*)array[1];
      array[1] = new PetscScalar[num_dof];
      for (int i=0; i<num_dof; i++)
      {
        array[1][i]=array_org[i];
      }
    }
#endif
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[1],&c); CHKERRQ(ierr);
    ierr=VecAssemblyBegin(b);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(b);  CHKERRQ(ierr);
    ierr=VecAssemblyBegin(c);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(c);  CHKERRQ(ierr);
    MatMult(A, b, c);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double *datapt;
      // FIXME: for blocked DOF's
      m3dc1_field_getdataptr(&out_field, &vid, &datapt);
      for (int i=0; i<num_dof; i++)
        datapt[i]=std::real(array[1][i]);
      delete []array[0];
      delete []array[1];
    }
#endif
    m3dc1_field_sum(&out_field);
  }
  return M3DC1_SUCCESS;
}

// ***********************************
//              MATRIX_SOLVE
// ***********************************

matrix_solve::matrix_solve(int i, int s, FieldID fld)
  : m3dc1_matrix(i,s,fld)
  , kspSet(0)
{
  is_par = 1;
  m3dc1_mesh * msh = m3dc1_mesh::instance(); // external variable
  int num_own_nds = msh->num_own_ent[0]; // assumes only vertices hold dofs
  m3dc1_field * mf = (*msh->field_container)[fld];
  int blk_sz = mf->get_dof_per_value();
  int dof_per_nd = mf->get_num_value() * blk_sz;
  int num_own_dof = num_own_nds * dof_per_nd;
  MatCreate(PETSC_COMM_WORLD,&A);
  const char * par_mat_tp = (blk_sz == 1 ? MATMPIAIJ : MATMPIBAIJ);
  MatSetType(A,par_mat_tp);
  MatSetSizes(A,num_own_dof,num_own_dof,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetBlockSize(A,blk_sz);
  allocateMatrix(A,msh,mf);
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  if(!(blk_sz-1)) // only supported for AIJ not BAIJ
    MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
  // TODO: remove remoteA
  int num_lcl_nds = msh->num_local_ent[0];
  int num_lcl_dof = num_lcl_nds * dof_per_nd;
  MatCreate(PETSC_COMM_SELF,&remoteA);
  const char * seq_mat_tp = (blk_sz == 1 ? MATSEQAIJ : MATSEQBAIJ);
  MatSetType(remoteA,seq_mat_tp);
  MatSetSizes(remoteA,num_lcl_dof,num_lcl_dof,PETSC_DETERMINE,PETSC_DETERMINE);
  MatSetBlockSize(remoteA,blk_sz);
  setUpRemoteAStruct();
  MatSetOption(remoteA,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  if(!(blk_sz-1)) // only supported for AIJ not BAIJ
    MatSetOption(remoteA,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
#ifdef DEBUG
  if(!PCU_Comm_Self())
    describeMatrix(A);
#endif
}

matrix_solve::~matrix_solve()
{
  MatDestroy(&remoteA);
  MatDestroy(&A);
}

void matrix_solve::reset_values()
{
  MatZeroEntries(A);
  MatZeroEntries(remoteA);
  set_status(M3DC1_NOT_FIXED); // allow matrix value modification
#ifdef DEBUG_
  PetscInt rstart, rend, r_rstart, r_rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  //MatGetSize(A, &n, NULL); -- this returns global matrix size
  MatGetOwnershipRange(A, &rstart, &rend);
  MatGetOwnershipRange(remoteA, &r_rstart, &r_rend);

  for (PetscInt row=rstart; row<rend; ++row)
  {
    MatGetRow(A, row, &ncols, &cols, &vals);
    for (int i=0; i<ncols; ++i)
      assert(m3dc1_double_isequal(vals[i],0.0));
    MatRestoreRow(A, row, &ncols, &cols, &vals); // prevent memory leak
  }

  for (PetscInt row=r_rstart; row<r_rend; ++row)
  {
    MatGetRow(remoteA, row, &ncols, &cols, &vals);
    for (int i=0; i<ncols; ++i)
      assert(m3dc1_double_isequal(vals[i],0.0));
    MatRestoreRow(remoteA, row, &ncols, &cols, &vals); // prevent memory leak
  }
#endif
};


int matrix_solve::add_blockvalues(int rbsize, int * rows, int cbsize, int * columns, double* values)
{
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  int bs;
  MatGetBlockSize(remoteA, &bs);
  vector<PetscScalar> petscValues(rbsize*cbsize*bs*bs);

  for (int i=0; i<rbsize*bs; i++)
  {
    for (int j=0; j<cbsize*bs; j++)
    {
      if (scalar_type==M3DC1_REAL) petscValues.at(i*cbsize*bs+j)=values[i*cbsize*bs+j];
      else
      {
#ifdef PETSC_USE_COMPLEX
        petscValues.at(i*cbsize*bs+j)=complex<double>(values[2*i*cbsize*bs+2*j], values[2*i*cbsize*bs+2*j+1]);
#else
        if (!PCU_Comm_Self())
        std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
  }
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, &petscValues[0], ADD_VALUES);
#else
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, (PetscScalar*)values, ADD_VALUES);
#endif
  return M3DC1_SUCCESS;
}

int matrix_solve::assemble()
{
  PetscErrorCode ierr;
  double t1 = MPI_Wtime(), t2=t1;
  if (!m3dc1_solver::instance()->assembleOption)
  {
    ierr = MatAssemblyBegin(remoteA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(remoteA, MAT_FINAL_ASSEMBLY);
    t2 = MPI_Wtime();
    //pass remoteA to ownnering process
    int brgType = m3dc1_mesh::instance()->mesh->getDimension();

    int dofPerVar = 6;
    char field_name[256];
    int num_values, value_type, total_num_dof, vertex_type=0;
    m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
    dofPerVar=total_num_dof/num_values;

    int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];
    PetscInt firstRow, lastRowPlusOne;
    ierr = MatGetOwnershipRange(A, &firstRow, &lastRowPlusOne);

    std::map<int, std::vector<int> > idxSendBuff, idxRecvBuff;
    std::map<int, std::vector<PetscScalar> > valuesSendBuff, valuesRecvBuff;
    int blockMatSize = total_num_dof*total_num_dof;
    for (std::map<int, std::map<int, int> > ::iterator it = remoteNodeRow.begin(); it!=remoteNodeRow.end(); it++)
    {
      idxSendBuff[it->first].resize(it->second.size()+remoteNodeRowSize[it->first]);
      valuesSendBuff[it->first].resize(remoteNodeRowSize[it->first]*blockMatSize);
      int idxOffset=0;
      int valueOffset=0;
      for (std::map<int, int> ::iterator it2 =it->second.begin(); it2!=it->second.end();it2++)
      {
        idxSendBuff[it->first].at(idxOffset++)=it2->second;
        apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, it2->first);

        std::vector<apf::MeshEntity*> vecAdj;
        apf::Adjacent elements;
        getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);
        for (int i=0; i<elements.getSize(); ++i)
        {
          if (!m3dc1_mesh::instance()->mesh->isGhost(elements[i]))
            vecAdj.push_back(elements[i]);
        }
        vecAdj.push_back(ent);
        int numAdj = vecAdj.size();
        assert(numAdj==it2->second);
        std::vector<int> localNodeId(numAdj);
        std::vector<int> columns(total_num_dof*numAdj);
        for (int i=0; i<numAdj; i++)
        {
          int local_id = get_ent_localid(m3dc1_mesh::instance()->mesh, vecAdj.at(i));
          localNodeId.at(i)=local_id;
          int start_global_dof_id, end_global_dof_id_plus_one;
          // FIXME: for blocked DOF's
          m3dc1_ent_getglobaldofid (&vertex_type, &local_id, &fieldOrdering, &start_global_dof_id, &end_global_dof_id_plus_one);
          idxSendBuff[it->first].at(idxOffset++)=start_global_dof_id;
        }
        int offset=0;
        for (int i=0; i<numAdj; i++)
        {
          int startColumn = localNodeId.at(i)*total_num_dof;
          for (int j=0; j<total_num_dof; j++)
            columns.at(offset++)=startColumn+j;
        }
        ierr = MatGetValues(remoteA, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &valuesSendBuff[it->first].at(valueOffset));
        valueOffset+=it2->second*blockMatSize;
      }
      assert(idxOffset==idxSendBuff[it->first].size());
      assert(valueOffset==valuesSendBuff[it->first].size());
    }
    // ierr = MatDestroy(&remoteA); // seol: shall destroy in destructor

    //send and receive message size
    int sendTag=2020;
    MPI_Request my_request[256];
    MPI_Status my_status[256];
    int requestOffset=0;
    std::map<int, std::pair<int, int> > msgSendSize;
    std::map<int, std::pair<int, int> > msgRecvSize;
    for (std::map<int, int >::iterator it = remoteNodeRowSize.begin(); it!=remoteNodeRowSize.end(); it++)
    {
      int destPid=it->first;
      msgSendSize[destPid].first=idxSendBuff[it->first].size();
      msgSendSize[destPid].second = valuesSendBuff[it->first].size();
      MPI_Isend(&(msgSendSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for (std::set<int>::iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      MPI_Irecv(&(msgRecvSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request,my_status);
    //set up receive buff
    for (std::map<int, std::pair<int, int> >::iterator it = msgRecvSize.begin(); it!= msgRecvSize.end(); it++)
    {
      idxRecvBuff[it->first].resize(it->second.first);
      valuesRecvBuff[it->first].resize(it->second.second);
    }
    // now get data
    sendTag=9999;
    requestOffset=0;
    for (std::map<int, int >::iterator it = remoteNodeRowSize. begin(); it!=remoteNodeRowSize.end(); it++)
    {
      int destPid=it->first;
      MPI_Isend(&(idxSendBuff[destPid].at(0)),idxSendBuff[destPid].size(),MPI_INT,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Isend(&(valuesSendBuff[destPid].at(0)),sizeof(PetscScalar)*valuesSendBuff[destPid].size(),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for (std::set<int>::iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      MPI_Irecv(&(idxRecvBuff[destPid].at(0)),idxRecvBuff[destPid].size(),MPI_INT,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Irecv(&(valuesRecvBuff[destPid].at(0)),sizeof(PetscScalar)*valuesRecvBuff[destPid].size(),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request,my_status);

    for ( std::map<int, std::vector<int> >::iterator it =idxSendBuff.begin(); it!=idxSendBuff.end(); it++)
      std::vector<int>().swap(it->second);
    for ( std::map<int, std::vector<PetscScalar> >::iterator it =valuesSendBuff.begin(); it!=valuesSendBuff.end(); it++)
      std::vector<PetscScalar>().swap(it->second);
    valuesSendBuff.clear();
    idxSendBuff.clear();

    // now assemble the matrix
    for (std::set<int>::iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      int valueOffset=0;
      int idxOffset=0;
      vector<int> & idx = idxRecvBuff[destPid];
      vector<PetscScalar> & values = valuesRecvBuff[destPid];
      int numValues=values.size();
      while (valueOffset<numValues)
      {
        int numAdj = idx.at(idxOffset++);
        std::vector<int> columns(total_num_dof*numAdj);
        int offset=0;
        for (int i=0; i<numAdj; i++, idxOffset++)
        {
          for (int j=0; j<total_num_dof; j++)
          {
            columns.at(offset++)=idx.at(idxOffset)+j;
          }
        }
        assert (columns.at(total_num_dof*(numAdj-1))>=firstRow && *columns.rbegin()<lastRowPlusOne);
        ierr = MatSetValues(A, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &values.at(valueOffset),ADD_VALUES);

        valueOffset+=blockMatSize*numAdj;
      }
      std::vector<int>().swap(idxRecvBuff[destPid]);
      std::vector<PetscScalar>().swap(valuesRecvBuff[destPid]);
    }
    valuesRecvBuff.clear();
    idxRecvBuff.clear();
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  mat_status=M3DC1_FIXED;
  return M3DC1_SUCCESS;
}

int matrix_solve:: set_bc( int row)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  MatSetValue(A, row, row, 1.0, ADD_VALUES);
  return M3DC1_SUCCESS;
}

int matrix_solve:: set_row( int row, int numVals, int* columns, double * vals)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  for (int i=0; i<numVals; i++)
  {
    if (get_scalar_type() == M3DC1_REAL) set_value(row, columns[i], 1, vals[i], 0);
    else set_value(row, columns[i], 1, vals[2*i], vals[2*i+1]);
  }
  return M3DC1_SUCCESS;
}

int matrix_solve::solve(FieldID field_id)
{
  Vec x, b;
  copyField2PetscVec(field_id, b, get_scalar_type());
  int ierr = VecDuplicate(b, &x);CHKERRQ(ierr);
  ksp = new KSP;
  setKspType();
  ierr = KSPSolve(*ksp, b, x);
  CHKERRQ(ierr);
  PetscInt its;
  ierr = KSPGetIterationNumber(*ksp, &its);
  CHKERRQ(ierr);
  int iter_num=its;
  if (PCU_Comm_Self() == 0)
    std::cout <<"\t-- # solver iterations " << its << std::endl;
  iterNum = its;
  //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  copyPetscVec2Field(x, field_id, get_scalar_type());
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  // delete ksp
  KSPDestroy(ksp);
  delete ksp;
  return M3DC1_SUCCESS;
}

int matrix_solve:: setKspType()
{
  PetscErrorCode ierr;
  KSPCreate(MPI_COMM_WORLD, ksp);
  ierr = KSPSetOperators(*ksp, A, A /*, SAME_PRECONDITIONER DIFFERENT_NONZERO_PATTERN*/);CHKERRQ(ierr);
  ierr = KSPSetTolerances(*ksp, .000001, .000000001,
                          PETSC_DEFAULT, 1000);CHKERRQ(ierr);
  int num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  assert(total_num_dof/num_values==C1TRIDOFNODE*(m3dc1_mesh::instance()->mesh->getDimension()-1));
  // if 2D problem use superlu
  if (m3dc1_mesh::instance()->mesh->getDimension()==2)
  {
    ierr=KSPSetType(*ksp, KSPPREONLY);CHKERRQ(ierr);
    PC pc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr=PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);  CHKERRQ(ierr);
  }

  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
 // ++kspSet;
  return M3DC1_SUCCESS;
}

#endif //#ifndef M3DC1_MESHGEN
