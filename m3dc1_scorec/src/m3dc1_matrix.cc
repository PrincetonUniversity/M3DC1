/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifdef M3DC1_PETSC
#include "m3dc1_matrix.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh.h"
#include <vector>
#include "PCU.h"
#include "m3dc1_mesh.h"
#include <assert.h>
#include <iostream>

#ifdef PETSC_USE_COMPLEX
#include "petscsys.h" // for PetscComplex
#include <complex>
using std::complex;
#endif

using std::vector;

// ***********************************
// 		HELPER
// ***********************************

void printMemStat()
{
  PetscLogDouble mem, mem_max;
  PetscMemoryGetCurrentUsage(&mem);
  PetscMemoryGetMaximumUsage(&mem_max);
  std::cout<<"\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank "<<PCU_Comm_Self()<<" current "<<mem/1e6<<std::endl;
}

int copyField2PetscVec(FieldID field_id, Vec& petscVec, int scalar_type)
{
  int num_own_ent= m3dc1_mesh::instance()->num_own_ent[0];
  int num_own_dof=0, vertex_type=0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

/*int ierr = VecCreateMPI(MPI_COMM_WORLD, num_own_dof, PETSC_DECIDE, &petscVec); */
  int ierr = VecCreate(MPI_COMM_WORLD, &petscVec); CHKERRQ(ierr);
  ierr = VecSetSizes(petscVec, num_own_dof, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = VecSetFromOptions(petscVec);CHKERRQ(ierr);
  VecAssemblyBegin(petscVec);

  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0];

  double dof_data[FIXSIZEBUFF];
  assert(sizeof(dof_data)>=dofPerEnt*2*sizeof(double));
  int nodeCounter=0;

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* ent;
  int inode;
  apf::MeshIterator* ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it)))
  {
    inode = getMdsIndex(mesh, ent);
    if (!is_ent_original(mesh,ent)) continue;
    nodeCounter+=1;
    int num_dof;
    m3dc1_ent_getdofdata (&vertex_type, &inode, &field_id, &num_dof, dof_data);
    assert(num_dof*(1+scalar_type)<=sizeof(dof_data)/sizeof(double));
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx=0;
    for (int i=0; i<dofPerEnt; ++i)
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
        startIdx+=1;
      } 
      ierr = VecSetValue(petscVec, start_global_dof_id+i, value, INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }
  mesh->end(ent_it);

  assert(nodeCounter==num_own_ent);
  ierr=VecAssemblyEnd(petscVec);
  CHKERRQ(ierr);
  return 0;
}

int copyPetscVec2Field(Vec& petscVec, FieldID field_id, int scalar_type)
{
  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0], num_own_dof=0, vertex_type=0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  std::vector<PetscInt> ix(dofPerEnt);
  std::vector<PetscScalar> values(dofPerEnt);
  std::vector<double> dof_data(dofPerEnt*(1+scalar_type));
  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0];

  int ierr, inode;

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* ent;
  apf::MeshIterator* ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it)))
  {
    inode = getMdsIndex(mesh, ent);
    if (!is_ent_original(mesh, ent)) continue;
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    
    for (int i=0; i<dofPerEnt; ++i)
      ix.at(i)=startIdx+i;
    ierr=VecGetValues(petscVec, dofPerEnt, &ix[0], &values[0]); CHKERRQ(ierr);
    startIdx=0;
    for (int i=0; i<dofPerEnt; ++i)
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
        ++startIdx;
#else 
        if (!PCU_Comm_Self())
          std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, &field_id, &dofPerEnt, &dof_data[0]);
  }
  mesh->end(ent_it);

  synchronize_field((*m3dc1_mesh::instance()->field_container)[field_id]->get_field());
  return 0;
}

// ***********************************
// 		M3DC1_SOLVER
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
  std::map<int, m3dc1_matrix*>::iterator mit=matrix_container->find(matrix_id);
  if (mit == matrix_container->end()) 
    return (m3dc1_matrix*)NULL;
  return mit->second;
}

// ***********************************
// 		M3DC1_MATRIX
// ***********************************

m3dc1_matrix::m3dc1_matrix(int i, int s, FieldID f)
: mesh(m3dc1_mesh::instance()->mesh), id(i), scalar_type(s), fieldOrdering(f)
{
  mat_status = M3DC1_NOT_FIXED;
  A=new Mat;
}

int m3dc1_matrix::destroy()
{
  PetscErrorCode ierr = MatDestroy(A);
  CHKERRQ(ierr);    
  return M3DC1_SUCCESS;
}

m3dc1_matrix::~m3dc1_matrix()
{
  destroy();
  delete A;
} 

int m3dc1_matrix::get_values(vector<int>& rows, vector<int>& n_columns, vector<int>& columns, vector<double>& values)
{
  if (!mat_status)  // matrix is not fixed
  {
    if (!PCU_Comm_Self())
      std::cout <<__func__<<" failed: matrix "<<id<<" is not fixed\n";
    return M3DC1_FAILURE;
  }

#ifdef PETSC_USE_COMPLEX
   if (!PCU_Comm_Self())
     std::cout<<"[M3DC1 ERROR] "<<__func__<<": not supported for complex\n";
   return M3DC1_FAILURE;
#else
  PetscErrorCode ierr;
  PetscInt rstart, rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  ierr = MatGetOwnershipRange(*A, &rstart, &rend);
  CHKERRQ(ierr);
  for (PetscInt row=rstart; row<rend; ++row)
  { 
    ierr = MatGetRow(*A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
    rows.push_back(row);
    n_columns.push_back(ncols);
    for (int i=0; i<ncols; ++i)
    {  
      columns.push_back(cols[i]);
      values.push_back(vals[i]);
    }
    ierr = MatRestoreRow(*A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
  }
  assert(rows.size()==rend-rstart);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_matrix::set_value(int row, int col, int operation, double real_val, double imag_val) //insertion/addition with global numbering
{
  if (mat_status) // matrix is fixed
  {
    if (!PCU_Comm_Self())
      std::cout <<__func__<<" failed: matrix "<<id<<" is fixed\n";
    return M3DC1_FAILURE;
  }

  PetscErrorCode ierr;
  
  if (scalar_type==M3DC1_REAL) // real
  {
    if (operation)
      ierr = MatSetValue(*A, row, col, real_val, ADD_VALUES);
    else
      ierr = MatSetValue(*A, row, col, real_val, INSERT_VALUES);
  }
  else // complex
  {
#ifdef PETSC_USE_COMPLEX
    PetscScalar value = complex<double>(real_val,imag_val);
    if (operation)
      ierr = MatSetValue(*A, row, col, value, ADD_VALUES);
    else
      ierr = MatSetValue(*A, row, col, value, INSERT_VALUES);
#else
    if (!PCU_Comm_Self())
      std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is not configured with --with-scalar-type=complex\n";
      abort();
#endif
  }
  CHKERRQ(ierr);
}

int m3dc1_matrix::add_values(int rsize, int * rows, int csize, int * columns, double* values)
{
  if (mat_status) // matrix is fixed
  {
    if (!PCU_Comm_Self())
      std::cout <<__func__<<" failed: matrix "<<id<<" is fixed\n";
    return M3DC1_FAILURE;
  }

  PetscErrorCode ierr;
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  vector<PetscScalar> petscValues(rsize*csize);
  for (int i=0; i<rsize; ++i)
  {
    for (int j=0; j<csize; ++j)
    {
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

  ierr = MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
#else
  ierr = MatSetValues(*A, rsize, rows, csize, columns, (PetscScalar*)values, ADD_VALUES);
#endif
  CHKERRQ(ierr);
}

int matrix_mult::setupMat()
{
  if (localMat) setupSeqMat();
  else setupParaMat();
}

int matrix_mult::preAllocate ()
{
  if (localMat) preAllocateSeqMat();
  else preAllocateParaMat();
}


int  m3dc1_matrix::preAllocateParaMat()
{
  int bs=1;
  MatType type;
  MatGetType(*A, &type);

  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0],num_own_dof=0, vertex_type=0;
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  if (strcmp(type, MATSEQAIJ)==0 || strcmp(type, MATMPIAIJ)==0) 
    bs=1;
  else 
    bs=dofPerEnt;
  int numBlocks = num_own_dof / bs;
  int numBlockNode = dofPerEnt / bs;
  std::vector<PetscInt> dnnz(numBlocks), onnz(numBlocks);
  int startDof, endDofPlusOne;
  m3dc1_field_getowndofid (&fieldOrdering, &startDof, &endDofPlusOne);

  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0], inode;

  int nnzStash=0;
  int brgType = mesh->getDimension();

  apf::MeshEntity* ent;
  apf::MeshIterator* ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it)))
  // for (int inode=0; inode<num_vtx; ++inode)
  {
    inode=getMdsIndex(mesh, ent); //  ent = getMdsEntity(mesh, vertex_type, inode);
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &fieldOrdering, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    if (start_global_dof_id<startDof || start_global_dof_id>=endDofPlusOne)
    {
      apf::Adjacent elements;
      getBridgeAdjacent(mesh, ent, brgType, 0, elements);
      int num_elem=0;
      for (int i=0; i<elements.getSize(); ++i)
      {
        if (!mesh->isGhost(elements[i]))
          ++num_elem;
      }

      nnzStash+=dofPerEnt*dofPerEnt*(num_elem+1);
      continue;
    }
    startIdx -= startDof;
    startIdx /=bs; 

    int adjNodeOwned, adjNodeGlb;
    mesh->getIntTag(ent, m3dc1_mesh::instance()->num_global_adj_node_tag, &adjNodeGlb);
    mesh->getIntTag(ent, m3dc1_mesh::instance()->num_own_adj_node_tag, &adjNodeOwned);
    assert(adjNodeGlb>=adjNodeOwned);

    for (int i=0; i<numBlockNode; ++i)
    {
      dnnz.at(startIdx+i)=(1+adjNodeOwned)*numBlockNode;
      onnz.at(startIdx+i)=(adjNodeGlb-adjNodeOwned)*numBlockNode;
    }
  }
  mesh->end(ent_it);

  if (bs==1) 
    MatMPIAIJSetPreallocation(*A, 0, &dnnz[0], 0, &onnz[0]);
  else  
    MatMPIBAIJSetPreallocation(*A, bs, 0, &dnnz[0], 0, &onnz[0]);
} 

int matrix_solve::setUpRemoteAStruct()
{
  assert(remotePidOwned==NULL && remoteNodeRow==NULL && remoteNodeRowSize==NULL);

  remotePidOwned = new std::set<int>;
  remoteNodeRow = new std::map<int, std::map<int, int> >;
  remoteNodeRowSize = new std::map<int, int>;

  int dofPerVar = 6, vertex_type=0;
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];

  std::vector<int> nnz_remote(num_values*num_vtx);
  int brgType = mesh->getDimension();
  
  apf::MeshEntity* ent;
  apf::MeshIterator* ent_it = mesh->begin(0);
  int inode;
  while ((ent = mesh->iterate(ent_it)))
  {
    inode = getMdsIndex(mesh, ent);
    int owner=get_ent_ownpartid(mesh, ent);
    if (owner!=PCU_Comm_Self())
    {
      apf::Adjacent elements;
      getBridgeAdjacent(mesh, ent, brgType, 0, elements);
      int num_elem=0;
      for (int i=0; i<elements.getSize(); ++i)
      {
        if (!mesh->isGhost(elements[i]))
          ++num_elem;
      }

      (*remoteNodeRow)[owner][inode]=num_elem+1;
      (*remoteNodeRowSize)[owner]+=num_elem+1;
      for (int i=0; i<num_values; ++i)
        nnz_remote[inode*num_values+i]=(num_elem+1)*num_values;
    }
    else 
    {
      apf::Copies remotes;
      mesh->getRemotes(ent,remotes);
      APF_ITERATE(apf::Copies, remotes, it)
        remotePidOwned->insert(it->first);
    }
  }
  mesh->end(ent_it);

  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF,&remoteA);
  CHKERRQ(ierr);
  ierr = MatSetType(remoteA, MATSEQBAIJ); CHKERRQ(ierr);
  ierr = MatSetBlockSize(remoteA, dofPerVar); CHKERRQ(ierr);
  ierr = MatSetSizes(remoteA, total_num_dof*num_vtx, total_num_dof*num_vtx, PETSC_DECIDE, PETSC_DECIDE); 
  CHKERRQ(ierr);
  MatSeqBAIJSetPreallocation(remoteA, dofPerVar, 0, &nnz_remote[0]);
  ierr = MatSetUp (remoteA); CHKERRQ(ierr);
}

int  m3dc1_matrix::preAllocateSeqMat()
{
  int bs=1, vertex_type=0;
  MatType type;
  MatGetType(*A, &type);

  int num_vtx=m3dc1_mesh::instance()->num_local_ent[0];
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[fieldOrdering];
  int num_dof = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();

  int dofPerEnt=0;
  if (num_vtx) dofPerEnt = num_dof/num_vtx;

  if (strcmp(type, MATSEQAIJ)==0 || strcmp(type, MATMPIAIJ)==0) 
    bs=1;
  else 
    bs=dofPerEnt;
  int numBlocks = num_dof / bs;
  int numBlockNode = dofPerEnt / bs;
  std::vector<PetscInt> nnz(numBlocks);
  int brgType = 2;
  if (mesh->getDimension()==3) brgType = 3;

  apf::MeshEntity* ent;
  apf::MeshIterator* ent_it = mesh->begin(0);
  int inode;
  while ((ent = mesh->iterate(ent_it)))
  {
    int start_dof, end_dof_plus_one;
    inode = getMdsIndex(mesh, ent);
    m3dc1_ent_getlocaldofid (&vertex_type, &inode, &fieldOrdering, &start_dof, &end_dof_plus_one);
    int startIdx = start_dof;
    assert(startIdx<num_dof);

    apf::Adjacent elements;
    getBridgeAdjacent(mesh, ent, brgType, 0, elements);
    int numAdj=0;
    for (int i=0; i<elements.getSize(); ++i)
    {
      if (!mesh->isGhost(elements[i]))
        ++numAdj;
    }

    startIdx /=bs; 
    for (int i=0; i<numBlockNode; ++i)
    {
      nnz.at(startIdx+i)=(1+numAdj)*numBlockNode;
    }
  }
  mesh->end(ent_it);

  if (bs==1) 
    MatSeqAIJSetPreallocation(*A, 0, &nnz[0]);
  else  
    MatSeqBAIJSetPreallocation(*A, bs, 0, &nnz[0]);
} 

int m3dc1_matrix::setupParaMat()
{
  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0], num_own_dof;
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;
  PetscInt mat_dim = num_own_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(MPI_COMM_WORLD, A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(*A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);

  ierr = MatSetType(*A, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
}

int m3dc1_matrix::setupSeqMat()
{
  int num_ent=m3dc1_mesh::instance()->num_local_ent[0];

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[fieldOrdering];
  int num_dof = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();

  int dofPerEnt=0;
  if (num_ent) dofPerEnt = num_dof/num_ent;

  PetscInt mat_dim = num_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF, A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(*A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
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
  ierr = MatView(*A, lab); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&lab); CHKERRQ(ierr);
}

int m3dc1_matrix::printInfo()
{
  MatInfo info;
  MatGetInfo(*A, MAT_LOCAL,&info);
  std::cout<<"Matrix "<<id<<" info "<<std::endl;
  std::cout<<"\t nz_allocated,nz_used,nz_unneeded "<<info.nz_allocated<<" "<<info.nz_used<<" "<<info.nz_unneeded<<std::endl;
  std::cout<<"\t memory mallocs "<<info.memory<<" "<<info.mallocs<<std::endl; 
  PetscInt nstash, reallocs, bnstash, breallocs;
  MatStashGetInfo(*A,&nstash,&reallocs,&bnstash,&breallocs);
  std::cout<<"\t nstash, reallocs, bnstash, breallocs "<<nstash<<" "<<reallocs<<" "<<bnstash<<" "<<breallocs<<std::endl;
}
// ***********************************
// 		MATRIX_MULTIPLY
// ***********************************
matrix_mult::matrix_mult(int i, int s, FieldID field)
: m3dc1_matrix(i,s,field), localMat(1)
{ 
  initialize();
}

int matrix_mult::initialize()
{
  // initialize matrix
  setupMat();
  preAllocate();
  int ierr = MatSetUp (*A); // "MatSetUp" sets up internal matrix data structure for the later use
  //disable error when preallocate not enough
  ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);
  ierr = MatSetOption(*A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE); CHKERRQ(ierr);
}

int matrix_mult::assemble()
{
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); 
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  mat_status = M3DC1_FIXED;
}

int matrix_mult::multiply(FieldID in_field, FieldID out_field)
{
  if (!localMat)
  {
    Vec b, c;
    copyField2PetscVec(in_field, b, get_scalar_type());
    int ierr = VecDuplicate(b, &c); CHKERRQ(ierr);
    MatMult(*A, b, c);
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
    MatGetBlockSize(*A, &bs);
    PetscScalar * array[2];
    m3dc1_field_getdataptr(&in_field, (double**)array);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double * array_org = (double*)array[0];
      array[0] = new PetscScalar[num_dof];
      for (int i=0; i<num_dof; ++i)
        array[0][i]=array_org[i];
    }
#endif
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[0],&b); CHKERRQ(ierr);
    m3dc1_field_getdataptr(&out_field, (double**)array+1);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double * array_org = (double*)array[1];
      array[1] = new PetscScalar[num_dof];
      for (int i=0; i<num_dof; ++i)
        array[1][i]=array_org[i];
    }
#endif
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[1],&c); CHKERRQ(ierr);
    ierr=VecAssemblyBegin(b);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(b);  CHKERRQ(ierr);
    ierr=VecAssemblyBegin(c);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(c);  CHKERRQ(ierr);
    MatMult(*A, b, c);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type())
    {
      double *datapt;
      m3dc1_field_getdataptr(&out_field, &datapt);
      for (int i=0; i<num_dof; ++i)
        datapt[i]=std::real(array[1][i]); 
      delete [] array[0];
      delete [] array[1];
    }
#endif
    m3dc1_field_sum(&out_field);
  }
}

// ***********************************
// 		MATRIX_SOLVE
// ***********************************

matrix_solve::matrix_solve(int i, int s, FieldID f): m3dc1_matrix(i,s,f) 
{  
  ksp = new KSP;
  kspSet=0;
  remotePidOwned=NULL;
  remoteNodeRow=NULL; // <pid, <locnode>, numAdj>
  remoteNodeRowSize=NULL;
  initialize();
}

matrix_solve::~matrix_solve()
{
  if (kspSet)
    KSPDestroy(ksp);
  delete ksp;
  MatDestroy(&remoteA);
}

int matrix_solve::initialize()
{
  // initialize matrix
  setupMat();
  preAllocate();
  if (!m3dc1_solver::instance()->assembleOption) setUpRemoteAStruct();
  int ierr = MatSetUp (*A); // "MatSetUp" sets up internal matrix data structure for the later use
  //disable error when preallocate not enough 
  ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);
  //commented per Jin's request on Nov 9, 2017
  //ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(*A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE); CHKERRQ(ierr); 
  CHKERRQ(ierr);
}

int matrix_solve::setupMat()
{
  setupParaMat();
}

int matrix_solve::preAllocate ()
{
  preAllocateParaMat();
}

int matrix_solve::reset_values() 
{ 
  int ierr = MatZeroEntries(*A); 
  //MatZeroEntries(remoteA); 
    delete remotePidOwned;
    delete remoteNodeRow;
    delete remoteNodeRowSize;
    ierr =MatDestroy(&remoteA);
    if (!m3dc1_solver::instance()->assembleOption) setUpRemoteAStruct();

  mat_status = M3DC1_NOT_FIXED; // allow matrix value modification
  //start second solve
  if(kspSet==1) 
  {
    kspSet=2;
    // Set operators, keeping the identical preconditioner matrix for
    // all linear solves.  This approach is often effective when the
    // linear systems do not change very much between successive steps.
    ierr= KSPSetReusePreconditioner(*ksp,PETSC_TRUE); CHKERRQ(ierr);
  }
  
  if (!PCU_Comm_Self())
    std::cout<<"[M3DC1 ERROR] "<<__func__<<": mat_status=M3DC1_NOT_FIXED "<<mat_status<<" kspSet="<<kspSet<<"\n";
#ifdef DEBUG_
  PetscInt rstart, rend, r_rstart, r_rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  //MatGetSize(*A, &n, NULL); -- this returns global matrix size
  MatGetOwnershipRange(*A, &rstart, &rend);
  MatGetOwnershipRange(remoteA, &r_rstart, &r_rend);

  for (PetscInt row=rstart; row<rend; ++row)
  { 
    MatGetRow(*A, row, &ncols, &cols, &vals);
    for (int i=0; i<ncols; ++i)
      assert(m3dc1_double_isequal(vals[i],0.0));
    MatRestoreRow(*A, row, &ncols, &cols, &vals); // prevent memory leak
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

  for (int i=0; i<rbsize*bs; ++i)
  {
    for (int j=0; j<cbsize*bs; ++j)
    {
      if (scalar_type==M3DC1_REAL) petscValues.at(i*cbsize*bs+j)=values[i*cbsize*bs+j];
      else
      {
#ifdef PETSC_USE_COMPLEX
        petscValues.at(i*cbsize*bs+j)=complex<double>(values[2*i*cbsize*bs+2*j], values[2*i*cbsize*bs+2*j+1]);
#else
        if (!PCU_Comm_Self())
        std::cout<<"[M3DC1 ERROR] "<<__func__<<": PETSc is configured with --with-scalar-type=real\n";
        abort();
#endif
      }
    }
  }
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, &petscValues[0], ADD_VALUES);
#else
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, (PetscScalar*)values, ADD_VALUES);
#endif
}

int matrix_solve::assemble()
{
  PetscErrorCode ierr;
  if (!m3dc1_solver::instance()->assembleOption)
  {
    ierr = MatAssemblyBegin(remoteA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(remoteA, MAT_FINAL_ASSEMBLY);
    //pass remoteA to ownnering process
    int brgType = mesh->getDimension();

    int dofPerVar = 6;
    char field_name[256];
    int num_values, value_type, total_num_dof, vertex_type=0;
    m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
    dofPerVar=total_num_dof/num_values;
 
    int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];
    PetscInt firstRow, lastRowPlusOne;
    ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);

    std::map<int, std::vector<int> >* idxSendBuff = new  std::map<int, std::vector<int> >;
    std::map<int, std::vector<int> >* idxRecvBuff = new  std::map<int, std::vector<int> >;

    std::map<int, std::vector<PetscScalar> >* valuesSendBuff = new std::map<int, std::vector<PetscScalar> >;
    std::map<int, std::vector<PetscScalar> >* valuesRecvBuff = new std::map<int, std::vector<PetscScalar> >;

    int blockMatSize = total_num_dof*total_num_dof, idxOffset, valueOffset;
    int numAdj, local_id, offset, startColumn;
    int start_global_dof_id, end_global_dof_id_plus_one;
    for (std::map<int, std::map<int, int> >::iterator it=remoteNodeRow->begin(); it!=remoteNodeRow->end(); ++it)
    {
      (*idxSendBuff)[it->first].resize(it->second.size()+(*remoteNodeRowSize)[it->first]);
      (*valuesSendBuff)[it->first].resize((*remoteNodeRowSize)[it->first]*blockMatSize);
      idxOffset=0;
      valueOffset=0;
      for (std::map<int, int>::iterator it2 =it->second.begin(); it2!=it->second.end(); ++it2)
      {
        (*idxSendBuff)[it->first].at(idxOffset++)=it2->second;
        apf::MeshEntity* ent = getMdsEntity(mesh, 0, it2->first);

        std::vector<apf::MeshEntity*> vecAdj;
        apf::Adjacent elements;
        getBridgeAdjacent(mesh, ent, brgType, 0, elements);
        for (int i=0; i<elements.getSize(); ++i)
        {
          if (!mesh->isGhost(elements[i]))
            vecAdj.push_back(elements[i]);
        }
        vecAdj.push_back(ent);
        numAdj = vecAdj.size();
        assert(numAdj==it2->second);
        std::vector<int> localNodeId(numAdj);
        std::vector<int> columns(total_num_dof*numAdj);
        for (int i=0; i<numAdj; ++i)
        {
          local_id = get_ent_localid(mesh, vecAdj.at(i));
          localNodeId.at(i)=local_id;
          
          m3dc1_ent_getglobaldofid (&vertex_type, &local_id, &fieldOrdering, &start_global_dof_id, 
               &end_global_dof_id_plus_one);
          (*idxSendBuff)[it->first].at(idxOffset++)=start_global_dof_id;
        }
        offset=0;
        for (int i=0; i<numAdj; ++i)
        {
          startColumn = localNodeId.at(i)*total_num_dof;
          for (int j=0; j<total_num_dof; ++j)
            columns.at(offset++)=startColumn+j;
        }
        ierr = MatGetValues(remoteA, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), 
               total_num_dof*numAdj, &columns[0], &(*valuesSendBuff)[it->first].at(valueOffset));
        valueOffset+=it2->second*blockMatSize;
      }
      assert(idxOffset==(*idxSendBuff)[it->first].size());
      assert(valueOffset==(*valuesSendBuff)[it->first].size());
    }
    // ierr = MatDestroy(&remoteA); // seol: shall destroy in destructor

    //send and receive message size
    int sendTag=2020;
    MPI_Request my_request[256];
    MPI_Status my_status[256];
    int requestOffset=0;
    std::map<int, std::pair<int, int> > msgSendSize;
    std::map<int, std::pair<int, int> > msgRecvSize;
    for (std::map<int, int>::iterator it=remoteNodeRowSize->begin(); it!=remoteNodeRowSize->end(); ++it)
    {
      int destPid=it->first;
      msgSendSize[destPid].first=(*idxSendBuff)[it->first].size();
      msgSendSize[destPid].second = (*valuesSendBuff)[it->first].size();
      MPI_Isend(&(msgSendSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,
                MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for (std::set<int>::iterator it=remotePidOwned->begin(); it!=remotePidOwned->end(); ++it)
    {
      int destPid=*it;
      MPI_Irecv(&(msgRecvSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,
                MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request, my_status);
    //set up receive buff
    for (std::map<int, std::pair<int, int> >::iterator it=msgRecvSize.begin(); it!= msgRecvSize.end(); ++it)
    {
      (*idxRecvBuff)[it->first].resize(it->second.first);
      (*valuesRecvBuff)[it->first].resize(it->second.second); 
    }
    msgSendSize.clear();
    msgRecvSize.clear();

    // now get data
    sendTag=9999;
    requestOffset=0;
    for (std::map<int, int>::iterator it=remoteNodeRowSize->begin(); it!=remoteNodeRowSize->end(); ++it)
    {
      int destPid=it->first;
      MPI_Isend(&((*idxSendBuff)[destPid].at(0)),(*idxSendBuff)[destPid].size(),MPI_INT,destPid,sendTag,
                MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Isend(&((*valuesSendBuff)[destPid].at(0)),sizeof(PetscScalar)*(*valuesSendBuff)[destPid].size(),
                MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for (std::set<int>::iterator it=remotePidOwned->begin(); it!=remotePidOwned->end(); ++it)
    {
      int destPid=*it;
      MPI_Irecv(&((*idxRecvBuff)[destPid].at(0)),(*idxRecvBuff)[destPid].size(),
                MPI_INT,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Irecv(&((*valuesRecvBuff)[destPid].at(0)),sizeof(PetscScalar)*(*valuesRecvBuff)[destPid].size(),
                MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request,my_status);

    for (std::map<int, std::vector<int> >::iterator it =idxSendBuff->begin(); it!=idxSendBuff->end(); ++it)
      std::vector<int>().swap(it->second);
    for (std::map<int, std::vector<PetscScalar> >::iterator it =valuesSendBuff->begin(); it!=valuesSendBuff->end(); ++it)
      std::vector<PetscScalar>().swap(it->second);

    // clean up auxiliary std container
    valuesSendBuff->clear();
    valuesSendBuff=NULL;
    idxSendBuff->clear();
    idxSendBuff = NULL;

    // now assemble the matrix
    for (std::set<int>::iterator it=remotePidOwned->begin(); it!=remotePidOwned->end(); ++it)
    {
      int destPid=*it;
      int valueOffset=0;
      int idxOffset=0;
      vector<int> & idx = (*idxRecvBuff)[destPid];
      vector<PetscScalar> & values = (*valuesRecvBuff)[destPid];
      int numValues=values.size();
      while (valueOffset<numValues)
      {
        int numAdj = idx.at(idxOffset++); 
        std::vector<int> columns(total_num_dof*numAdj);
        int offset=0;
        for (int i=0; i<numAdj; ++i, ++idxOffset)
          for (int j=0; j<total_num_dof; ++j)
            columns.at(offset++)=idx.at(idxOffset)+j;

        ierr = MatSetValues(*A, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), 
               total_num_dof*numAdj, &columns[0], &values.at(valueOffset),ADD_VALUES);

        valueOffset+=blockMatSize*numAdj;
      }
      std::vector<int>().swap((*idxRecvBuff)[destPid]);
      std::vector<PetscScalar>().swap((*valuesRecvBuff)[destPid]);
    }
    valuesRecvBuff->clear();
    valuesRecvBuff=NULL;
    idxRecvBuff->clear();
    idxRecvBuff=NULL;
  }

  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); 
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // clean up auxiliary data
  remotePidOwned->clear();
  remoteNodeRow->clear();
  remoteNodeRowSize->clear();

  remotePidOwned=NULL;
  remoteNodeRow=NULL; // <pid, <locnode>, numAdj>
  remoteNodeRowSize=NULL;

  mat_status=M3DC1_FIXED;
}

int matrix_solve:: set_bc(int row)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  MatSetValue(*A, row, row, 1.0, ADD_VALUES);
}

int matrix_solve:: set_row(int row, int numVals, int* columns, double * vals)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  for (int i=0; i<numVals; ++i)
  {
    if (get_scalar_type() == M3DC1_REAL) set_value(row, columns[i], 1, vals[i], 0);
    else set_value(row, columns[i], 1, vals[2*i], vals[2*i+1]); 
  }
}

int matrix_solve::solve(FieldID field_id)
{
  Vec x, b;
  copyField2PetscVec(field_id, b, get_scalar_type());
  int ierr = VecDuplicate(b, &x); CHKERRQ(ierr);

  if(!kspSet) setKspType();
  if(kspSet==2) {
         ierr= KSPSetOperators(*ksp,*A,*A); CHKERRQ(ierr);
         if (!PCU_Comm_Self())
           std::cout <<"\t-- Update A, Reuse Preconditioner" << std::endl;
  }


  if (!PCU_Comm_Self())
    std::cout << "~~~~~ checking number of jocobi blocks in Petsc ~~~~~" << std::endl;
  PC pc;
  PetscErrorCode pe;
  pe = KSPGetPC(*ksp, &pc); CHKERRQ(pe);
  PetscInt numBlocks;
  PCBJacobiGetTotalBlocks(pc, &numBlocks, 0);

  if (!PCU_Comm_Self())
  {
    std::cout << PCU_Comm_Self() << " total number of blocks is " << numBlocks << std::endl;
    std::cout << "~~~~~ eof checking number of jocobi blocks in Petsc ~~~~~" << std::endl;
  }


  //KSPSetUp(*ksp);
 // KSPSetUpOnBlocks(*ksp); CHKERRQ(ierr);

  ierr = KSPSolve(*ksp, b, x); 
  CHKERRQ(ierr);
//  PetscInt its;
  ierr = KSPGetIterationNumber(*ksp, &its);
  CHKERRQ(ierr);

  if (PCU_Comm_Self() == 0)
    std::cout <<"\t-- # solver iterations " << its << std::endl;
  //iterNum = its;
 
  copyPetscVec2Field(x, field_id, get_scalar_type());

  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  mat_status = M3DC1_SOLVED;
}

int matrix_solve:: setKspType()
{
  PetscErrorCode ierr;
  ierr = KSPCreate(MPI_COMM_WORLD, ksp);
  CHKERRQ(ierr);
         // Set operators, keeping the identical preconditioner matrix for
         // all linear solves.  This approach is often effective when the
         // linear systems do not change very much between successive steps.
         //ierr= KSPSetReusePreconditioner(*ksp,PETSC_TRUE); CHKERRQ(ierr);
         //if (!PCU_Comm_Self())
         //  std::cout <<"\t-- Reuse Preconditioner" << std::endl;
  ierr = KSPSetOperators(*ksp, *A, *A /*, SAME_PRECONDITIONER DIFFERENT_NONZERO_PATTERN*/); 
  CHKERRQ(ierr);
  ierr = KSPSetTolerances(*ksp, .000001, .000000001, PETSC_DEFAULT, 1000);
  CHKERRQ(ierr);

  int num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  assert(total_num_dof/num_values==C1TRIDOFNODE*(mesh->getDimension()-1));

  // if 2D problem use superlu
  if (mesh->getDimension()==2)
  {
    ierr=KSPSetType(*ksp, KSPPREONLY); CHKERRQ(ierr);
    PC pc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU); CHKERRQ(ierr);
#ifndef PETSCMASTER
// petsc-3.8.3 and older
    ierr=PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
#else
    ierr=PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST);
#endif
    CHKERRQ(ierr);
  }

  ierr = KSPSetFromOptions(*ksp); CHKERRQ(ierr);
  kspSet=1;
}

#endif //#ifdef M3DC1_PETSC
