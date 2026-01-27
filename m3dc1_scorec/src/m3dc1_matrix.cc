/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifdef M3DC1_PETSC
#include "m3dc1_matrix.h"
#include "PCU.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh.h"
#include "m3dc1_mesh.h"
#include <assert.h>
#include <iostream>
#include <vector>

#ifdef PETSC_USE_COMPLEX
#include "petscsys.h" // for PetscComplex
#include <complex>
using std::complex;
#endif

using std::vector;

// ***********************************
// 		HELPER
// ***********************************
PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy)
{
  Vec x,r;
  PetscReal norms[3];
  int stride=12,ierr;

  PetscFunctionBeginUser;
  /* Build the solution vector */
  ierr=KSPBuildSolution(ksp, NULL, &x);
  ierr=KSPBuildResidual(ksp, NULL, NULL, &r);

  ierr=VecStrideNorm(r, 0*stride, NORM_2, &norms[0]);
  ierr=VecStrideNorm(r, 1*stride, NORM_2, &norms[1]);
  ierr=VecStrideNorm(r, 2*stride, NORM_2, &norms[2]);

  ierr=PetscPrintf(PETSC_COMM_WORLD, 
		  "solve 5 %" PetscInt_FMT "th iteration: KSP Residual norm [ %1.12e  %1.12e  %1.12e  %1.12e ]\n", 
		  n, (double)norms[0], (double)norms[1], (double)norms[2], rnorm);
  ierr=VecDestroy(&r);
  PetscFunctionReturn(0); // PETSC_SUCCESS (0) to indicate success
}

void printMemStat() {
  PetscLogDouble mem, mem_max;
  PetscMemoryGetCurrentUsage(&mem);
  PetscMemoryGetMaximumUsage(&mem_max);
  std::cout
      << "\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank "
      << PCU_Comm_Self() << " current " << mem / 1e6 << std::endl;
}

int copyField2PetscVec_5(FieldID field_id, Vec petscVec, int scalar_type)
{
  int num_own_ent= m3dc1_mesh::instance()->num_own_ent[0];
  int num_own_dof=0, vertex_type=0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;

  PetscFunctionBeginUser;

  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

        if (!PCU_Comm_Self())
          std::cout<<"[M3DC1 INFO] "<<__func__<<": MatCreateVecs\n";
/*int ierr = VecCreateMPI(MPI_COMM_WORLD, num_own_dof, PETSC_DECIDE, &petscVec); 
  int ierr = VecCreate(MPI_COMM_WORLD, &petscVec); CHKERRQ(ierr);
  ierr=VecSetBlockSize(petscVec, dofPerEnt); // the blocksize=dofPerEnt, but the stridesize=dofPerEnt/3
  ierr = VecSetSizes(petscVec, num_own_dof, PETSC_DECIDE); CHKERRQ(ierr);
 * */
  int ierr = VecSetFromOptions(petscVec);CHKERRQ(ierr);
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
  PetscFunctionReturn(PETSC_SUCCESS);
//  return 0;
}

int copyField2PetscVec(FieldID field_id, Vec petscVec, int scalar_type) {
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0];
  int num_own_dof = 0, vertex_type = 0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt = 0;
  if (num_own_ent)
    dofPerEnt = num_own_dof / num_own_ent;

  int ierr;/*= VecCreateMPI(MPI_COMM_WORLD, num_own_dof, PETSC_DECIDE,
   * &petscVec);
  int ierr = VecCreate(MPI_COMM_WORLD, &petscVec);
  CHKERRQ(ierr);
  ierr = VecSetBlockSize(
      petscVec,
      dofPerEnt); // the blocksize=dofPerEnt, but the stridesize=dofPerEnt/3
  ierr = VecSetSizes(petscVec, num_own_dof, PETSC_DECIDE);
  CHKERRQ(ierr); */
  ierr = VecSetFromOptions(petscVec);
  CHKERRQ(ierr);
  VecAssemblyBegin(petscVec);

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];

  double dof_data[FIXSIZEBUFF];
  assert(sizeof(dof_data) >= dofPerEnt * 2 * sizeof(double));
  int nodeCounter = 0;

  apf::Mesh2 *mesh = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity *ent;
  int inode;
  apf::MeshIterator *ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it))) {
    inode = getMdsIndex(mesh, ent);
    if (!is_ent_original(mesh, ent))
      continue;
    nodeCounter += 1;
    int num_dof;
    m3dc1_ent_getdofdata(&vertex_type, &inode, &field_id, &num_dof, dof_data);
    assert(num_dof * (1 + scalar_type) <= sizeof(dof_data) / sizeof(double));
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid(&vertex_type, &inode, &field_id,
                             &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = 0;
    for (int i = 0; i < dofPerEnt; ++i) {
      PetscScalar value;
      if (scalar_type == M3DC1_REAL)
        value = dof_data[startIdx++];
      else {
#ifdef PETSC_USE_COMPLEX
        value =
            complex<double>(dof_data[startIdx * 2], dof_data[startIdx * 2 + 1]);
#else
        if (!PCU_Comm_Self())
          std::cout
              << "[M3DC1 ERROR] " << __func__
              << ": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
        startIdx += 1;
      }
      ierr =
          VecSetValue(petscVec, start_global_dof_id + i, value, INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }
  mesh->end(ent_it);

  assert(nodeCounter == num_own_ent);
  ierr = VecAssemblyEnd(petscVec);
  CHKERRQ(ierr);
  return 0;
}

int copyPetscVec2Field(Vec &petscVec, FieldID field_id, int scalar_type) {
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0], num_own_dof = 0,
      vertex_type = 0;
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt = 0;
  if (num_own_ent)
    dofPerEnt = num_own_dof / num_own_ent;

  std::vector<PetscInt> ix(dofPerEnt);
  std::vector<PetscScalar> values(dofPerEnt);
  std::vector<double> dof_data(dofPerEnt * (1 + scalar_type));
  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];

  int ierr, inode;

  apf::Mesh2 *mesh = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity *ent;
  apf::MeshIterator *ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it))) {
    inode = getMdsIndex(mesh, ent);
    if (!is_ent_original(mesh, ent))
      continue;
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid(&vertex_type, &inode, &field_id,
                             &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;

    for (int i = 0; i < dofPerEnt; ++i)
      ix.at(i) = startIdx + i;
    ierr = VecGetValues(petscVec, dofPerEnt, &ix[0], &values[0]);
    CHKERRQ(ierr);
    startIdx = 0;
    for (int i = 0; i < dofPerEnt; ++i) {
      if (scalar_type == M3DC1_REAL) {
#ifdef PETSC_USE_COMPLEX
        dof_data.at(startIdx++) = values.at(i).real();
#else
        dof_data.at(startIdx++) = values.at(i);
#endif
      } else {
#ifdef PETSC_USE_COMPLEX
        dof_data.at(2 * startIdx) = values.at(i).real();
        dof_data.at(2 * startIdx + 1) = values.at(i).imag();
        ++startIdx;
#else
        if (!PCU_Comm_Self())
          std::cout
              << "[M3DC1 ERROR] " << __func__
              << ": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
    m3dc1_ent_setdofdata(&vertex_type, &inode, &field_id, &dofPerEnt,
                         &dof_data[0]);
  }
  mesh->end(ent_it);

  synchronize_field(
      (*m3dc1_mesh::instance()->field_container)[field_id]->get_field());
  return 0;
}

// ***********************************
// 		M3DC1_SOLVER
// ***********************************

m3dc1_solver *m3dc1_solver::_instance = NULL;
m3dc1_solver *m3dc1_solver::instance() {
  if (_instance == NULL)
    _instance = new m3dc1_solver();
  return _instance;
}

m3dc1_solver::~m3dc1_solver() {
  if (matrix_container != NULL)
    matrix_container->clear();
  matrix_container = NULL;
  delete _instance;
}

void m3dc1_solver::add_matrix(int matrix_id, m3dc1_matrix *matrix) {
  assert(matrix_container->find(matrix_id) == matrix_container->end());
  matrix_container->insert(
      std::map<int, m3dc1_matrix *>::value_type(matrix_id, matrix));
}

m3dc1_matrix *m3dc1_solver::get_matrix(int matrix_id) {
  std::map<int, m3dc1_matrix *>::iterator mit =
      matrix_container->find(matrix_id);
  if (mit == matrix_container->end())
    return (m3dc1_matrix *)NULL;
  return mit->second;
}

// ***********************************
// 		M3DC1_MATRIX
// ***********************************

m3dc1_matrix::m3dc1_matrix(int i, int s, FieldID f)
    : mesh(m3dc1_mesh::instance()->mesh), id(i), scalar_type(s),
      fieldOrdering(f) {
  mat_status = M3DC1_NOT_FIXED;
  _A=NULL;
}

int m3dc1_matrix::destroy() {
  PetscErrorCode ierr = MatDestroy(&_A);
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

m3dc1_matrix::~m3dc1_matrix() {
  destroy();
}

int m3dc1_matrix::get_values(vector<PetscInt> &rows, vector<int> &n_columns,
                             vector<PetscInt> &columns,
                             vector<double> &values) {
  if (!mat_status) // matrix is not fixed
  {
    if (!PCU_Comm_Self())
      std::cout << __func__ << " failed: matrix " << id << " is not fixed\n";
    return M3DC1_FAILURE;
  }

#ifdef PETSC_USE_COMPLEX
  if (!PCU_Comm_Self())
    std::cout << "[M3DC1 ERROR] " << __func__
              << ": not supported for complex\n";
  return M3DC1_FAILURE;
#else
  PetscErrorCode ierr;
  PetscInt rstart, rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  ierr = MatGetOwnershipRange(_A, &rstart, &rend);
  CHKERRQ(ierr);
  for (PetscInt row = rstart; row < rend; ++row) {
    ierr = MatGetRow(_A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
    rows.push_back(row);
    n_columns.push_back(ncols);
    for (int i = 0; i < ncols; ++i) {
      columns.push_back(cols[i]);
      values.push_back(vals[i]);
    }
    ierr = MatRestoreRow(_A, row, &ncols, &cols, &vals);
    CHKERRQ(ierr);
  }
  assert(rows.size() == rend - rstart);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_matrix::set_value(
    int row, int col, int operation, double real_val,
    double imag_val) // insertion/addition with global numbering
{
  if (mat_status) // matrix is fixed
  {
    if (!PCU_Comm_Self())
      std::cout << __func__ << " failed: matrix " << id << " is fixed\n";
    return M3DC1_FAILURE;
  }

  PetscErrorCode ierr;

  if (scalar_type == M3DC1_REAL) // real
  {
    if (operation)
      ierr = MatSetValue(_A, row, col, real_val, ADD_VALUES);
    else
      ierr = MatSetValue(_A, row, col, real_val, INSERT_VALUES);
  } else // complex
  {
#ifdef PETSC_USE_COMPLEX
    PetscScalar value = complex<double>(real_val, imag_val);
    if (operation)
      ierr = MatSetValue(_A, row, col, value, ADD_VALUES);
    else
      ierr = MatSetValue(_A, row, col, value, INSERT_VALUES);
#else
    if (!PCU_Comm_Self())
      std::cout
          << "[M3DC1 ERROR] " << __func__
          << ": PETSc is not configured with --with-scalar-type=complex\n";
    abort();
#endif
  }
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::add_values(int rsize, PetscInt *rows, int csize,
                             PetscInt *columns, double *values) {
  if (mat_status) // matrix is fixed
  {
    if (!PCU_Comm_Self())
      std::cout << __func__ << " failed: matrix " << id << " is fixed\n";
    return M3DC1_FAILURE;
  }

  PetscErrorCode ierr;
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  vector<PetscScalar> petscValues(rsize * csize);
  for (int i = 0; i < rsize; ++i) {
    for (int j = 0; j < csize; ++j) {
      if (scalar_type == M3DC1_REAL)
        petscValues.at(i * csize + j) = values[i * csize + j];
      else {
#ifdef PETSC_USE_COMPLEX
        petscValues.at(i * csize + j) = complex<double>(
            values[2 * i * csize + 2 * j], values[2 * i * csize + 2 * j + 1]);
#else
        if (!PCU_Comm_Self())
          std::cout
              << "[M3DC1 ERROR] " << __func__
              << ": PETSc is not configured with --with-scalar-type=complex\n";
        abort();
#endif
      }
    }
  }

  ierr = MatSetValues(_A, rsize, rows, csize, columns, &petscValues[0],
                      ADD_VALUES);
#else
  ierr = MatSetValues(_A, rsize, rows, csize, columns, (PetscScalar *)values,
                      ADD_VALUES);
#endif
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int matrix_mult::setupMat() {
  if (localMat)
    setupSeqMat();
  else
    setupParaMat();
  return M3DC1_SUCCESS;
}

int matrix_mult::preAllocate() {
  if (localMat)
    preAllocateSeqMat();
  else
    preAllocateParaMat();
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::preAllocateParaMat() {
  int bs = 1;
  MatType type;
  MatGetType(_A, &type);

  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0], num_own_dof = 0,
      vertex_type = 0;
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  int dofPerEnt = 0;
  if (num_own_ent)
    dofPerEnt = num_own_dof / num_own_ent;

  if (strcmp(type, MATSEQAIJ) == 0 || strcmp(type, MATMPIAIJ) == 0)
    bs = 1;
  else
    bs = dofPerEnt;
  int numBlocks = num_own_dof / bs;
  int numBlockNode = dofPerEnt / bs;
  // cj    if (PCU_Comm_Self()==1) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // bs="<<bs<<" num_own_dof="<<num_own_dof<<" numBlocks="<<numBlocks<<"
  // dofPerEnt="<<dofPerEnt<<" numBlockNode="<<numBlockNode<<"\n";

  std::vector<PetscInt> dnnz(numBlocks), onnz(numBlocks);
  int startDof, endDofPlusOne;
  m3dc1_field_getowndofid(&fieldOrdering, &startDof, &endDofPlusOne);

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0], inode;
  // cj    if (PCU_Comm_Self()==1) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // startDof="<<startDof<<" endDofPlusOne="<<endDofPlusOne<<"
  // num_vtx="<<num_vtx<<"\n";

  int nnzStash = 0;
  int brgType = mesh->getDimension();

  apf::MeshEntity *ent;
  apf::MeshIterator *ent_it = mesh->begin(0);
  while ((ent = mesh->iterate(ent_it)))
  // for (int inode=0; inode<num_vtx; ++inode)
  {
    inode = getMdsIndex(mesh,
                        ent); //  ent = getMdsEntity(mesh, vertex_type, inode);
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid(&vertex_type, &inode, &fieldOrdering,
                             &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    if (start_global_dof_id < startDof ||
        start_global_dof_id >= endDofPlusOne) {
      apf::Adjacent elements;
      getBridgeAdjacent(mesh, ent, brgType, 0, elements);
      int num_elem = 0;
      for (int i = 0; i < elements.getSize(); ++i) {
        if (!mesh->isGhost(elements[i]))
          ++num_elem;
      }

      nnzStash += dofPerEnt * dofPerEnt * (num_elem + 1);
      continue;
    }
    startIdx -= startDof;
    startIdx /= bs;

    int adjNodeOwned, adjNodeGlb;
    mesh->getIntTag(ent, m3dc1_mesh::instance()->num_global_adj_node_tag,
                    &adjNodeGlb);
    mesh->getIntTag(ent, m3dc1_mesh::instance()->num_own_adj_node_tag,
                    &adjNodeOwned);
    assert(adjNodeGlb >= adjNodeOwned);

    for (int i = 0; i < numBlockNode; ++i) {
      dnnz.at(startIdx + i) = (1 + adjNodeOwned) * numBlockNode;
      onnz.at(startIdx + i) = (adjNodeGlb - adjNodeOwned) * numBlockNode;
    }
  }
  mesh->end(ent_it);

  if (bs == 1)
    MatMPIAIJSetPreallocation(_A, 0, &dnnz[0], 0, &onnz[0]);
  else
    MatMPIBAIJSetPreallocation(_A, bs, 0, &dnnz[0], 0, &onnz[0]);
  return M3DC1_SUCCESS;
}

int matrix_solve::setUpRemoteAStruct() {
  assert(remotePidOwned == NULL && remoteNodeRow == NULL &&
         remoteNodeRowSize == NULL);

  remotePidOwned = new std::set<int>;
  remoteNodeRow = new std::map<int, std::map<int, int>>;
  remoteNodeRowSize = new std::map<int, int>;

  int dofPerVar = 6, vertex_type = 0;
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type,
                      &total_num_dof);
  dofPerVar = total_num_dof / num_values;

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];

  std::vector<PetscInt> nnz_remote(num_values * num_vtx);
  int brgType = mesh->getDimension();

  apf::MeshEntity *ent;
  apf::MeshIterator *ent_it = mesh->begin(0);
  int inode;
  while ((ent = mesh->iterate(ent_it))) {
    inode = getMdsIndex(mesh, ent);
    int owner = get_ent_ownpartid(mesh, ent);
    if (owner != PCU_Comm_Self()) {
      apf::Adjacent elements;
      getBridgeAdjacent(mesh, ent, brgType, 0, elements);
      int num_elem = 0;
      for (int i = 0; i < elements.getSize(); ++i) {
        if (!mesh->isGhost(elements[i]))
          ++num_elem;
      }

      (*remoteNodeRow)[owner][inode] = num_elem + 1;
      (*remoteNodeRowSize)[owner] += num_elem + 1;
      for (int i = 0; i < num_values; ++i)
        nnz_remote[inode * num_values + i] = (num_elem + 1) * num_values;
    } else {
      apf::Copies remotes;
      mesh->getRemotes(ent, remotes);
      APF_ITERATE(apf::Copies, remotes, it)
      remotePidOwned->insert(it->first);
    }
  }
  mesh->end(ent_it);

  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF, &remoteA);
  CHKERRQ(ierr);
  ierr = MatSetType(remoteA, MATSEQBAIJ);
  CHKERRQ(ierr);
  ierr = MatSetBlockSize(remoteA, dofPerVar);
  CHKERRQ(ierr);
  ierr = MatSetSizes(remoteA, total_num_dof * num_vtx, total_num_dof * num_vtx,
                     PETSC_DECIDE, PETSC_DECIDE);
  CHKERRQ(ierr);
  MatSeqBAIJSetPreallocation(remoteA, dofPerVar, 0, &nnz_remote[0]);
  ierr = MatSetUp(remoteA);
  CHKERRQ(ierr);
  // cj  if (!PCU_Comm_Self()) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // MatCreate remoteA bs="<<dofPerVar<<" total_num_dof="<<total_num_dof<<"
  // num_vtx="<<num_vtx<<" mat_dim="<<total_num_dof*num_vtx<<"
  // num_values="<<num_values<<" total_num_dof="<<total_num_dof<<"\n";
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::preAllocateSeqMat() {
  int bs = 1, vertex_type = 0;
  MatType type;
  MatGetType(_A, &type);

  int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];
  m3dc1_field *mf = (*(m3dc1_mesh::instance()->field_container))[fieldOrdering];
  int num_dof = (m3dc1_mesh::instance()->num_local_ent[0]) *
                mf->get_num_value() * mf->get_dof_per_value();

  int dofPerEnt = 0;
  if (num_vtx)
    dofPerEnt = num_dof / num_vtx;

  if (strcmp(type, MATSEQAIJ) == 0 || strcmp(type, MATMPIAIJ) == 0)
    bs = 1;
  else
    bs = dofPerEnt;
  int numBlocks = num_dof / bs;
  int numBlockNode = dofPerEnt / bs;
  // cj//    if (!PCU_Comm_Self()) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // bs="<<bs<<" num_own_dof="<<num_dof<<" numBlocks="<<numBlocks<<"
  // dofPerEnt="<<dofPerEnt<<" numBlockNode="<<numBlockNode<<"\n";
  std::vector<PetscInt> nnz(numBlocks);
  int brgType = 2;
  if (mesh->getDimension() == 3)
    brgType = 3;

  apf::MeshEntity *ent;
  apf::MeshIterator *ent_it = mesh->begin(0);
  int inode;
  while ((ent = mesh->iterate(ent_it))) {
    int start_dof, end_dof_plus_one;
    inode = getMdsIndex(mesh, ent);
    m3dc1_ent_getlocaldofid(&vertex_type, &inode, &fieldOrdering, &start_dof,
                            &end_dof_plus_one);
    int startIdx = start_dof;
    assert(startIdx < num_dof);

    apf::Adjacent elements;
    getBridgeAdjacent(mesh, ent, brgType, 0, elements);
    int numAdj = 0;
    for (int i = 0; i < elements.getSize(); ++i) {
      if (!mesh->isGhost(elements[i]))
        ++numAdj;
    }

    startIdx /= bs;
    for (int i = 0; i < numBlockNode; ++i) {
      nnz.at(startIdx + i) = (1 + numAdj) * numBlockNode;
    }
  }
  mesh->end(ent_it);

  if (bs == 1)
    MatSeqAIJSetPreallocation(_A, 0, &nnz[0]);
  else
    MatSeqBAIJSetPreallocation(_A, bs, 0, &nnz[0]);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::setupParaMat() {
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0], num_own_dof;
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  int dofPerEnt = 0;
  if (num_own_ent)
    dofPerEnt = num_own_dof / num_own_ent;
  PetscInt mat_dim = num_own_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(MPI_COMM_WORLD, &_A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(_A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = MatSetBlockSize(_A, dofPerEnt);

  ierr = MatSetType(_A, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr= MatSetOptionsPrefix(_A,"mhard_");
  ierr = MatSetFromOptions(_A);
  CHKERRQ(ierr);
  // cj  if (!PCU_Comm_Self()) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // MatCreate A num_own_dof="<<num_own_dof<<" num_own_ent="<<num_own_ent<<"
  // dofPerEnt="<<dofPerEnt<<" mat_dim="<<mat_dim<<"\n";
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::setupSeqMat() {
  int num_ent = m3dc1_mesh::instance()->num_local_ent[0];

  m3dc1_field *mf = (*(m3dc1_mesh::instance()->field_container))[fieldOrdering];
  int num_dof = (m3dc1_mesh::instance()->num_local_ent[0]) *
                mf->get_num_value() * mf->get_dof_per_value();

  int dofPerEnt = 0;
  if (num_ent)
    dofPerEnt = num_dof / num_ent;

  PetscInt mat_dim = num_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF, &_A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(_A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = MatSetBlockSize(_A, dofPerEnt);
  ierr = MatSetFromOptions(_A);
  CHKERRQ(ierr);
  // cj  if (!PCU_Comm_Self()) std::cout<<"[M3DC1 INFO] "<<__func__<<":
  // MatCreate A num_local_dof="<<num_dof<<" num_local_ent="<<num_ent<<"
  // dofPerEnt="<<dofPerEnt<<" mat_dim="<<mat_dim<<"\n";
  return M3DC1_SUCCESS;
}

int m3dc1_matrix::write(const char *file_name) {
  PetscErrorCode ierr;
  PetscViewer lab;
  if (get_type() == 0) {
    char name_buff[256];
    sprintf(name_buff, "%s-%d.m", file_name, PCU_Comm_Self());
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name_buff, &lab);
    CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, file_name, &lab);
    CHKERRQ(ierr);
  }
  ierr = PetscViewerPushFormat(lab, PETSC_VIEWER_ASCII_MATLAB);
  CHKERRQ(ierr);
  ierr = MatView(_A, lab);
  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&lab);
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

void m3dc1_matrix::printInfo() {
  MatInfo info;
  MatGetInfo(_A, MAT_LOCAL, &info);
  std::cout << "Matrix " << id << " info " << std::endl;
  std::cout << "\t nz_allocated,nz_used,nz_unneeded " << info.nz_allocated
            << " " << info.nz_used << " " << info.nz_unneeded << std::endl;
  std::cout << "\t memory mallocs " << info.memory << " " << info.mallocs
            << std::endl;
  PetscInt nstash, reallocs, bnstash, breallocs;
  MatStashGetInfo(_A, &nstash, &reallocs, &bnstash, &breallocs);
  std::cout << "\t nstash, reallocs, bnstash, breallocs " << nstash << " "
            << reallocs << " " << bnstash << " " << breallocs << std::endl;
}
// ***********************************
// 		MATRIX_MULTIPLY
// ***********************************
matrix_mult::matrix_mult(int i, int s, FieldID field)
    : m3dc1_matrix(i, s, field), localMat(1) {
  mymatrix_id=i;
  initialize();
}

int matrix_mult::initialize() {
  // initialize matrix
  setupMat();
  preAllocate();
  int ierr = MatSetUp(_A); // "MatSetUp" sets up internal matrix data structure
                           // for the later use
  // disable error when preallocate not enough
  ierr = MatSetOption(_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRQ(ierr);
  ierr = MatSetOption(_A, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE);
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int matrix_mult::assemble() {
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  mat_status = M3DC1_FIXED;
  return M3DC1_SUCCESS;
}

int matrix_mult::multiply(FieldID in_field, FieldID out_field) {
  if (!localMat) {
    Vec b, c;
    int ierr;
    ierr=MatCreateVecs(_A, &c, &b);
    copyField2PetscVec(in_field, b, get_scalar_type());
    //int ierr;= VecDuplicate(b, &c);
    //CHKERRQ(ierr);
    MatMult(_A, b, c);
    copyPetscVec2Field(c, out_field, get_scalar_type());
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&c);
    CHKERRQ(ierr);
    return 0;
  } else {
    Vec b, c;
    m3dc1_field *mf = (*(m3dc1_mesh::instance()->field_container))[in_field];
    int num_dof = (m3dc1_mesh::instance()->num_local_ent[0]) *
                  mf->get_num_value() * mf->get_dof_per_value();

#ifdef DEBUG
    m3dc1_field *mf2 = (*(m3dc1_mesh::instance()->field_container))[out_field];
    int num_dof2 = (m3dc1_mesh::instance()->num_local_ent[0]) *
                   mf->get_num_value() * mf->get_dof_per_value();
    assert(num_dof == num_dof2);
#endif
    PetscInt bs;
    int ierr;
    MatGetBlockSize(_A, &bs);
    PetscScalar *array[2];
    m3dc1_field_getdataptr(&in_field, (double **)array);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type()) {
      double *array_org = (double *)array[0];
      array[0] = new PetscScalar[num_dof];
      for (int i = 0; i < num_dof; ++i)
        array[0][i] = array_org[i];
    }
#endif
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, bs, num_dof,
                                 (PetscScalar *)array[0], &b);
    CHKERRQ(ierr);
    m3dc1_field_getdataptr(&out_field, (double **)array + 1);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type()) {
      double *array_org = (double *)array[1];
      array[1] = new PetscScalar[num_dof];
      for (int i = 0; i < num_dof; ++i)
        array[1][i] = array_org[i];
    }
#endif
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, bs, num_dof,
                                 (PetscScalar *)array[1], &c);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(c);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(c);
    CHKERRQ(ierr);
    MatMult(_A, b, c);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&c);
    CHKERRQ(ierr);
#ifdef PETSC_USE_COMPLEX
    if (!get_scalar_type()) {
      double *datapt;
      m3dc1_field_getdataptr(&out_field, &datapt);
      for (int i = 0; i < num_dof; ++i)
        datapt[i] = std::real(array[1][i]);
      delete[] array[0];
      delete[] array[1];
    }
#endif
    m3dc1_field_sum(&out_field);
  }
  return M3DC1_SUCCESS;
}

// ***********************************
// 		MATRIX_SOLVE
// ***********************************

matrix_solve::matrix_solve(int i, int s, FieldID f) : m3dc1_matrix(i, s, f) {
  _ksp=NULL;
  _BgmgSet = 0;
  _kspSet = 0;
  remotePidOwned = NULL;
  remoteNodeRow = NULL; // <pid, <locnode>, numAdj>
  remoteNodeRowSize = NULL;
  mymatrix_id=i;
  initialize();
}

matrix_solve::~matrix_solve() {
  if (_BgmgSet) {
    int nlevels=mg_nlevels-1;
    for (int level = 0; level < nlevels; level++) {
      MatDestroy(&(mg_interp_mat[level]));
      // KSPDestroy(&(mg_level_ksp[level]));
      // PCDestroy(&(mg_level_pc[level]));
    }
    delete[] mg_interp_mat;
    delete[] mg_level_ksp;
    delete[] mg_level_pc;
    _BgmgSet = 0;
  }

  if (_kspSet) {
    KSPDestroy(&_ksp);
    _ksp=NULL;
  }

  MatDestroy(&remoteA);
}

int matrix_solve::initialize() {
  // initialize matrix
  setupMat();
  preAllocate();
  if (!m3dc1_solver::instance()->assembleOption)
    setUpRemoteAStruct();
  int ierr = MatSetUp(_A); // "MatSetUp" sets up internal matrix data structure
                           // for the later use
  // disable error when preallocate not enough
  ierr = MatSetOption(_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRQ(ierr);
  // commented per Jin's request on Nov 9, 2017
  // ierr = MatSetOption(_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
  // CHKERRQ(ierr);
  ierr = MatSetOption(_A, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE);
  CHKERRQ(ierr);
  CHKERRQ(ierr);
  return M3DC1_SUCCESS;
}

int matrix_solve::setupMat() {
  setupParaMat();
  return M3DC1_SUCCESS;
}

int matrix_solve::preAllocate() {
  preAllocateParaMat();
  return M3DC1_SUCCESS;
}

int matrix_solve::reset_values() {
  int ierr = MatZeroEntries(_A);
  // MatZeroEntries(remoteA);
  delete remotePidOwned;
  delete remoteNodeRow;
  delete remoteNodeRowSize;
  ierr = MatDestroy(&remoteA);
  if (!m3dc1_solver::instance()->assembleOption)
    setUpRemoteAStruct();

  mat_status = M3DC1_NOT_FIXED; // allow matrix value modification
  // start second solve
  if (_kspSet >= 1) {
    _kspSet = 2;
    // Set operators, keeping the identical preconditioner matrix for
    // all linear solves.  This approach is often effective when the
    // linear systems do not change very much between successive steps.
    ierr = KSPSetReusePreconditioner(_ksp, PETSC_TRUE);
    CHKERRQ(ierr);
  }

  if (!PCU_Comm_Self())
    std::cout << "[M3DC1 ERROR] " << __func__ << ": mat_status=M3DC1_NOT_FIXED "
              << mat_status << " kspSet=" << _kspSet << "\n";
#ifdef DEBUG_
  PetscInt rstart, rend, r_rstart, r_rend, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  // MatGetSize(_A, &n, NULL); -- this returns global matrix size
  MatGetOwnershipRange(_A, &rstart, &rend);
  MatGetOwnershipRange(remoteA, &r_rstart, &r_rend);

  for (PetscInt row = rstart; row < rend; ++row) {
    MatGetRow(_A, row, &ncols, &cols, &vals);
    for (int i = 0; i < ncols; ++i)
      assert(m3dc1_double_isequal(vals[i], 0.0));
    MatRestoreRow(_A, row, &ncols, &cols, &vals); // prevent memory leak
  }

  for (PetscInt row = r_rstart; row < r_rend; ++row) {
    MatGetRow(remoteA, row, &ncols, &cols, &vals);
    for (int i = 0; i < ncols; ++i)
      assert(m3dc1_double_isequal(vals[i], 0.0));
    MatRestoreRow(remoteA, row, &ncols, &cols, &vals); // prevent memory leak
  }
#endif
  return M3DC1_SUCCESS;
}

int matrix_solve::update_values() 
{ 
  int ierr = MatZeroEntries(_A); 
  //MatZeroEntries(remoteA); 
    delete remotePidOwned;
    delete remoteNodeRow;
    delete remoteNodeRowSize;
    ierr =MatDestroy(&remoteA);
    if (!m3dc1_solver::instance()->assembleOption) setUpRemoteAStruct();

  mat_status = M3DC1_NOT_FIXED; // allow matrix value modification
  //start second solve
  if(_kspSet>=1) 
  {
    _kspSet=3;
    // we set this manually (true) to get KSP to not setup, and need to turn on now
    ierr= KSPSetReusePreconditioner(_ksp,PETSC_FALSE); CHKERRQ(ierr);
    ierr= KSPSetOperators(_ksp,_A,_A); CHKERRQ(ierr);
         if (!PCU_Comm_Self())
           std::cout <<"\t-- Update A, Update Preconditioner" << std::endl;
  }
  
  if (!PCU_Comm_Self())
    std::cout<<"[M3DC1 INFO] "<<__func__<<": mat_status=M3DC1_NOT_FIXED "<<mat_status<<" kspSet="<<_kspSet<<"\n";
  return M3DC1_SUCCESS;
}

int matrix_solve::add_blockvalues(int rbsize, PetscInt *rows, int cbsize,
                                  PetscInt *columns, double *values) {
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  PetscInt bs;
  MatGetBlockSize(remoteA, &bs);
  vector<PetscScalar> petscValues(rbsize * cbsize * bs * bs);

  for (int i = 0; i < rbsize * bs; ++i) {
    for (int j = 0; j < cbsize * bs; ++j) {
      if (scalar_type == M3DC1_REAL)
        petscValues.at(i * cbsize * bs + j) = values[i * cbsize * bs + j];
      else {
#ifdef PETSC_USE_COMPLEX
        petscValues.at(i * cbsize * bs + j) =
            complex<double>(values[2 * i * cbsize * bs + 2 * j],
                            values[2 * i * cbsize * bs + 2 * j + 1]);
#else
        if (!PCU_Comm_Self())
          std::cout << "[M3DC1 ERROR] " << __func__
                    << ": PETSc is configured with --with-scalar-type=real\n";
        abort();
#endif
      }
    }
  }
  int ierr = MatSetValuesBlocked(remoteA, rbsize, rows, cbsize, columns,
                                 &petscValues[0], ADD_VALUES);
#else
  int ierr = MatSetValuesBlocked(remoteA, rbsize, rows, cbsize, columns,
                                 (PetscScalar *)values, ADD_VALUES);
#endif
  return M3DC1_SUCCESS;
}

int matrix_solve::assemble() {
  PetscErrorCode ierr;
  if (!m3dc1_solver::instance()->assembleOption) {
    ierr = MatAssemblyBegin(remoteA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(remoteA, MAT_FINAL_ASSEMBLY);
    // pass remoteA to ownnering process
    int brgType = mesh->getDimension();

    int dofPerVar = 6;
    char field_name[256];
    int num_values, value_type, total_num_dof, vertex_type = 0;
    m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type,
                        &total_num_dof);
    dofPerVar = total_num_dof / num_values;

    int num_vtx = m3dc1_mesh::instance()->num_local_ent[0];
    PetscInt firstRow, lastRowPlusOne;
    ierr = MatGetOwnershipRange(_A, &firstRow, &lastRowPlusOne);

    std::map<int, std::vector<int>> *idxSendBuff =
        new std::map<int, std::vector<int>>;
    std::map<int, std::vector<int>> *idxRecvBuff =
        new std::map<int, std::vector<int>>;

    std::map<int, std::vector<PetscScalar>> *valuesSendBuff =
        new std::map<int, std::vector<PetscScalar>>;
    std::map<int, std::vector<PetscScalar>> *valuesRecvBuff =
        new std::map<int, std::vector<PetscScalar>>;

    int blockMatSize = total_num_dof * total_num_dof, idxOffset, valueOffset;
    int numAdj, local_id, offset, startColumn;
    int start_global_dof_id, end_global_dof_id_plus_one;
    for (std::map<int, std::map<int, int>>::iterator it =
             remoteNodeRow->begin();
         it != remoteNodeRow->end(); ++it) {
      (*idxSendBuff)[it->first].resize(it->second.size() +
                                       (*remoteNodeRowSize)[it->first]);
      (*valuesSendBuff)[it->first].resize((*remoteNodeRowSize)[it->first] *
                                          blockMatSize);
      idxOffset = 0;
      valueOffset = 0;
      for (std::map<int, int>::iterator it2 = it->second.begin();
           it2 != it->second.end(); ++it2) {
        (*idxSendBuff)[it->first].at(idxOffset++) = it2->second;
        apf::MeshEntity *ent = getMdsEntity(mesh, 0, it2->first);

        std::vector<apf::MeshEntity *> vecAdj;
        apf::Adjacent elements;
        getBridgeAdjacent(mesh, ent, brgType, 0, elements);
        for (int i = 0; i < elements.getSize(); ++i) {
          if (!mesh->isGhost(elements[i]))
            vecAdj.push_back(elements[i]);
        }
        vecAdj.push_back(ent);
        numAdj = vecAdj.size();
        assert(numAdj == it2->second);
        std::vector<int> localNodeId(numAdj);
        std::vector<PetscInt> columns(total_num_dof * numAdj);
        for (int i = 0; i < numAdj; ++i) {
          local_id = get_ent_localid(mesh, vecAdj.at(i));
          localNodeId.at(i) = local_id;

          m3dc1_ent_getglobaldofid(&vertex_type, &local_id, &fieldOrdering,
                                   &start_global_dof_id,
                                   &end_global_dof_id_plus_one);
          (*idxSendBuff)[it->first].at(idxOffset++) = start_global_dof_id;
        }
        offset = 0;
        for (int i = 0; i < numAdj; ++i) {
          startColumn = localNodeId.at(i) * total_num_dof;
          for (int j = 0; j < total_num_dof; ++j)
            columns.at(offset++) = startColumn + j;
        }
        ierr = MatGetValues(remoteA, total_num_dof,
                            &columns.at(total_num_dof * (numAdj - 1)),
                            total_num_dof * numAdj, &columns[0],
                            &(*valuesSendBuff)[it->first].at(valueOffset));
        valueOffset += it2->second * blockMatSize;
      }
      assert(idxOffset == (*idxSendBuff)[it->first].size());
      assert(valueOffset == (*valuesSendBuff)[it->first].size());
    }
    // ierr = MatDestroy(&remoteA); // seol: shall destroy in destructor

    // send and receive message size
    int sendTag = 2020;
    MPI_Request my_request[256];
    MPI_Status my_status[256];
    int requestOffset = 0;
    std::map<int, std::pair<int, int>> msgSendSize;
    std::map<int, std::pair<int, int>> msgRecvSize;
    for (std::map<int, int>::iterator it = remoteNodeRowSize->begin();
         it != remoteNodeRowSize->end(); ++it) {
      int destPid = it->first;
      msgSendSize[destPid].first = (*idxSendBuff)[it->first].size();
      msgSendSize[destPid].second = (*valuesSendBuff)[it->first].size();
      MPI_Isend(&(msgSendSize[destPid]), sizeof(std::pair<int, int>), MPI_BYTE,
                destPid, sendTag, MPI_COMM_WORLD,
                &(my_request[requestOffset++]));
    }
    assert(requestOffset < 256);
    for (std::set<int>::iterator it = remotePidOwned->begin();
         it != remotePidOwned->end(); ++it) {
      int destPid = *it;
      MPI_Irecv(&(msgRecvSize[destPid]), sizeof(std::pair<int, int>), MPI_BYTE,
                destPid, sendTag, MPI_COMM_WORLD,
                &(my_request[requestOffset++]));
    }
    assert(requestOffset < 256);
    MPI_Waitall(requestOffset, my_request, my_status);
    // set up receive buff
    for (std::map<int, std::pair<int, int>>::iterator it = msgRecvSize.begin();
         it != msgRecvSize.end(); ++it) {
      (*idxRecvBuff)[it->first].resize(it->second.first);
      (*valuesRecvBuff)[it->first].resize(it->second.second);
    }
    msgSendSize.clear();
    msgRecvSize.clear();

    // now get data
    sendTag = 9999;
    requestOffset = 0;
    for (std::map<int, int>::iterator it = remoteNodeRowSize->begin();
         it != remoteNodeRowSize->end(); ++it) {
      int destPid = it->first;
      MPI_Isend(&((*idxSendBuff)[destPid].at(0)),
                (*idxSendBuff)[destPid].size(), MPI_INT, destPid, sendTag,
                MPI_COMM_WORLD, &(my_request[requestOffset++]));
      MPI_Isend(&((*valuesSendBuff)[destPid].at(0)),
                sizeof(PetscScalar) * (*valuesSendBuff)[destPid].size(),
                MPI_BYTE, destPid, sendTag, MPI_COMM_WORLD,
                &(my_request[requestOffset++]));
    }
    assert(requestOffset < 256);
    for (std::set<int>::iterator it = remotePidOwned->begin();
         it != remotePidOwned->end(); ++it) {
      int destPid = *it;
      MPI_Irecv(&((*idxRecvBuff)[destPid].at(0)),
                (*idxRecvBuff)[destPid].size(), MPI_INT, destPid, sendTag,
                MPI_COMM_WORLD, &(my_request[requestOffset++]));
      MPI_Irecv(&((*valuesRecvBuff)[destPid].at(0)),
                sizeof(PetscScalar) * (*valuesRecvBuff)[destPid].size(),
                MPI_BYTE, destPid, sendTag, MPI_COMM_WORLD,
                &(my_request[requestOffset++]));
    }
    assert(requestOffset < 256);
    MPI_Waitall(requestOffset, my_request, my_status);

    for (std::map<int, std::vector<int>>::iterator it = idxSendBuff->begin();
         it != idxSendBuff->end(); ++it)
      std::vector<int>().swap(it->second);
    for (std::map<int, std::vector<PetscScalar>>::iterator it =
             valuesSendBuff->begin();
         it != valuesSendBuff->end(); ++it)
      std::vector<PetscScalar>().swap(it->second);

    // clean up auxiliary std container
    valuesSendBuff->clear();
    valuesSendBuff = NULL;
    idxSendBuff->clear();
    idxSendBuff = NULL;

    // now assemble the matrix
    for (std::set<int>::iterator it = remotePidOwned->begin();
         it != remotePidOwned->end(); ++it) {
      int destPid = *it;
      int valueOffset = 0;
      int idxOffset = 0;
      vector<int> &idx = (*idxRecvBuff)[destPid];
      vector<PetscScalar> &values = (*valuesRecvBuff)[destPid];
      int numValues = values.size();
      while (valueOffset < numValues) {
        int numAdj = idx.at(idxOffset++);
        std::vector<PetscInt> columns(total_num_dof * numAdj);
        int offset = 0;
        for (int i = 0; i < numAdj; ++i, ++idxOffset)
          for (int j = 0; j < total_num_dof; ++j)
            columns.at(offset++) = idx.at(idxOffset) + j;

        ierr = MatSetValues(_A, total_num_dof,
                            &columns.at(total_num_dof * (numAdj - 1)),
                            total_num_dof * numAdj, &columns[0],
                            &values.at(valueOffset), ADD_VALUES);

        valueOffset += blockMatSize * numAdj;
      }
      std::vector<int>().swap((*idxRecvBuff)[destPid]);
      std::vector<PetscScalar>().swap((*valuesRecvBuff)[destPid]);
    }
    valuesRecvBuff->clear();
    valuesRecvBuff = NULL;
    idxRecvBuff->clear();
    idxRecvBuff = NULL;
  }

  ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  // clean up auxiliary data
  remotePidOwned->clear();
  remoteNodeRow->clear();
  remoteNodeRowSize->clear();

  remotePidOwned = NULL;
  remoteNodeRow = NULL; // <pid, <locnode>, numAdj>
  remoteNodeRowSize = NULL;

  mat_status = M3DC1_FIXED;
  return M3DC1_SUCCESS;
}

int matrix_solve::set_bc(int row) {
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(_A, &firstRow, &lastRowPlusOne);
  assert(row >= firstRow && row < lastRowPlusOne);
#endif
  MatSetValue(_A, row, row, 1.0, ADD_VALUES);
  return M3DC1_SUCCESS;
}

int matrix_solve::set_row(int row, int numVals, int *columns, double *vals) {
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(_A, &firstRow, &lastRowPlusOne);
  assert(row >= firstRow && row < lastRowPlusOne);
#endif
  for (int i = 0; i < numVals; ++i) {
    if (get_scalar_type() == M3DC1_REAL)
      set_value(row, columns[i], 1, vals[i], 0);
    else
      set_value(row, columns[i], 1, vals[2 * i], vals[2 * i + 1]);
  }
  return M3DC1_SUCCESS;
}

int matrix_solve::solve(FieldID field_id) {
  Vec x, b;
  int ierr;
  ierr=MatCreateVecs(_A, &x, &b);
  copyField2PetscVec_5(field_id, b, get_scalar_type());
  //int ierr;= VecDuplicate(b, &x);
  //CHKERRQ(ierr);

  if (!_kspSet)
    setKspType();
  if (_kspSet == 2) {
    ierr = KSPSetOperators(_ksp, _A, _A);
    CHKERRQ(ierr);
    if (!PCU_Comm_Self())
      std::cout << "\t-- Update A, Reuse Preconditioner" << std::endl;
  }

  // KSPSetUp(_ksp);
  // KSPSetUpOnBlocks(_ksp); CHKERRQ(ierr);

  ierr = KSPSolve(_ksp, b, x);
  CHKERRQ(ierr);
  //  PetscInt its;
  ierr = KSPGetIterationNumber(_ksp, &its);
  CHKERRQ(ierr);

  if (PCU_Comm_Self() == 0)
    std::cout << "\t-- # solver iterations " << its << std::endl;
  // iterNum = its;

  copyPetscVec2Field(x, field_id, get_scalar_type());

  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  ierr = VecDestroy(&x);
  CHKERRQ(ierr);
  mat_status = M3DC1_SOLVED;
  return M3DC1_SUCCESS;
}

// solve with non-zero initial guess
int matrix_solve::solve_with_guess(FieldID field_id, FieldID xVec_guess) {
  Vec x, b;
  copyField2PetscVec(field_id, b, get_scalar_type());
  copyField2PetscVec(xVec_guess, x, get_scalar_type());
  int ierr;
  KSPType ksptype;

  if (!_kspSet)
    setKspType();
  if (_kspSet == 2) {
    ierr = KSPSetOperators(_ksp, _A, _A);
    CHKERRQ(ierr);
    if (!PCU_Comm_Self())
      std::cout << "\t-- Update A, Reuse Preconditioner" << std::endl;
  }

  ierr = KSPGetType(_ksp, &ksptype);
  CHKERRQ(ierr);
  if (strcmp(ksptype, "preonly") == 0) {
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE);
    CHKERRQ(ierr);
    if (PCU_Comm_Self() == 0)
      std::cout << "\t Due to ksptype=\"preonly\", the initial guess is set to "
                   "be Zero."
                << std::endl;
  } else {
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
    CHKERRQ(ierr);
  }
  ierr = KSPSolve(_ksp, b, x);
  CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(_ksp, &its);
  CHKERRQ(ierr);

  if (PCU_Comm_Self() == 0)
    std::cout << "\t-- # solver iterations " << its << std::endl;

  copyPetscVec2Field(x, field_id, get_scalar_type());

  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  ierr = VecDestroy(&x);
  CHKERRQ(ierr);
  mat_status = M3DC1_SOLVED;
  return M3DC1_SUCCESS;
}

int matrix_solve::setKspType() {
  PetscErrorCode ierr;
  ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);
  CHKERRQ(ierr);

  // mgsolve is turned on only if the solve is 5 or 17
  PetscInt ss = 2, mgsolve[2];
  mgsolve[0] = -1;
  mgsolve[1] = -1;
  ierr = PetscOptionsGetIntArray(NULL, NULL, "-mgsolve", mgsolve, &ss, NULL);
  CHKERRQ(ierr);
  if (mymatrix_id == mgsolve[0] || mymatrix_id == mgsolve[1]) {
    if (!PCU_Comm_Self())
      std::cout << "[M3DC1 INFO] " << __func__ << ": matrix " << mymatrix_id
                << " is going to use BGMG preconditioner" << "\n";
    if (!_BgmgSet)
      setBgmgType();
  }

  // Set operators, keeping the identical preconditioner matrix for
  // all linear solves.  This approach is often effective when the
  // linear systems do not change very much between successive steps.
  // ierr= KSPSetReusePreconditioner(_ksp,PETSC_TRUE); CHKERRQ(ierr);
  // if (!PCU_Comm_Self())
  //  std::cout <<"\t-- Reuse Preconditioner" << std::endl;
  ierr = KSPSetOperators(
      _ksp, _A, _A /*, SAME_PRECONDITIONER DIFFERENT_NONZERO_PATTERN*/);
  CHKERRQ(ierr);
  ierr = KSPSetTolerances(_ksp, .000001, .000000001, PETSC_DEFAULT, 1000);
  CHKERRQ(ierr);

  PetscBool flg = PETSC_FALSE;
  ierr=PetscOptionsGetBool(NULL, NULL, "-ksp_monitor_5", &flg, NULL);
  if (flg && mymatrix_id==5) ierr=KSPMonitorSet(_ksp, MyKSPMonitor, NULL, 0);

  int num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type,
                      &total_num_dof);
  assert(total_num_dof / num_values ==
         C1TRIDOFNODE * (mesh->getDimension() - 1));

  // if 2D problem use superlu
  if (mesh->getDimension() == 2) {
    ierr = KSPSetType(_ksp, KSPPREONLY);
    CHKERRQ(ierr);
    PC pc;
    ierr = KSPGetPC(_ksp, &pc);
    CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU);
    CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);
    CHKERRQ(ierr);
  } else { /* conflict with mg settings, to be fixed later
           ierr= KSPSetType(_ksp,KSPFGMRES);
           PC pc;
           ierr= KSPGetPC(_ksp,&pc);
           ierr= PCSetType(pc,PCBJACOBI);

           int nplane, *blks;
           m3dc1_plane_getnum(&nplane);
           ierr=PetscMalloc1(nplane, &blks);

           int num_own_dof, global_dim, plane_dim;
           m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
           MPI_Allreduce(&num_own_dof, &global_dim, 1, MPI_INTEGER, MPI_SUM,
           MPI_COMM_WORLD ); plane_dim=global_dim/nplane;

           for (int i = 0; i < nplane; i++) blks[i] = plane_dim;
           ierr=PCBJacobiSetTotalBlocks(pc, nplane, blks);
           ierr=PetscFree(blks);
           */
    if (mymatrix_id == 5)
      ierr = KSPSetOptionsPrefix(_ksp, "hard_");
      ierr = MatViewFromOptions(_A, NULL, "-A_view");
  }

  ierr = KSPSetFromOptions(_ksp);
  CHKERRQ(ierr);
  _kspSet = 1;
  return M3DC1_SUCCESS;
}

int matrix_solve::setBgmgType() {
  // if (mesh->getDimension()!=3 || mymatrix_id!=5) return 0;
  if (mesh->getDimension() != 3)
    return 0;

  PetscInt ierr;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //          Setup Level Data
  // 0 is always the coarsest level; n-1 is the finest.  This is backward
  // compared to what some people do.
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int nplane, nlevels; //=16;
  m3dc1_plane_getnum(&nplane);
  if (nplane % 2) {
    if (!PCU_Comm_Self())
      std::cout << "[M3DC1 INFO] " << __func__ << ": odd number of planes "
                << nplane << " could not be coarsened " << "\n";
    return 0;
  }

  // set default bgmg levels to 2
  nlevels = 2;
  // or to many levels defined by the size of nplanes
  nlevels = PetscInt(log(PetscReal(nplane)) / log(2.) + 1.e-5);
  nlevels++;
  // or to levels given as the srun command line option, for example
  // "-mg_nlevels 3"
  ierr = PetscOptionsGetInt(NULL, NULL, "-mg_nlevels", &nlevels, NULL);
  if (!PCU_Comm_Self())
    std::cout << "[M3DC1 INFO] " << __func__
              << ": f requested total_mg_nlevels=" << nlevels << "\n";

  int *nplanes, *mg_nplanes; // number of planes per mg level
  // number of planes for each level
  ierr = PetscMalloc1(nlevels, &nplanes);
  // finest level
  mg_nlevels = 2;
  nplanes[nlevels - 1] = nplane;
  if (!PCU_Comm_Self())
    std::cout << "[M3DC1 INFO] " << __func__ << ": fine level " << nlevels - 1
              << " has " << nplanes[nlevels - 1] << " planes" << "\n";
  // rest of the levels
  for (int level = nlevels - 2; level >= 0; --level) {
    nplanes[level] = nplanes[level + 1] / 2;
    if ((mg_nlevels >= nlevels) || (nplanes[level] % 2)) {
      if (!PCU_Comm_Self())
        std::cout << "[M3DC1 INFO] " << __func__ << ": odd number of planes "
                  << nplanes[level] << " could not be coarsened " << "\n";
      break;
    } else {
      mg_nlevels++;
      if (!PCU_Comm_Self())
        std::cout << "[M3DC1 INFO] " << __func__ << ": level " << level
                  << " has " << nplanes[level] << " planes" << "\n";
    }
  }
  if (!PCU_Comm_Self())
    std::cout << "[M3DC1 INFO] " << __func__
              << ": f actual total_mg_nlevels=" << mg_nlevels << "\n";
  ierr = PetscMalloc1(mg_nlevels, &mg_nplanes);
  for (int level = mg_nlevels - 1; level >= 0; --level) {
    mg_nplanes[level] = nplanes[level + (nlevels - mg_nlevels)];
  }
  ierr = PetscFree(nplanes);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //          Create KSP and set multigrid options in PC
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PC pcmg;
#ifdef RICHARDSON
  PC pcksp;
  ierr= KSPGetPC(*ksp,&pcksp);
  ierr= PCSetType(pcksp,PCKSP);
  KSP ksprich;
  ierr=PCKSPGetKSP(pcksp, &ksprich);
  ierr=KSPSetType(ksprich, KSPRICHARDSON);
  ierr=KSPSetTolerances(ksprich, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
  ierr=KSPSetNormType(ksprich, KSP_NORM_NONE);
  ierr=KSPSetConvergenceTest(ksprich, KSPConvergedSkip, NULL, NULL);

  ierr=KSPGetPC(ksprich, &pcmg);
#else
    ierr=KSPGetPC(_ksp, &pcmg);
#endif
  //    ierr= KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetOptionsPrefix(_ksp, "hard_");
  ierr = KSPGetPC(_ksp, &pcmg);
  ierr = PCSetType(pcmg, PCMG);
  ierr = PCMGSetLevels(pcmg, mg_nlevels, NULL);
  ierr = PCMGSetType(pcmg, PC_MG_MULTIPLICATIVE);
  ierr = PCMGSetGalerkin(pcmg, PC_MG_GALERKIN_PMAT);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //          Create Interpolation Operators from level-1 to level
  //          mat_dim=endDofPlusOne-startDof=Iend-Istartnum_own_dof=Iendc-Istartc
  //          global_dim=mglobal=nglobal
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  PetscInt Istart, Iend;
  PetscInt Istartc, Iendc;
  PetscInt mglobal, nglobal;
  ierr = MatGetOwnershipRange(_A, &Istart, &Iend);
  ierr = MatGetOwnershipRangeColumn(_A, &Istartc, &Iendc);
  // debug ierr= MatGetSize(_A, &mglobal, &nglobal);
  // num_own_ent is the number of vertices owned by each process
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0], num_own_dof;
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  // dof pet ent of num_own_ent
  int dofPerEnt = 0; //=12,24,36
  if (num_own_ent)
    dofPerEnt = num_own_dof / num_own_ent;
  PetscInt mat_dim = num_own_dof, global_dim, plane_dim;
  int startDof, endDofPlusOne;
  m3dc1_field_getowndofid(&fieldOrdering, &startDof, &endDofPlusOne);
  // equivalent to mat_dim=endDofPlusOne-startDof=Iend-Istart (local dim)
  MPI_Allreduce(&mat_dim, &global_dim, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
  plane_dim = global_dim / nplane;

  PetscInt myrank, maxrank, npart, planeid, partitionid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &maxrank);
  npart = maxrank / nplane;
  // planeid equivalent to PetscInt(myrank/npart);
  m3dc1_plane_getid(&planeid);
  partitionid = myrank % npart;

  // reset runtime solver options on the finest level
  // this is the first place to check
  char mg_pcbj[64], mg_pcbjblocknumber[8];
  {
    int level = mg_nlevels - 1;
    sprintf(mg_pcbj, "%s%d%s", "-hard_mg_levels_", level, "_pc_bjacobi_blocks");
    sprintf(mg_pcbjblocknumber, "%d", mg_nplanes[level]);
    PetscOptionsSetValue(NULL, mg_pcbj, mg_pcbjblocknumber);
  }

  int irow, icol, icol2, irow_end, icol_end, icol2_end;

  mg_interp_mat = new Mat[mg_nlevels - 1];
  mg_level_ksp = new KSP[mg_nlevels - 1];
  mg_level_pc = new PC[mg_nlevels - 1];

  int *mg_planeid, factor;
  ierr = PetscMalloc1(mg_nlevels - 1, &mg_planeid);
  for (int level = mg_nlevels - 2; level >= 0; --level) {
    factor = pow(2, mg_nlevels - 2 - level);
    mg_planeid[level] = planeid / factor;
  }

  int *mg_offset;
  ierr = PetscMalloc1(mg_nlevels - 1, &mg_offset);
  for (int level = mg_nlevels - 2; level >= 0; --level) {
    mg_offset[level] = PetscInt((mg_planeid[level] + 1) / 2);
  }

  int *mg_num_own_ent;              /*example 5 levels: 0 1 2 3 4 */
  int *mg_start_ent, mg_start_entx; /*example 5 levels: 0 1 2 3 */
  ierr = PetscMalloc1(mg_nlevels, &mg_num_own_ent);
  ierr = PetscMalloc1(mg_nlevels - 1, &mg_start_ent);
  mg_num_own_ent[mg_nlevels - 1] = num_own_ent;
  for (int level = mg_nlevels - 2; level >= 0; --level) {
    mapping(maxrank, myrank, mg_nlevels - 1 - level, mg_num_own_ent[level + 1],
            mg_nplanes[level], planeid, npart, partitionid,
            &mg_num_own_ent[level], &mg_start_ent[level], &mg_start_entx);
  }

  for (int level = 0; level < mg_nlevels - 1; level++) {
    ierr = MatCreate(PETSC_COMM_WORLD, &mg_interp_mat[level]);
    ierr= MatSetOptionsPrefix(mg_interp_mat[level],"mhard_");

    ierr = MatSetSizes(
        mg_interp_mat[level], mg_num_own_ent[level + 1] * dofPerEnt,
        mg_num_own_ent[level] * dofPerEnt, plane_dim * mg_nplanes[level + 1],
        plane_dim * mg_nplanes[level]);

    ierr = MatSetType(mg_interp_mat[level], MATMPIAIJ);
    ierr = MatSetBlockSize(mg_interp_mat[level], dofPerEnt);
    ierr = MatSetFromOptions(mg_interp_mat[level]);
    ierr = MatSetUp(mg_interp_mat[level]);
    ierr = MatZeroEntries(mg_interp_mat[level]);
    ierr = MatGetOwnershipRangeColumn(mg_interp_mat[level], &Istartc, &Iendc);

    // hermite cubic extra term 1/8 delta (the span of elements on the coarse
    // mesh)
    PetscReal hc = M_PI / nplane / 2.;
    int iv, idof;
    for (iv = 0; iv < mg_num_own_ent[level + 1]; iv++) {
      for (idof = 0; idof < dofPerEnt; idof++) {
        irow =
            mg_start_ent[level] * dofPerEnt + iv * dofPerEnt + idof /*Istart*/;
        irow_end = mg_start_ent[level] * dofPerEnt +
                   mg_num_own_ent[level + 1] * dofPerEnt;
        icol = irow - plane_dim * mg_offset[level];
        icol_end = irow_end - plane_dim * mg_offset[level];
        icol2 = icol + plane_dim;
        icol2_end = icol_end + plane_dim;
        if ((mg_planeid[level] + 1) == mg_nplanes[level + 1])
          icol2 = icol2 % plane_dim;
        if ((mg_planeid[level] + 1) == mg_nplanes[level + 1])
          icol2_end = 1 + (icol2_end - 1) % plane_dim;
        if (!(mg_planeid[level] % 2)) {
          ierr = MatSetValue(mg_interp_mat[level], irow, icol, 1., ADD_VALUES);
        } else {
          ierr = MatSetValue(mg_interp_mat[level], irow, icol, .5, ADD_VALUES);
          ierr = MatSetValue(mg_interp_mat[level], irow, icol2, .5, ADD_VALUES);
          /*
          if(irow%36>=0  && irow%36<=5 ||
                          irow%36>=12 && irow%36<=17 ||
                          irow%36>=24 && irow%36<=29)  {
                  ierr= MatSetValue(mg_interp_mat[level],irow, 6+icol,hc,
          ADD_VALUES); ierr= MatSetValue(mg_interp_mat[level],irow, 6+icol2,-hc,
          ADD_VALUES);
          }
           */
        }
      }
    }

    ierr = MatAssemblyBegin(mg_interp_mat[level], MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(mg_interp_mat[level], MAT_FINAL_ASSEMBLY);

    //   runtime options:
    //   -A_view ascii:stdout
    //   -A_view ascii[:[filename][:[ascii_info][:append]]]
    //   -A_view ascii[:[filename][:[ascii_info_detail][:append]]]
    //   -A_view ascii[:[filename][:[ascii_matlab][:append]]]
    //   -A_view binary[:[filename][:[ascii_info][:append]]]
    //   -A_view binary[:[filename][:[ascii_info_detail][:append]]]
    //   -A_view binary[:[filename][:[ascii_matlab][:append]]]
    if (level == 0)
      ierr = MatViewFromOptions(mg_interp_mat[level], NULL, "-I0_view");
    if (level == 1)
      ierr = MatViewFromOptions(mg_interp_mat[level], NULL, "-I1_view");

    // Set Interpolation Operators

    int ilevel = level + 1;
    ierr = PCMGSetInterpolation(pcmg, ilevel, mg_interp_mat[level]);

    // Set Smoothers on each level

    ierr = PCMGGetSmoother(pcmg, level, &(mg_level_ksp[level]));
    ierr = KSPGetPC(mg_level_ksp[level], &(mg_level_pc[level]));
    ierr = KSPSetType(mg_level_ksp[level], KSPGMRES);
    ierr = KSPSetTolerances(mg_level_ksp[level], PETSC_DEFAULT, PETSC_DEFAULT,
                            PETSC_DEFAULT, 3);
    ierr = PCSetType(mg_level_pc[level], PCBJACOBI);

    int *blks;
    ierr = PetscMalloc1(mg_nplanes[level], &blks);
    for (int i = 0; i < mg_nplanes[level]; i++)
      blks[i] = plane_dim;
    ierr = PCBJacobiSetTotalBlocks(mg_level_pc[level], mg_nplanes[level], blks);
    ierr = PetscFree(blks);

    // todo
    // KSP        *subksp;  /* array of local KSP contexts on this processor */
    // PC          subpc;   /* PC context for subdomain */
    // ierr=KSPSetUp(mg_level_ksp[level]);
    // ierr=PCBJacobiGetSubKSP(mg_level_pc[level], NULL, NULL, &subksp);
    /*for (int i = 0; i < mg_nplanes[level]; i++) {
        ierr=KSPGetPC(subksp[i], &subpc);
        ierr=PCSetType(subpc, PCLU);
        ierr=PCFactorSetMatSolverType(subpc, MATSOLVERMUMPS);
        ierr=KSPSetType(subksp[i], KSPPREONLY);
    }*/
  }

  ierr = PetscFree(mg_start_ent);
  ierr = PetscFree(mg_num_own_ent);
  ierr = PetscFree(mg_offset);
  ierr = PetscFree(mg_planeid);
  ierr = PetscFree(mg_nplanes);
  _BgmgSet = 1;
  return M3DC1_SUCCESS;
}

int matrix_solve::mapping(int maxrank, int myrank, int factor, int num_own_ent,
                          int nplane, int planeid, int npartition,
                          int partitionid, int *num_own_ent_level,
                          int *start_ent, int *start_ent_level) {
  int i, *num_all_ent, *num_all_ent_level;

  num_all_ent = (int *)malloc(maxrank * sizeof(int));

  MPI_Allgather(&num_own_ent, 1, MPI_INT, num_all_ent, 1, MPI_INT,
                PETSC_COMM_WORLD);

  int iloc, imod = pow(2, factor);

  for (int i = 0; i < imod / 2; i++) {
    if ((planeid % imod) ==
        (0 +
         2 * i)) { // from the even planeid, keep the first part on this plane
      if (!(myrank % 2)) {
        iloc = myrank - i * npartition - partitionid / 2;
        *num_own_ent_level = num_all_ent[iloc] / 2;
      }
      if ((myrank % 2)) {
        iloc = (myrank - 1) - i * npartition - partitionid / 2;
        *num_own_ent_level = num_all_ent[iloc] - num_all_ent[iloc] / 2;
      }
    } else if ((planeid % imod) ==
               (1 + 2 * i)) { // from the old planeid, sent the second part on
                              // this plane
      if (!(myrank % 2)) {
        iloc = myrank - i * npartition - npartition / 2 - partitionid / 2;
        *num_own_ent_level = num_all_ent[iloc] / 2;
      }
      if ((myrank % 2)) {
        iloc = (myrank - 1) - i * npartition - npartition / 2 - partitionid / 2;
        *num_own_ent_level = num_all_ent[iloc] - num_all_ent[iloc] / 2;
      }
    }
  }
  int all_ent = 0;
  for (i = 0; i < maxrank; i++)
    all_ent = all_ent + num_all_ent[i];

  num_all_ent_level = (int *)malloc(maxrank * sizeof(int));
  MPI_Allgather(num_own_ent_level, 1, MPI_INT, num_all_ent_level, 1, MPI_INT,
                PETSC_COMM_WORLD);

  *start_ent = 0;
  for (i = 0; i < myrank; i++)
    *start_ent = *start_ent + num_all_ent[i];

  *start_ent_level = 0;
  for (i = 0; i < myrank; i++)
    *start_ent_level = *start_ent_level + num_all_ent_level[i];

  free(num_all_ent_level);
  free(num_all_ent);
  return 0;
}

#endif // #ifdef M3DC1_PETSC
