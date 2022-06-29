/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifdef M3DC1_PETSC
#include "m3dc1_matrix.h"
#include "m3dc1_scorec.h"
#include "apf.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "apfMesh.h"
#include "apfMDS.h"
#include <vector>
#include <set>
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
using std::set;

void printMemStat()
{
  PetscLogDouble mem, mem_max;
  PetscMemoryGetCurrentUsage(&mem);
  PetscMemoryGetMaximumUsage(&mem_max);
  std::cout<<"\tMemory usage (MB) reported by PetscMemoryGetCurrentUsage: Rank "<<PCU_Comm_Self()<<" current "<<mem/1e6<<std::endl;
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

// ***********************************
// 		M3DC1_MATRIX
// ***********************************

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


int matrix_solve::initialize()
{
  // initialize matrix
  setupMat();
  preAllocate();
  if(!m3dc1_solver::instance()->assembleOption) setUpRemoteAStruct();
  int ierr = MatSetUp (*A); // "MatSetUp" sets up internal matrix data structure for the later use
  //disable error when preallocate not enough
  //check later
  ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);
  //ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(*A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE); CHKERRQ(ierr);
  CHKERRQ(ierr);
}

int matrix_mult::initialize()
{
  // initialize matrix
  setupMat();
  preAllocate();
  int ierr = MatSetUp (*A); // "MatSetUp" sets up internal matrix data structure for the later use
  //disable error when preallocate not enough
  //check later
  ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);
  //ierr = MatSetOption(*A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(*A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE); CHKERRQ(ierr);
  CHKERRQ(ierr);
}


m3dc1_matrix::m3dc1_matrix(int i, int s, FieldID f): id(i), scalar_type(s), fieldOrdering(f)
{
  mat_status = M3DC1_NOT_FIXED;
  A=new Mat;
}

int m3dc1_matrix::destroy()
{
  PetscErrorCode ierr = MatDestroy(A);
  CHKERRQ(ierr);    
}

m3dc1_matrix::~m3dc1_matrix()
{
  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<" deleting matrix "<<id<<"\n";
  destroy();
  delete A;
} 

int m3dc1_matrix::set_value(int row, int col, int operation, double real_val, double imag_val) //insertion/addition with global numbering
{
  if (mat_status == M3DC1_FIXED)
    return M3DC1_FAILURE;
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
  if (mat_status == M3DC1_FIXED)
    return M3DC1_FAILURE;
  PetscErrorCode ierr;
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  vector<PetscScalar> petscValues(rsize*csize);
  for(int i=0; i<rsize; i++)
  {
    //if(id==22)
      //std::cout<<std::endl<<"id "<<id<<" row "<<rows[i]<<std::endl;
    for(int j=0; j<csize; j++)
    {
      //if(id==22)
        //std::cout<<" colum "<<columns[j]<<" "<<values[i*csize+j]<<" ";
      if(scalar_type==M3DC1_REAL) petscValues.at(i*csize+j)=values[i*csize+j];
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
int matrix_solve :: add_blockvalues(int rbsize, int * rows, int cbsize, int * columns, double* values)
{
#if defined(DEBUG) || defined(PETSC_USE_COMPLEX)
  int bs;
  MatGetBlockSize(remoteA, &bs);
  vector<PetscScalar> petscValues(rbsize*cbsize*bs*bs);
  //std::cout<<PCU_Comm_Self()<<" bs "<<bs<<std::endl;
  //for(int i=0; i<rbsize; i++) std::cout<<" row "<<rows[i]<<" "<<std::endl;
  //for(int i=0; i<cbsize; i++) std::cout<<" columns "<<columns[i]<<" "<<std::endl;
  for(int i=0; i<rbsize*bs; i++)
  {
    for(int j=0; j<cbsize*bs; j++)
    {
      if(scalar_type==M3DC1_REAL) petscValues.at(i*cbsize*bs+j)=values[i*cbsize*bs+j];
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
      //std::cout<<PCU_Comm_Self()<<" "<<i*cbsize*bs+j<<" "<<petscValues.at(i*cbsize*bs+j)<<std::endl;
    }
  }
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, &petscValues[0], ADD_VALUES);
#else
  int ierr = MatSetValuesBlocked(remoteA,rbsize, rows, cbsize, columns, (PetscScalar*)values, ADD_VALUES);
#endif
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


// ***********************************
// 		matrix_mult
// ***********************************

int matrix_mult::multiply(FieldID in_field, FieldID out_field)
{
  if(!localMat)
  {
    Vec b, c;
    copyField2PetscVec(in_field, b, get_scalar_type());
    //std::cout<<" before mult "<<std::endl;
    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    int ierr = VecDuplicate(b, &c);CHKERRQ(ierr);
    MatMult(*A, b, c);
    copyPetscVec2Field(c, out_field, get_scalar_type());
    //std::cout<<" after mult "<<std::endl;
    //VecView(c, PETSC_VIEWER_STDOUT_WORLD);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
    return 0;
  }
  else
  {
    Vec b, c;
    int num_dof;
    m3dc1_field_getnumlocaldof(&in_field, &num_dof);
#ifdef DEBUG
    int num_dof2;
    m3dc1_field_getnumlocaldof(&out_field, &num_dof2);
    assert(num_dof==num_dof2);
#endif
    int bs;
    int ierr;
    MatGetBlockSize(*A, &bs);
    PetscScalar * array[2];
    m3dc1_field_getdataptr(&in_field, (double**)array);
#ifdef PETSC_USE_COMPLEX
    if(!get_scalar_type())
    {
      double * array_org = (double*)array[0];
      array[0] = new PetscScalar[num_dof];
      for(int i=0; i<num_dof; i++)
      {
        array[0][i]=array_org[i];
      }
    }
#endif
#ifdef OLDPETSC
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, num_dof, (PetscScalar*) array[0],&b); CHKERRQ(ierr);
#else
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[0],&b); CHKERRQ(ierr);
#endif
    m3dc1_field_getdataptr(&out_field, (double**)array+1);
#ifdef PETSC_USE_COMPLEX
    if(!get_scalar_type())
    {
      double * array_org = (double*)array[1];
      array[1] = new PetscScalar[num_dof];
      for(int i=0; i<num_dof; i++)
      {
        array[1][i]=array_org[i];
      }
    }
#endif
#ifdef OLDPETSC
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, num_dof, (PetscScalar*) array[1],&c); CHKERRQ(ierr);
#else
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, bs, num_dof, (PetscScalar*) array[1],&c); CHKERRQ(ierr);
#endif
    ierr=VecAssemblyBegin(b);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(b);  CHKERRQ(ierr);
    ierr=VecAssemblyBegin(c);  CHKERRQ(ierr);
    ierr=VecAssemblyEnd(c);  CHKERRQ(ierr);
    MatMult(*A, b, c);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
#ifdef PETSC_USE_COMPLEX
    if(!get_scalar_type())
    {
      double *datapt;
      m3dc1_field_getdataptr(&out_field, &datapt);
      for(int i=0; i<num_dof; i++)
        datapt[i]=std::real(array[1][i]); 
      delete []array[0];
      delete []array[1];
    }
#endif
    m3dc1_field_sum(&out_field);
  }
}
int matrix_mult::assemble()
{
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); 
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  set_status(M3DC1_FIXED);
}

// ***********************************
// 		matrix_solve
// ***********************************
matrix_solve::matrix_solve(int i, int s, FieldID f): m3dc1_matrix(i,s,f) 
{  
  ksp = new KSP;
  kspSet=0;
  initialize();
}

matrix_solve::~matrix_solve()
{
  if(kspSet)
    KSPDestroy(ksp);
  delete ksp;
  if(mat_status==M3DC1_NOT_FIXED && m3dc1_solver::instance()->assembleOption==0) MatDestroy(&remoteA);
}

int matrix_solve::assemble()
{
  PetscErrorCode ierr;
  double t1 = MPI_Wtime(), t2=t1;
 // if (!PCU_Comm_Self())
 // {
 //   std::cout<<"Before assemble"<<std::endl;
 //   printMemStat();
  //}
  if(!m3dc1_solver::instance()->assembleOption)
  {
    ierr = MatAssemblyBegin(remoteA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(remoteA, MAT_FINAL_ASSEMBLY);
    t2 = MPI_Wtime();
    //if (!PCU_Comm_Self()) std::cout<<"\t Assembly remoteA time "<<t2-t1<<std::endl;
    //ierr = MatView(remoteA, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    //pass remoteA to ownnering process
    int brgType = 2;
    if (m3dc1_mesh::instance()->mesh->getDimension()==3) brgType =3;

    int dofPerVar = 6;
    char field_name[256];
    int num_values, value_type, total_num_dof;
    m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
    dofPerVar=total_num_dof/num_values;
 
    int num_vtx, vertex_type=0;
    m3dc1_mesh_getnument (&vertex_type, &num_vtx);
    PetscInt firstRow, lastRowPlusOne;
    ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);

#ifdef USEPCUASSEMBLEMAT
    PCU_Comm_Begin();
    for(int inode=0; inode<num_vtx; inode++)
    {
      int owner=-1;
      m3dc1_ent_getownpartid (&vertex_type, &inode, &owner);
      if(owner==PCU_Comm_Self()) continue;

      apf::MeshEntity* ent =apf::getMdsEntity(m3dc1_mesh::instance()->mesh, vertex_type, inode);
      apf::MeshEntity* ownerEnt = get_ent_owncopy(m3dc1_mesh::instance()->mesh, ent);

      apf::Adjacent elements;
      getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);

      std::vector<apf::MeshEntity*> vecAdj;
      for (int i=0; i<elements.getSize(); ++i)
        vecAdj.push_back(elements[i]);
      vecAdj.push_back(ent);

      int numAdj = vecAdj.size();
      std::vector<int> globalId(numAdj);
      std::vector<int> localNodeId(numAdj);
      std::vector<int> columns(total_num_dof*numAdj);
      for(int i=0; i<numAdj; i++)
      {
        int local_id = get_ent_localid(m3dc1_mesh::instance()->mesh, vecAdj.at(i));
        localNodeId.at(i)=local_id;
        int start_global_dof_id, end_global_dof_id_plus_one;
        m3dc1_ent_getglobaldofid (&vertex_type, &local_id, &fieldOrdering, &start_global_dof_id, &end_global_dof_id_plus_one);
        globalId.at(i)=start_global_dof_id;
      }
      PCU_COMM_PACK(owner, ownerEnt);
      PCU_Comm_Pack(owner, &numAdj, sizeof(int));
      PCU_Comm_Pack(owner, &globalId.at(0), globalId.size()*sizeof(int));
      int offset=0;
      for(int i=0; i<numAdj; i++)
      {
        for(int j=0; j<total_num_dof; j++)
          columns.at(offset++)=localNodeId.at(i)*total_num_dof+j;
      }
      std::vector<PetscScalar> values (total_num_dof*numAdj*total_num_dof);
      ierr = MatGetValues(remoteA, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &values[0]);
      PCU_Comm_Pack(owner, &values[0], values.size()*sizeof(PetscScalar));
    }
    ierr = MatDestroy(&remoteA);
    CHKERRQ(ierr);
    PCU_Comm_Send();

    while (PCU_Comm_Listen())
    {
      while ( ! PCU_Comm_Unpacked())
      {
        apf::MeshEntity* ent;
        PCU_COMM_UNPACK(ent);
        int numAdj=-1;
        PCU_Comm_Unpack(&numAdj,sizeof(int));
        std::vector<int> globalId(numAdj);
        PCU_Comm_Unpack(&globalId[0],sizeof(int)*globalId.size());
        std::vector<PetscScalar> values (total_num_dof*numAdj*total_num_dof);
        PCU_Comm_Unpack(&values[0],sizeof(PetscScalar)*values.size());
        std::vector<int> columns(total_num_dof*numAdj);
        int offset=0;
        for(int i=0; i<numAdj; i++)
        {
          for(int j=0; j<total_num_dof; j++)
          {
            columns.at(offset++)=globalId.at(i)+j;
            //std::cout<<" columns "<<columns.at(offset-1)<<std::endl;
          }
        }
        assert (columns.at(total_num_dof*(numAdj-1))>=firstRow && *columns.rbegin()<lastRowPlusOne);
        ierr = MatSetValues(*A, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &values[0],ADD_VALUES);
      }
    }
#else
    std::map<int, std::vector<int> > idxSendBuff, idxRecvBuff;
    std::map<int, std::vector<PetscScalar> > valuesSendBuff, valuesRecvBuff;
    int blockMatSize = total_num_dof*total_num_dof;
    for(std::map<int, std::map<int, int> > ::iterator it = remoteNodeRow.begin(); it!=remoteNodeRow.end(); it++)
    {
      idxSendBuff[it->first].resize(it->second.size()+remoteNodeRowSize[it->first]);
      valuesSendBuff[it->first].resize(remoteNodeRowSize[it->first]*blockMatSize);
      int idxOffset=0;
      int valueOffset=0;
      for(std::map<int, int> ::iterator it2 =it->second.begin(); it2!=it->second.end();it2++)
      {
        idxSendBuff[it->first].at(idxOffset++)=it2->second;
        apf::MeshEntity* ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh, 0, it2->first);

        std::vector<apf::MeshEntity*> vecAdj;
        apf::Adjacent elements;
        getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);
        for (int i=0; i<elements.getSize(); ++i)
          vecAdj.push_back(elements[i]);

        vecAdj.push_back(ent);
        int numAdj = vecAdj.size();
        assert(numAdj==it2->second);
        std::vector<int> localNodeId(numAdj);
        std::vector<int> columns(total_num_dof*numAdj);
        for(int i=0; i<numAdj; i++)
        {
          int local_id = get_ent_localid(m3dc1_mesh::instance()->mesh, vecAdj.at(i));
          localNodeId.at(i)=local_id;
          int start_global_dof_id, end_global_dof_id_plus_one;
          m3dc1_ent_getglobaldofid (&vertex_type, &local_id, &fieldOrdering, &start_global_dof_id, &end_global_dof_id_plus_one);
          idxSendBuff[it->first].at(idxOffset++)=start_global_dof_id;
        }
        int offset=0;
        for(int i=0; i<numAdj; i++)
        {
          int startColumn = localNodeId.at(i)*total_num_dof;
          for(int j=0; j<total_num_dof; j++)
            columns.at(offset++)=startColumn+j;
        }
        ierr = MatGetValues(remoteA, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &valuesSendBuff[it->first].at(valueOffset));
        //for(int i=0; i<total_num_dof*numAdj; i++)
          //std::cout<<" get values indx "<<columns.at(i)<<std::endl;
        //for(int i=0; i<it2->second*blockMatSize; i++)
          //std::cout<<"values "<<i<<" "<<valuesSendBuff[it->first].at(valueOffset+i)<<std::endl;
        valueOffset+=it2->second*blockMatSize;
      }
      assert(idxOffset==idxSendBuff[it->first].size());
      assert(valueOffset==valuesSendBuff[it->first].size());
    }
    ierr = MatDestroy(&remoteA);

    //send and receive message size
    int sendTag=2020;
    MPI_Request my_request[256];
    MPI_Status my_status[256];
    int requestOffset=0;
    std::map<int, std::pair<int, int> > msgSendSize;
    std::map<int, std::pair<int, int> > msgRecvSize;
    for(std::map<int, int > :: iterator it = remoteNodeRowSize.begin(); it!=remoteNodeRowSize.end(); it++)
    {
      int destPid=it->first;
      msgSendSize[destPid].first=idxSendBuff[it->first].size();
      msgSendSize[destPid].second = valuesSendBuff[it->first].size();
      MPI_Isend(&(msgSendSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for(std::set<int> :: iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      MPI_Irecv(&(msgRecvSize[destPid]),sizeof(std::pair<int, int>),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request,my_status);
    //set up receive buff
    for(std::map<int, std::pair<int, int> > :: iterator it = msgRecvSize.begin(); it!= msgRecvSize.end(); it++)
    {
      idxRecvBuff[it->first].resize(it->second.first);
      valuesRecvBuff[it->first].resize(it->second.second); 
    }
    // now get data
    sendTag=9999;
    requestOffset=0;
    for(std::map<int, int > :: iterator it = remoteNodeRowSize. begin(); it!=remoteNodeRowSize.end(); it++)
    {
      int destPid=it->first;
      MPI_Isend(&(idxSendBuff[destPid].at(0)),idxSendBuff[destPid].size(),MPI_INT,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Isend(&(valuesSendBuff[destPid].at(0)),sizeof(PetscScalar)*valuesSendBuff[destPid].size(),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    for(std::set<int> :: iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      MPI_Irecv(&(idxRecvBuff[destPid].at(0)),idxRecvBuff[destPid].size(),MPI_INT,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
      MPI_Irecv(&(valuesRecvBuff[destPid].at(0)),sizeof(PetscScalar)*valuesRecvBuff[destPid].size(),MPI_BYTE,destPid,sendTag,MPI_COMM_WORLD,&(my_request[requestOffset++]));
    }
    assert(requestOffset<256);
    MPI_Waitall(requestOffset,my_request,my_status);
    //if (!PCU_Comm_Self())
    //{
    //  std::cout<<"Peak in scorec assemble"<<std::endl;
    //  printMemStat();
    //}

    for( std::map<int, std::vector<int> > :: iterator it =idxSendBuff.begin(); it!=idxSendBuff.end(); it++)
      std::vector<int>().swap(it->second);
    for( std::map<int, std::vector<PetscScalar> > :: iterator it =valuesSendBuff.begin(); it!=valuesSendBuff.end(); it++)
      std::vector<PetscScalar>().swap(it->second);
    valuesSendBuff.clear();
    idxSendBuff.clear();

    // now assemble the matrix
    for(std::set<int> :: iterator it = remotePidOwned.begin(); it!=remotePidOwned.end(); it++)
    {
      int destPid=*it;
      int valueOffset=0;
      int idxOffset=0;
      vector<int> & idx = idxRecvBuff[destPid];
      vector<PetscScalar> & values = valuesRecvBuff[destPid];
      int numValues=values.size();
      while(valueOffset<numValues)
      {
        int numAdj = idx.at(idxOffset++); 
        std::vector<int> columns(total_num_dof*numAdj);
        int offset=0;
        for(int i=0; i<numAdj; i++, idxOffset++)
        {
          for(int j=0; j<total_num_dof; j++)
          {
            columns.at(offset++)=idx.at(idxOffset)+j;
          }
        }
        assert (columns.at(total_num_dof*(numAdj-1))>=firstRow && *columns.rbegin()<lastRowPlusOne);
        ierr = MatSetValues(*A, total_num_dof, &columns.at(total_num_dof*(numAdj-1)), total_num_dof*numAdj, &columns[0], &values.at(valueOffset),ADD_VALUES);
        /*int start_row=total_num_dof*(numAdj-1);
        for(int i=0; i<total_num_dof; i++)
        {
          int row = columns.at(start_row++);
          for(int j=0; j<total_num_dof*numAdj; j++)
             MatSetValue(*A, row, columns.at(j), values.at(valueOffset++),ADD_VALUES);
        }*/
        valueOffset+=blockMatSize*numAdj;
      }
      std::vector<int>().swap(idxRecvBuff[destPid]);
      std::vector<PetscScalar>().swap(valuesRecvBuff[destPid]);
    }
    valuesRecvBuff.clear();
    idxRecvBuff.clear();
#endif
  }
  //double t3 = MPI_Wtime();
  //MPI_Barrier(MPI_COMM_WORLD);
  //double t5 = MPI_Wtime();
  //if (!PCU_Comm_Self()) std::cout<<"\tfill matrix "<<t5-t2<<std::endl;
  //if (!PCU_Comm_Self())
  //{
  //  std::cout<<"Before call MatAssemblyBegin"<<std::endl;
  //  printMemStat();
  //}

  //printInfo();
  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); 
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //if (!PCU_Comm_Self())
  //{
  //  std::cout<<"After call MatAssembly"<<std::endl;
  //  printMemStat();
  //}

  //double t6 = MPI_Wtime();
  //if (!PCU_Comm_Self()) std::cout<<"\tfinal assembly "<<t6-t5<<std::endl;
  //for(set<int>::iterator it =remoteProc.begin(); it!=remoteProc.end(); it++)
    //std::cout<<"proc "<<PCU_Comm_Self()<<" send to "<<*it<<std::endl; 
  mat_status=M3DC1_FIXED;
}
int m3dc1_matrix :: flushAssembly()
{
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(*A, MAT_FLUSH_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A, MAT_FLUSH_ASSEMBLY);
  CHKERRQ(ierr);
}
int matrix_solve:: set_bc( int row)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  MatSetValue(*A, row, row, 1.0, ADD_VALUES);
}

int matrix_solve:: set_row( int row, int numVals, int* columns, double * vals)
{
#ifdef DEBUG
  PetscInt firstRow, lastRowPlusOne;
  int ierr = MatGetOwnershipRange(*A, &firstRow, &lastRowPlusOne);
  assert (row>=firstRow && row<lastRowPlusOne);
#endif
  for(int i=0; i<numVals; i++)
  {
    if(get_scalar_type() == M3DC1_REAL) set_value(row, columns[i], 1, vals[i], 0);
    else set_value(row, columns[i], 1, vals[2*i], vals[2*i+1]); 
  }
}
int  m3dc1_matrix :: preAllocateParaMat()
{
  int bs=1;
  MatType type;
  MatGetType(*A, &type);

  int num_own_ent,num_own_dof=0, vertex_type=0;
  m3dc1_mesh_getnumownent (&vertex_type, &num_own_ent);
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
  int num_vtx=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

  int nnzStash=0;
  int brgType = 2;
  if (m3dc1_mesh::instance()->mesh->getDimension()==3) brgType =3;

  apf::MeshEntity* ent;
  for(int inode=0; inode<num_vtx; inode++)
  {
    ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh, vertex_type, inode);
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &fieldOrdering, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    if(start_global_dof_id<startDof || start_global_dof_id>=endDofPlusOne)
    {
      apf::Adjacent elements;
      getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);
      nnzStash+=dofPerEnt*dofPerEnt*(elements.getSize()+1);
      continue;
    }
    startIdx -= startDof;
    startIdx /=bs; 

    int adjNodeOwned, adjNodeGlb;
    m3dc1_mesh::instance()->mesh->getIntTag(ent, m3dc1_mesh::instance()->num_global_adj_node_tag, &adjNodeGlb);
m3dc1_mesh::instance()->mesh->getIntTag(ent, m3dc1_mesh::instance()->num_own_adj_node_tag, &adjNodeOwned);
    assert(adjNodeGlb>=adjNodeOwned);

    for(int i=0; i<numBlockNode; i++)
    {
      dnnz.at(startIdx+i)=(1+adjNodeOwned)*numBlockNode;
      onnz.at(startIdx+i)=(adjNodeGlb-adjNodeOwned)*numBlockNode;
    }
  }
  if (bs==1) 
    MatMPIAIJSetPreallocation(*A, 0, &dnnz[0], 0, &onnz[0]);
  else  
    MatMPIBAIJSetPreallocation(*A, bs, 0, &dnnz[0], 0, &onnz[0]);
  
  // print avg_num_non_zero, max_num_non_zero
  int max_dnnz_non_zero=dnnz[0], max_onnz_non_zero=onnz[0];
  int sum_dnnz_non_zero=0, sum_onnz_non_zero=0;
  for (int i=0; i<dnnz.size(); ++i)
  {
    sum_dnnz_non_zero+=dnnz[i];
    if (max_dnnz_non_zero<dnnz[i]) max_dnnz_non_zero=dnnz[i];
  }
  for (int i=0; i<onnz.size(); ++i)
  {
    sum_onnz_non_zero+=onnz[i];
    if (max_onnz_non_zero<onnz[i]) max_onnz_non_zero=onnz[i];
  }

  int avg_dnnz_non_zero=sum_dnnz_non_zero/dnnz.size();
  int avg_onnz_non_zero=sum_onnz_non_zero/onnz.size();

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] creating matrix "<<id<<": field_id "<<fieldOrdering<<", #dof "<<dofPerEnt
             <<", avg_dnnz_non_zero "<<avg_dnnz_non_zero<<", max_dnnz_non_zero "<<max_dnnz_non_zero<<"\n";
} 

int matrix_solve :: setUpRemoteAStruct()
{
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;

  int num_vtx, vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

  std::vector<int> nnz_remote(num_values*num_vtx);
  int brgType = 2;
  if (m3dc1_mesh::instance()->mesh->getDimension()==3) brgType =3;
  
  apf::MeshEntity* ent;
  for(int inode=0; inode<num_vtx; inode++)
  {
    ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh, vertex_type, inode);
    int owner=get_ent_ownpartid(m3dc1_mesh::instance()->mesh, ent);
    if (owner!=PCU_Comm_Self())
    {
      apf::Adjacent elements;
      getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);
      remoteNodeRow[owner][inode]=elements.getSize()+1;
      remoteNodeRowSize[owner]+=elements.getSize()+1;
      for(int i=0; i<num_values; i++)
        nnz_remote[inode*num_values+i]=(elements.getSize()+1)*num_values;
    }
    else 
    {
      apf::Copies remotes;
      m3dc1_mesh::instance()->mesh->getRemotes(ent,remotes);
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

}
int  m3dc1_matrix :: preAllocateSeqMat()
{
  int bs=1;
  MatType type;
  MatGetType(*A, &type);

  int num_vtx,num_dof=0, vertex_type=0;
  m3dc1_field_getnumlocaldof(&fieldOrdering, &num_dof);
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

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
  if (m3dc1_mesh::instance()->mesh->getDimension()==3) brgType = 3;

  apf::MeshEntity* ent;
  for(int inode=0; inode<num_vtx; inode++)
  {
    ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh, vertex_type, inode);
    int start_dof, end_dof_plus_one;
    m3dc1_ent_getlocaldofid (&vertex_type, &inode, &fieldOrdering, &start_dof, &end_dof_plus_one);
    int startIdx = start_dof;
    assert(startIdx<num_dof);

    apf::Adjacent elements;
    getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);

    int numAdj = elements.getSize();
    startIdx /=bs; 
    for(int i=0; i<numBlockNode; i++)
    {
      nnz.at(startIdx+i)=(1+numAdj)*numBlockNode;
    }
  }
  if (bs==1) 
    MatSeqAIJSetPreallocation(*A, 0, &nnz[0]);
  else  
    MatSeqBAIJSetPreallocation(*A, bs, 0, &nnz[0]);
} 

int m3dc1_matrix :: setupParaMat()
{
  int num_own_ent, vertex_type=0, num_own_dof;
  m3dc1_mesh_getnumownent (&vertex_type, &num_own_ent); 
  m3dc1_field_getnumowndof(&fieldOrdering, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;
  PetscInt mat_dim = num_own_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(MPI_COMM_WORLD, A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(*A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);

  ierr = MatSetType(*A, MATMPIAIJ);CHKERRQ(ierr);
}

int m3dc1_matrix :: setupSeqMat()
{
  int num_ent, vertex_type=0, num_dof;
  m3dc1_field_getnumlocaldof(&fieldOrdering, &num_dof);
  m3dc1_mesh_getnument (&vertex_type, &num_ent);
  int dofPerEnt=0;
  if (num_ent) dofPerEnt = num_dof/num_ent;

  PetscInt mat_dim = num_dof;

  // create matrix
  PetscErrorCode ierr = MatCreate(PETSC_COMM_SELF, A);
  CHKERRQ(ierr);
  // set matrix size
  ierr = MatSetSizes(*A, mat_dim, mat_dim, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);

  //use block mpi iij as default
  //ierr = MatSetType(*A, MATSEQBAIJ);CHKERRQ(ierr);
  //ierr = MatSetBlockSize(*A, dofPerEnt); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A);CHKERRQ(ierr);
}

int matrix_solve :: setupMat()
{
  setupParaMat();
}

int matrix_mult :: setupMat()
{
  if(localMat) setupSeqMat();
  else setupParaMat();
}

int matrix_solve :: preAllocate ()
{
  preAllocateParaMat();
}

int matrix_mult :: preAllocate ()
{
  if(localMat) preAllocateSeqMat();
  else preAllocateParaMat();
}

int copyField2PetscVec(FieldID field_id, Vec& petscVec, int scalar_type)
{
  int num_own_ent,num_own_dof=0, vertex_type=0;
  m3dc1_mesh_getnumownent (&vertex_type, &num_own_ent);
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  int ierr = VecCreateMPI(MPI_COMM_WORLD, num_own_dof, PETSC_DECIDE, &petscVec);
  CHKERRQ(ierr);
  VecAssemblyBegin(petscVec);

  int num_vtx=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  double dof_data[FIXSIZEBUFF];
  assert(sizeof(dof_data)>=dofPerEnt*2*sizeof(double));
  int nodeCounter=0;

  apf::MeshEntity* ent;
  for(int inode=0; inode<num_vtx; inode++)
  {
    ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh,vertex_type,inode);
    if (!is_ent_owned(m3dc1_mesh::instance()->mesh,ent)) continue;
      nodeCounter++;
    int num_dof;
    m3dc1_ent_getdofdata (&vertex_type, &inode, &field_id, &num_dof, dof_data);
    assert(num_dof*(1+scalar_type)<=sizeof(dof_data)/sizeof(double));
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx=0;
    for(int i=0; i<dofPerEnt; i++)
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
    //if(!PCU_Comm_Self()) 
    //for(int i=0; i<num_dof; i++) std::cout<<PCU_Comm_Self()<<" copy Field "<<field_id<<" to PetscVec inode "<<inode<<" "<<start_global_dof_id<<" "<<i<<" dof "<<dof_data[i]<<std::endl;

  }
  assert(nodeCounter==num_own_ent);
  ierr=VecAssemblyEnd(petscVec);
  CHKERRQ(ierr);
  return 0;
}

int copyPetscVec2Field(Vec& petscVec, FieldID field_id, int scalar_type)
{
  int num_own_ent,num_own_dof=0, vertex_type=0;
  m3dc1_mesh_getnumownent (&vertex_type, &num_own_ent);
  m3dc1_field_getnumowndof(&field_id, &num_own_dof);
  int dofPerEnt=0;
  if (num_own_ent) dofPerEnt = num_own_dof/num_own_ent;

  std::vector<PetscInt> ix(dofPerEnt);
  std::vector<PetscScalar> values(dofPerEnt);
  std::vector<double> dof_data(dofPerEnt*(1+scalar_type));
  int num_vtx=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

  int ierr;

  apf::MeshEntity* ent;
  for(int inode=0; inode<num_vtx; inode++)
  {
    ent = apf::getMdsEntity(m3dc1_mesh::instance()->mesh,vertex_type,inode);
    if (!is_ent_owned(m3dc1_mesh::instance()->mesh, ent)) continue;
    int start_global_dof_id, end_global_dof_id_plus_one;
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    
    for(int i=0; i<dofPerEnt; i++)
      ix.at(i)=startIdx+i;
    ierr=VecGetValues(petscVec, dofPerEnt, &ix[0], &values[0]); CHKERRQ(ierr);
    startIdx=0;
    for(int i=0; i<dofPerEnt; i++)
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
    //if(!PCU_Comm_Self()) 
    //for(int i=0; i<num_dof; i++) std::cout<<PCU_Comm_Self()<<" copy PetscVec2Field "<<field_id<<" inode "<<inode<<" "<<start_global_dof_id<<" "<<i<<" dof "<<dof_data[i]<<std::endl;
  }
  m3dc1_field_sync(&field_id);
  return 0;
}

int matrix_solve :: solve(FieldID field_id, int solver_type, double solver_tol)
{
  Vec x, b;
  copyField2PetscVec(field_id, b, get_scalar_type());
  int ierr = VecDuplicate(b, &x);CHKERRQ(ierr);
  //std::cout<<" before solve "<<std::endl;
  //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
  if (!kspSet) setKspType(solver_type, solver_tol);
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
}

int matrix_solve:: setKspType(int solver_type, double solver_tol)
{
  PetscErrorCode ierr;
  ierr = KSPCreate(MPI_COMM_WORLD, ksp);CHKERRQ(ierr);
#if defined(OLDPETSC) || defined(PETSC345)
  ierr = KSPSetOperators(*ksp, *A, *A , SAME_PRECONDITIONER/*, DIFFERENT_NONZERO_PATTERN */);CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(*ksp, *A, *A /*, SAME_PRECONDITIONER DIFFERENT_NONZERO_PATTERN*/);CHKERRQ(ierr);
#endif
  ierr = KSPSetTolerances(*ksp, .000001, solver_tol,
                          PETSC_DEFAULT, 1000);CHKERRQ(ierr);
  int num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(&fieldOrdering, field_name, &num_values, &value_type, &total_num_dof);
  assert(total_num_dof/num_values==C1TRIDOFNODE*(m3dc1_mesh::instance()->mesh->getDimension()-1));
  // if 2D problem use superlu
  if (m3dc1_mesh::instance()->mesh->getDimension()==2)
  {
#ifdef PETSC_USE_COMPLEX 
    ierr=KSPSetType(*ksp, KSPPREONLY);CHKERRQ(ierr);
    PC pc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr=PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);  CHKERRQ(ierr);
#else
    
  if (!solver_type) // iterative solver - Added by KK
  {
    ierr=KSPSetType(*ksp, KSPGMRES);CHKERRQ(ierr);
    PC pc;
    PC *subpc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCASM); CHKERRQ(ierr);

    PetscInt overlap = 1;    // width of subdomain overlap
    // Ref: http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex8.c.html
    ierr = PCASMSetOverlap(pc, overlap); CHKERRQ(ierr);
  }
  else    // Added by KK: If not iterative PETSc solver, use direct.
  {
    /* Commented out by KK
       if(1||num_values==1)
       {
    */
    ierr=KSPSetType(*ksp, KSPPREONLY);CHKERRQ(ierr);
    PC pc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr=PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU_DIST);  CHKERRQ(ierr);
  }
    /* Commented out by KK
       }
       else
       {

       // Commented out by KK
       ierr=KSPSetType(*ksp, KSPFGMRES);CHKERRQ(ierr);
       //int n;
       //double rnorm;
       //KSPMonitorDefault(*ksp,n, rnorm, NULL);
       PC pc;
       PC *subpc;
       ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
       ierr=PCSetType(pc,PCFIELDSPLIT); CHKERRQ(ierr);
       ierr =  PCFieldSplitSetBlockSize(pc,total_num_dof/num_values); CHKERRQ(ierr);
       for(int i=0; i<num_values; i++)
       {
       sprintf(field_name, "%dth",i);
       ierr =  PCFieldSplitSetFields(pc, field_name, 1, &i, &i);
       }
       if(num_values==2) PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR);
       ierr =  KSPSetUp(*ksp); CHKERRQ(ierr);
       KSP * subksp;
       int numSplit=-1;
       ierr = PCFieldSplitGetSubKSP(pc, &numSplit, &subksp);  CHKERRQ(ierr);
       assert(numSplit==num_values);
       for(int i=0; i<numSplit; i++)
       {
       ierr=KSPSetType(subksp[i], KSPPREONLY);CHKERRQ(ierr);
       PC pc;
       ierr=KSPGetPC(subksp[i], &pc); CHKERRQ(ierr);
       ierr=PCSetType(pc,PCLU); CHKERRQ(ierr); 
       ierr=PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU_DIST);  CHKERRQ(ierr);
       }
       PetscFree(subksp);
       }
    */
#endif    // PETSC_USE_COMPLEX
  }
  else // 3D mesh use bjacobi as default
    {
    /*
    ierr=KSPSetType(*ksp, KSPFGMRES);CHKERRQ(ierr);
    PC pc;
    ierr=KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
    ierr=PCBJacobiSetTotalBlocks(pc, m3dc1_model::instance()->num_plane, NULL); CHKERRQ(ierr);
    ierr=PCSetUp(pc);
    int n_local;
    KSP* subksps;
    ierr=PCBJacobiGetSubKSP(pc, &n_local, NULL, &subksps);
    assert(n_local<=1);
    for(int i=0; i<n_local; i++)
    {
      ierr=KSPSetType(subksps[i], KSPPREONLY);CHKERRQ(ierr);
    }
    ierr=KSPGetPC(subksps[0], &pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr=PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU_DIST);  CHKERRQ(ierr);*/
  } 
  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
  kspSet=1;
}

int m3dc1_matrix :: write (const char* file_name)
{
  PetscErrorCode ierr;
  PetscViewer lab;
  if(get_type()==0)
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
int m3dc1_matrix :: printInfo()
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
#endif //#ifdef M3DC1_PETSC
