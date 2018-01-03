/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifdef M3DC1_TRILINOS
#include "m3dc1_ls.h"
#include "m3dc1_scorec.h"
#include "apf.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "apfMesh.h"
#include "apfMDS.h"
#include "PCU.h"
#include "m3dc1_mesh.h"
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

// ***********************************
// 		M3DC1_LINEAR SYSTEM
// ***********************************

m3dc1_ls* m3dc1_ls::_instance=NULL;
m3dc1_ls* m3dc1_ls::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_ls();
  return _instance;
}

m3dc1_ls::~m3dc1_ls()
{
  if (matrix_container!=NULL)
    matrix_container->clear();
  matrix_container=NULL;
  delete _instance;
}

void m3dc1_ls::add_matrix(int matrix_id, m3dc1_epetra* matrix)
{
  assert(matrix_container->find(matrix_id)==matrix_container->end());
  matrix_container->insert(std::map<int, m3dc1_epetra*>::value_type(matrix_id, matrix));
}

m3dc1_epetra* m3dc1_ls::get_matrix(int matrix_id)
{
  std::map<int, m3dc1_epetra*>::iterator mit = matrix_container->find(matrix_id);
  if (mit == matrix_container->end()) 
    return (m3dc1_epetra*)NULL;
  return mit->second;
}


//*******************************************************
apf::Numbering* get_owned_numbering()
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;
  int start=0;
  apf::Numbering* n = createNumbering(m,"m3dc1_owned_numbering",apf::getConstant(0),1);
  
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (get_ent_ownpartid(m,e)!=PCU_Comm_Self()) continue;
    apf::number(n,e,0,0, start++);
  }
  m->end(it);
  std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": #owned_nodes="<<start<<"\n";
  return n;
}



// ***********************************
// 		M3DC1_EPETRA
// ***********************************

Epetra_Map* createEpetraMap(global_ordinal_type num_dof, bool owned)
{
  global_ordinal_type num_node, global_id;
  if (owned) 
    num_node = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_own_ent[0]);
  else // overlap
    num_node = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_local_ent[0]);

  apf::DynamicArray<global_ordinal_type> dofIndices(num_node*num_dof);

  // loop over owned node
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  int i=0;
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    if (owned && get_ent_ownpartid(m3dc1_mesh::instance()->mesh,e)!=PCU_Comm_Self()) continue;
    global_id = get_ent_globalid(m3dc1_mesh::instance()->mesh, e);
    for (global_ordinal_type j=0; j < num_dof; ++j)
      dofIndices[i*num_dof + j] = global_id*num_dof + j;
    ++i;
  }
  m3dc1_mesh::instance()->mesh->end(it);
  return new Epetra_Map(-1,dofIndices.getSize(),&dofIndices[0],0,Epetra_MpiComm(MPI_COMM_WORLD));
}

m3dc1_epetra::m3dc1_epetra(int i, int t, int s, FieldID f_id): id(i), matrix_type(t), scalar_type(s), field_id(f_id), num_solver_iter(0)
{

  _field = (*(m3dc1_mesh::instance()->field_container))[f_id]->get_field();
  global_ordinal_type num_dof = static_cast<global_ordinal_type>(apf::countComponents(_field));
  _owned_map = createEpetraMap(num_dof,true);
  _overlap_map = createEpetraMap(num_dof,false);

  // compute #non_zero_in_row in global/local matrix
//  nge = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_local_ent[0]*num_dof);
  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0];
  int num_own_dof=num_own_ent*num_dof;

  std::vector<int> local_dnnz(num_own_dof), global_dnnz(num_own_dof);

  int startDof, endDofPlusOne, vertex_type=0;
  m3dc1_field_getowndofid (&field_id, &startDof, &endDofPlusOne);

  int brgType = m3dc1_mesh::instance()->mesh->getDimension();
  int start_global_dof_id, end_global_dof_id_plus_one;

  apf::MeshEntity* ent;
  for(int inode=0; inode<m3dc1_mesh::instance()->num_local_ent[0]; inode++)
  {
    ent = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, inode);
    m3dc1_ent_getglobaldofid (&vertex_type, &inode, &field_id, &start_global_dof_id, &end_global_dof_id_plus_one);
    int startIdx = start_global_dof_id;
    if(start_global_dof_id<startDof || start_global_dof_id>=endDofPlusOne)
      continue;

    startIdx -= startDof;

    int local_num_adj, global_num_adj;
    m3dc1_mesh::instance()->mesh->getIntTag(ent, m3dc1_mesh::instance()->num_global_adj_node_tag, &global_num_adj);
    apf::Adjacent elements;
    getBridgeAdjacent(m3dc1_mesh::instance()->mesh, ent, brgType, 0, elements);
    local_num_adj = elements.getSize();

    for(int i=0; i<num_dof; i++)
    {
      local_dnnz.at(startIdx+i)=(1+local_num_adj)*num_dof;
      global_dnnz.at(startIdx+i)=(1+global_num_adj)*num_dof;
    }
  }
  
  // print avg_num_non_zero, max_num_non_zero
  int max_local_dnnz_non_zero=local_dnnz[0];
  int max_global_dnnz_non_zero=global_dnnz[0];
  for (int i=0; i<local_dnnz.size(); ++i)
  {
    if (max_local_dnnz_non_zero<local_dnnz[i]) max_local_dnnz_non_zero=local_dnnz[i];
    if (max_global_dnnz_non_zero<global_dnnz[i]) max_global_dnnz_non_zero=global_dnnz[i];
  }
  nge = max_global_dnnz_non_zero;
  if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 INFO] m3dc1_epetra_create: ID "<<id<<", type "<<matrix_type
      <<", field "<<field_id<<", #dof "<<num_dof<<", #non_zero global "<<max_global_dnnz_non_zero
      <<", local "<<max_local_dnnz_non_zero<<"\n";

  epetra_mat = new Epetra_CrsMatrix(Copy,*_overlap_map,max_local_dnnz_non_zero,false);
} 

m3dc1_epetra::~m3dc1_epetra() {}

void m3dc1_epetra::destroy()
{
  delete epetra_mat;
  delete _overlap_map;
  delete _owned_map;

}
//=========================================================================
void write_matrix(Epetra_CrsMatrix* A, const char* matrix_filename, bool skip_zero, int start_index)
{
  int MyPID = PCU_Comm_Self();
  int NumProc = PCU_Comm_Peers();

  FILE * fp =fopen(matrix_filename, "w");

  int NumMyRows1 = A->NumMyRows();
  int MaxNumIndices = A->MaxNumEntries();

  int * Indices_int = 0;
  long long * Indices_LL = 0;

  if(A->RowMap().GlobalIndicesInt())
     Indices_int = new int[MaxNumIndices];
  else if(A->RowMap().GlobalIndicesLongLong()) 
     Indices_LL = new long long[MaxNumIndices];
  else
    if (!PCU_Comm_Self()) std::cout<<"[M3DC1_SCOREC FATAL] "<<__func__<<": Unable to determine source global index type\n";

  double * values  = new double[MaxNumIndices];
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  int NumIndices, j;
#endif
  int i, num_nonzero=0;
  if (!skip_zero)
    num_nonzero=A->NumMyNonzeros();
  else // count nonzero
  {
    for (i=0; i<NumMyRows1; i++) 
    {
      if(A->RowMap().GlobalIndicesInt()) 
      {
         int Row = (int) A->GRID64(i); // Get global row number
         A->ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_int);

         for (j = 0; j < NumIndices ; j++) 
           if (values[j]!=0.0) ++num_nonzero;
      }
      else if(A->RowMap().GlobalIndicesLongLong()) 
      {
         long long Row = A->GRID64(i); // Get global row number
         A->ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_LL);

         for (j = 0; j < NumIndices ; j++) 
           if (values[j]!=0.0) ++num_nonzero;
      }
    }
  }

  fprintf(fp,  "%d\t%d\t%d\n", A->NumMyRows(), A->NumMyCols(), num_nonzero);

  for (i=0; i<NumMyRows1; i++) 
  {
    if(A->RowMap().GlobalIndicesInt()) 
    {
       int Row = (int) A->GRID64(i); // Get global row number
       A->ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_int);

       for (j = 0; j < NumIndices ; j++) 
       {
          if (skip_zero && values[j]==0.0) continue;
          fprintf(fp, "%d\t%d\t%E\n", Row+start_index,Indices_int[j]+start_index,values[j]);
       }
    }
    else if(A->RowMap().GlobalIndicesLongLong()) 
    {
       long long Row = A->GRID64(i); // Get global row number
       A->ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_LL);

       for (j = 0; j < NumIndices ; j++) 
       {
         if (skip_zero && values[j]==0.0) continue;
         fprintf(fp, "%d\t%d\t%E\n", Row+start_index,Indices_LL[j]+start_index,values[j]);
       }
    }
  }

  if(A->RowMap().GlobalIndicesInt()) 
    delete [] Indices_int;
  else if(A->RowMap().GlobalIndicesLongLong()) 
    delete [] Indices_LL;
  delete [] values;
  fclose(fp);  
  return;
}

#endif
