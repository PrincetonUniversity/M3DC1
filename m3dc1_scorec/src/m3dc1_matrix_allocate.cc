#include "m3dc1_matrix_allocate.h"
#include "m3dc1_field.h"
#include <apfMDS.h>
#include <petscmat.h>
#include <algorithm>
#include <cassert>
// this works for standard dof ordering and rank-aggregate ordering (might not work on group-aggregate or global-aggregate)
//  need to operate on the fld->numbering instead...
void allocateMatrix(Mat A, m3dc1_mesh * msh, m3dc1_field * fld)
{
  MatType tp;
  MatGetType(A,&tp);
  MPI_Comm cms[] = {PETSC_COMM_SELF,PETSC_COMM_WORLD};
  int num_nds[] = {msh->num_local_ent[0],msh->num_own_ent[0]}; // assuming only verts hold nodes
  int num_vrts = msh->mesh->count(0);
  int is_par = (strcmp(tp,MATMPIBAIJ) == 0 || strcmp(tp,MATMPIAIJ) == 0);
  int is_lcl = !is_par;
  int bs = fld->get_dof_per_value();
  int dof_per_nd = fld->get_num_value() * bs;
  int blk_per_nd = dof_per_nd / bs;
  int num_own_dof = num_nds[is_par] * dof_per_nd;
  int frst_blk_row = 0;
  MPI_Exscan(&num_own_dof,&frst_blk_row,1,MPI_INTEGER,MPI_SUM,cms[is_par]);
  frst_blk_row /= dof_per_nd;
  int blk_row_cnt = num_own_dof / bs;
  std::vector<int> dnnz(blk_row_cnt);
  std::vector<int> onnz(blk_row_cnt);
  int vrt_dim = 0;
  //FieldID fld_id = fld->get_id();
  int * dofs = new int[dof_per_nd];
  int lcl_nd = 0;
  for(int nd = 0; nd < num_nds[0]; ++nd) //  have to iterate over all nodes, but we only process those we need
  {
    //DBG(memset(&gbl_dofs[0],0,sizeof(int)*dof_per_nd));
    //m3dc1_ent_getlocaldofid(&vrt_dim,&nd,&fld_id,&dofs[0],&dof_cnt);
    apf::MeshEntity * ent = apf::getMdsEntity(msh->mesh,vrt_dim,nd);
    int dof_cnt = 0;
    if(is_par)
    {
      int gbl_ent_id = get_ent_globalid(msh->mesh,ent);
      get_ent_globaldofid(fld,gbl_ent_id,&dofs[0],&dof_cnt);
    }
    else
    {
      int lcl_ent_id = get_ent_localid(msh->mesh,ent);
      get_ent_localdofid(fld,lcl_ent_id,&dofs[0],&dof_cnt);
    }
    lcl_nd = is_ent_original(msh->mesh,ent);
    if(!lcl_nd && is_par)
      continue;
    int adj_own = 0;
    int adj_gbl = 0;
    msh->mesh->getIntTag(ent, msh->num_own_adj_node_tag, &adj_own);
    msh->mesh->getIntTag(ent, msh->num_global_adj_node_tag, &adj_gbl);
    assert(adj_gbl >= adj_own);
    for(int ii = 0; ii < dof_cnt; ii+=bs)
    {
      int blk_row = (dofs[ii] / bs) - frst_blk_row;
      assert(blk_row < blk_row_cnt);
      dnnz[blk_row] = blk_per_nd * ((adj_own + 1) + (is_lcl * (adj_gbl - adj_own))); // if the matrix is local we only use dnnz but need the sum
      onnz[blk_row] = blk_per_nd * (adj_gbl - adj_own);
    }
  }
  int blk = (bool)(bs-1);
  if(!blk && is_par)
    MatMPIAIJSetPreallocation(A,0,&dnnz[0],0,&onnz[0]);
  else if(blk && is_par)
    MatMPIBAIJSetPreallocation(A,bs,0,&dnnz[0],0,&onnz[0]);
  else if(!blk && is_lcl)
    MatSeqAIJSetPreallocation(A,0,&dnnz[0]);
  else // must be local and block
    MatSeqBAIJSetPreallocation(A,bs,0,&dnnz[0]);
  delete [] dofs;
}
