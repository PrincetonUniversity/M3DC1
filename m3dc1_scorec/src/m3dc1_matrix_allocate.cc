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
  int num_nds[] = {(int)msh->get_mesh()->count(0), msh->get_own_count(0)};//{msh->num_local_ent[0],msh->num_own_ent[0]}; // assuming only verts hold nodes
  //int num_vrts = msh->get_mesh()->count(0);
  int is_par = (strcmp(tp,MATMPIBAIJ) == 0 || strcmp(tp,MATMPIAIJ) == 0);
  int is_lcl = !is_par;
  int dof_per_blk = fld->get_dof_per_value();
  int blk_per_nd = fld->get_num_value();
  int dof_per_nd = dof_per_blk * blk_per_nd;
  int num_lcl_blk = num_nds[is_par] * blk_per_nd;
  int frst_blk_row = 0;
  MPI_Exscan(&num_lcl_blk,&frst_blk_row,1,MPI_INTEGER,MPI_SUM,cms[is_par]);
  std::vector<int> dnnz(num_lcl_blk);
  std::vector<int> onnz(num_lcl_blk);
  int vrt_dim = 0;
  //FieldID fld_id = fld->get_id();
  int * lcl_dofs = new int[dof_per_nd];
  int * gbl_dofs = new int[dof_per_nd];
  int * dof_ids[] = {&lcl_dofs[0],&gbl_dofs[0]};
  int lcl_nd = 0;
  int blk = (bool)(dof_per_blk-1);
  int dnnz_max = 0;
  int onnz_max = 0;
  for(int nd = 0; nd < num_nds[0]; ++nd) //  have to iterate over all nodes, but we only process those we need
  {
    apf::MeshEntity * ent = apf::getMdsEntity(msh->get_mesh(),vrt_dim,nd);
    lcl_nd = is_ent_original(ent);
    if(!lcl_nd && is_par)
      continue;
    int dof_cnt = 0;
    int lcl_ent_id = get_ent_localid(msh->get_mesh(),ent);
    get_ent_localdofid(fld,lcl_ent_id,&(dof_ids[0][0]),&dof_cnt);
    get_ent_globaldofid(fld,lcl_ent_id,&(dof_ids[1][0]),&dof_cnt);
    int adj_own = 0;
    int adj_gbl = 0;
    msh->get_mesh()->getIntTag(ent, msh->own_bridge_adj_tag(), &adj_own);
    msh->get_mesh()->getIntTag(ent, msh->global_bridge_adj_tag(), &adj_gbl);
    assert(adj_gbl >= adj_own);
    for(int ii = 0; ii < dof_cnt; ii+=dof_per_blk)
    {
      int blk_row = (dof_ids[is_par][ii] / dof_per_blk) - frst_blk_row;
      assert(blk_row < num_lcl_blk && blk_row >= 0);
      dnnz[blk_row] = blk_per_nd * (blk ? 1 : dof_per_blk) * ((adj_own + 1) + (is_lcl * (adj_gbl - adj_own))); // if the matrix is local we only use dnnz but need the sum
      onnz[blk_row] = blk_per_nd * (blk ? 1 : dof_per_blk) * (adj_gbl - adj_own);
      dnnz_max = dnnz_max > dnnz[blk_row] ? dnnz_max : dnnz[blk_row];
      onnz_max = onnz_max > onnz[blk_row] ? onnz_max : onnz[blk_row];
    }
  }
  if(!blk && is_par)
    MatMPIAIJSetPreallocation(A,dnnz_max,&dnnz[0],onnz_max,&onnz[0]);
  else if(blk && is_par)
    MatMPIBAIJSetPreallocation(A,dof_per_blk,dnnz_max,&dnnz[0],onnz_max,&onnz[0]);
  else if(!blk && is_lcl)
    MatSeqAIJSetPreallocation(A,dnnz_max,&dnnz[0]);
  else // must be local and block
    MatSeqBAIJSetPreallocation(A,dof_per_blk,dnnz_max,&dnnz[0]);
  delete [] lcl_dofs;
  delete [] gbl_dofs;
}
