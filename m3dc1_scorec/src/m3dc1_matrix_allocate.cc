#include "m3dc1_matrix_allocate.h"
#include "m3dc1_matrix.h"
#include "m3dc1_field.h"
#include <apfMDS.h>
#include <petscmat.h>
#include <algorithm>
#include <numeric>
#include <cassert>
void allocateMatrix(Mat A, m3dc1_matrix * m3_mat, m3dc1_mesh * m3_msh, m3dc1_field * fld)
{
  MPI_Comm cm;
  PetscObjectGetComm((PetscObject)A,&cm);
  int rnk = -1;
  int sz = 0;
  MPI_Comm_rank(cm,&rnk);
  MPI_Comm_size(cm,&sz);
  bool lcl = (sz == 1);
  MatType tp;
  MatGetType(A,&tp);
  apf::Mesh2 * msh = m3_msh->get_mesh();
  int dim = msh->getDimension();
  int num_ents = (int)msh->count(dim);
  int num_vrts = (int)msh->count(0);
  int dof_per_blk = fld->getDofsPerBlock();
  int blk_per_nd = fld->getBlocksPerNode();
  int num_gst_blk = 0;
  // the mesh entitty ownership has NOTHING to do with row ownership anymore..
  // so we have to explicitly count the owned and unowned rows
  std::vector<int> blk_ids(blk_per_nd,-1);
  for(int vid = 0; vid < num_vrts; ++vid)
  {
    apf::MeshEntity * vtx = apf::getMdsEntity(msh,0,vid);
    int vtx_lcl_id = get_ent_localid(msh,vtx);
    get_block_ids(fld,lcl,0,vtx_lcl_id,&blk_ids[0]);
    for(int blk_id = 0; blk_id < blk_per_nd; ++blk_id)
    {
      int i_rw = blk_ids[blk_id];
      int own_i_rw = m3_mat->willOwn(i_rw);
      if(!own_i_rw)
        num_gst_blk++;
    }
  }
  int frst_rw = m3_mat->calcFirstRow();
  int last_rw = m3_mat->calcLastRowP1();
  int num_rws = last_rw - frst_rw;
  int rws = 0;
  MPI_Allreduce(&num_rws,&rws,1,MPI_INTEGER,MPI_SUM,cm);
  // processing info for owned and ghost rows
  std::vector<std::vector<int> > prd_own(num_rws,std::vector<int>(rws+1,0));
  std::vector<std::vector<int> > prd_gst(num_gst_blk,std::vector<int>(rws+1,0));
  std::vector<int> dnnz(num_rws);
  std::vector<int> onnz(num_rws);
  // index intro prd_gst, so the gst idx is the first time we see gst rw
  std::vector<int> gbl_2_gst(rws,-1);
  int gst_idx = 0;
  std::vector<int> blk_ids_j(blk_per_nd);
  for(int vid = 0; vid < num_vrts; ++vid)
  {
    apf::MeshEntity * vtx = apf::getMdsEntity(msh,0,vid);
    int vtx_lcl_id = get_ent_localid(msh,vtx);
    get_block_ids(fld,lcl,0,vtx_lcl_id,&blk_ids[0]);
    for(int blk_id = 0; blk_id < blk_per_nd; ++blk_id)
    {
      int i_rw = blk_ids[blk_id];
      int own_i_rw = m3_mat->willOwn(i_rw);
      if(!own_i_rw)
      {
        gbl_2_gst[i_rw] = gst_idx;
        prd_gst[gst_idx][rws] = i_rw;
        gst_idx++;
      }
    }
  }
  gst_idx = 0;
  //std::cout << rnk << "\n";
  //
  for(int eid = 0; eid < num_ents; ++eid)
  {
    apf::MeshEntity * ent = apf::getMdsEntity(msh,dim,eid);
    //int geid = get_ent_globalid(msh,ent);
    bool ownd = is_ent_original(ent);
    //std::cout << "  " << geid << " , " << (ownd ? "own" : "gst" ) << " , ";
    if(lcl || ownd)
    {
      int adj_vrts = apf::Mesh::adjacentCount[msh->getType(ent)][0];
      std::vector<apf::MeshEntity*> vrts(adj_vrts,NULL);
      msh->getDownward(ent,0,&vrts[0]);
      for(int vrt_id = 0; vrt_id < adj_vrts; ++vrt_id)
      {
        int vtx_lcl_id = get_ent_localid(msh,vrts[vrt_id]);
        get_block_ids(fld,lcl,0,vtx_lcl_id,&blk_ids[0]);
        for(int blk_id = 0; blk_id < blk_per_nd; ++blk_id)
        {
          int i_rw = blk_ids[blk_id];
          int own_i_rw = m3_mat->willOwn(i_rw);
          //std::cout << i_rw << " , " << (own_i_rw == 1 ? "own" : "gst")
          //<< (vrt_id == adj_vrts-1 && blk_id == blk_per_nd-1 ? "" : " , ");
          //if(lcl || ownd)
          {
            for(int vrt_id_j = 0; vrt_id_j < adj_vrts; ++vrt_id_j)
            {
              int vtx_lcl_id_j = get_ent_localid(msh,vrts[vrt_id_j]);
              get_block_ids(fld,lcl,0,vtx_lcl_id_j,&blk_ids_j[0]);
              for(int blk_id_j = 0; blk_id_j < blk_per_nd; ++blk_id_j)
              {
                int j_rw = blk_ids_j[blk_id_j];
                //int own_j_rw = m3_mat->willOwn(j_rw);
                if(own_i_rw)
                {
                  int adj_i_rw = i_rw - frst_rw;
                  int & seen = prd_own[adj_i_rw][j_rw];
                  if(seen == 0)
                    seen = 1;
                }
                else
                {
                  int gst_idx = gbl_2_gst[i_rw];
                  int & seen = prd_gst[gst_idx][j_rw];
                  if(seen == 0)
                    seen = 1;
                }
              }
            }
          }
        }
      }
      //std::cout << std::endl;
    }
  }
// send non-local row info
  std::vector<int> snds(sz,0);
  for(int vid = 0; vid < num_vrts; ++vid)
  {
    apf::MeshEntity * vtx = apf::getMdsEntity(msh,0,vid);
    int vtx_lcl_id = get_ent_localid(msh,vtx);
    get_block_ids(fld,lcl,0,vtx_lcl_id,&blk_ids[0]);
    for(int blk_id = 0; blk_id < blk_per_nd; ++blk_id)
    {
      int i_rw = blk_ids[blk_id];
      int own_i_rw = m3_mat->willOwn(i_rw);
      if(!own_i_rw)
      {
        int ownr = m3_mat->whoOwns(i_rw);
        int gst_idx = gbl_2_gst[i_rw];
        MPI_Request rqst;
        MPI_Isend(&prd_gst[gst_idx][0],rws+1,MPI_INTEGER,ownr,0,cm,&rqst);
        snds[ownr]++;
      }
    }
  }
  std::vector<int> rcvs(sz,0);
  MPI_Allreduce(&snds[0],&rcvs[0],sz,MPI_INTEGER,MPI_SUM,cm);
  std::vector<MPI_Request> rqsts(rcvs[rnk]);
  std::vector<std::vector<int> > recv_rws(rcvs[rnk],std::vector<int>(rws+1,0));
  auto or_op = [] (int a, int b) -> int { return (a || b ? 1 : 0); };
  for(int recv_idx = 0; recv_idx < rcvs[rnk]; ++recv_idx)
    MPI_Irecv(&recv_rws[recv_idx][0],rws+1,MPI_INTEGER,MPI_ANY_SOURCE,0,cm,&rqsts[recv_idx]);
  int recvd = 0;
  int rqst_idx = 0;
  while(recvd != rcvs[rnk])
  {
    int flg = 0;
    MPI_Status sts;
    MPI_Test(&rqsts[rqst_idx],&flg,&sts);
    if(flg)
    {
      std::vector<int> & rw = recv_rws[rqst_idx];
      int gbl_rw = rw[rws];
      assert(m3_mat->willOwn(gbl_rw));
      int adj_lcl_rw = gbl_rw - frst_rw;
      // merge the processed rows
      std::transform(rw.begin(),rw.end()-1,
                     prd_own[adj_lcl_rw].begin(),
                     prd_own[adj_lcl_rw].begin(),
                     or_op);
      recvd++;
    }
    rqst_idx = (rqst_idx + 1) % rcvs[rnk];
  }
  for(int rw = 0; rw < num_rws; ++rw)
  {
    onnz[rw] += std::accumulate(prd_own[rw].begin(),
                                prd_own[rw].begin()+frst_rw,0);
    dnnz[rw] += std::accumulate(prd_own[rw].begin()+frst_rw,
                                prd_own[rw].begin()+frst_rw+num_rws,0);
    // don't add the last element (col) as it indicates global row_id
    onnz[rw] += std::accumulate(prd_own[rw].begin()+frst_rw+num_rws,
                                prd_own[rw].end()-1,0);

  }
  int dnnz_max = *std::max_element(dnnz.begin(),dnnz.end());
  int onnz_max = *std::max_element(onnz.begin(),onnz.end());
  bool is_par = sz > 1;
  int blk_sz = 0;
  MatGetBlockSize(A,&blk_sz);
  bool blk = blk_sz > 1;
  if(!blk && is_par)
    MatMPIAIJSetPreallocation(A,dnnz_max,&dnnz[0],onnz_max,&onnz[0]);
  else if(blk && is_par)
    MatMPIBAIJSetPreallocation(A,dof_per_blk,dnnz_max,&dnnz[0],onnz_max,&onnz[0]);
  else if(!blk && !is_par)
    MatSeqAIJSetPreallocation(A,dnnz_max,&dnnz[0]);
  else // must be local and block
    MatSeqBAIJSetPreallocation(A,dof_per_blk,dnnz_max,&dnnz[0]);
}
