/******************************************************************************
  (c) 2005-2016 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h" // debugging purpose
#include "m3dc1_matrix.h"
#include "m3dc1_field.h"
#include <PCU.h>
#include <pumi.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
int main(int argc, char * argv[])
{
  int result = 0;
  m3dc1_scorec_init(&argc,&argv);
  if (argc != 3 && !PCU_Comm_Self())
  {
    std::cout << "Usage: ./" << argv[0] << " model mesh" << std::endl;
    return M3DC1_FAILURE;
  }
  int num_plane = 2;
  assert(PCU_Comm_Peers() % num_plane == 0);
  m3dc1_model_setnumplane(&num_plane);
  m3dc1_model_load(argv[1]);
  m3dc1_model_print();
  m3dc1_mesh_load(argv[2]);
  printStats(m3dc1_mesh::instance()->get_mesh());
  int zero = 0;
  m3dc1_mesh_build3d(&zero, &zero, &zero);
    // set/get field dof values
  int num_dofs = 12;
  int dofs_per_blk[] = {12, 1};
  int blks_per_nd[] = {1, 12};
  int agg_scps[] = {m3dc1_field::LOCAL_AGGREGATION,
                    m3dc1_field::PLANE_AGGREGATION,
                    m3dc1_field::GLOBAL_AGGREGATION};
  const char * agg_str[] = {"lcl","pln","gbl"};
  std::stringstream fnm;
  int fld_id = 0;
  for(int agg_idx = 0; agg_idx < 3; ++agg_idx)
  {
    for(int blks_idx = 0; blks_idx < 2; ++blks_idx)
    {
      fnm << "agg_fld_" << dofs_per_blk[blks_idx]
          << "_" << blks_per_nd[blks_idx] << "_";
      fnm << agg_str[agg_idx];
      m3dc1_field_create(&fld_id,
                         fnm.str().c_str(),
                         &blks_per_nd[blks_idx],
                         &dofs_per_blk[blks_idx],
                         &agg_scps[agg_idx]);
      fnm.str("");
      fld_id++;
    }
  }
  int fld_cnt = fld_id;
  // get params to check field numberings
  int vrt_dim = 0;
  int lcl_vrt_cnt = 0;
  m3dc1_mesh_getnumownent(&vrt_dim,&lcl_vrt_cnt);
  int gbl_vrt_cnt = 0;
  m3dc1_mesh_getnumglobalent(&vrt_dim,&gbl_vrt_cnt);
  // develop first node dof ids for all cases
  int frst_gbl_vrt = 0;
  MPI_Exscan(&lcl_vrt_cnt,
             &frst_gbl_vrt,
             1,
             MPI_INTEGER,
             MPI_SUM,
             M3DC1_COMM_WORLD);
  // the first local vertex in the plane comm
  MPI_Comm pln_cm = m3dc1_model::instance()->getPlaneComm();
  int num_pln_vrts = 0;
  MPI_Allreduce(&lcl_vrt_cnt,
                &num_pln_vrts,
                1,
                MPI_INTEGER,
                MPI_SUM,
                pln_cm);
  int frst_pln_vrt = 0;
  MPI_Exscan(&lcl_vrt_cnt,
             &frst_pln_vrt,
             1,
             MPI_INTEGER,
             MPI_SUM,
             pln_cm);
  // the first plane vertex in the global comm
  int pln_rnk = -1;
  MPI_Comm_rank(pln_cm,&pln_rnk);
  int frst_gbl_pln_vrt = 0;
  int pln_cntrb = pln_rnk == 0 ? num_pln_vrts : 0;
  MPI_Scan(&pln_cntrb,
           &frst_gbl_pln_vrt,
           1,
           MPI_INTEGER,
           MPI_SUM,
           M3DC1_COMM_WORLD);
  frst_gbl_pln_vrt -= num_pln_vrts;
  // local aggregation
  int frst_gbl_dof = frst_gbl_vrt * num_dofs;
  int frst_pln_dof = frst_gbl_pln_vrt * num_dofs;
  int chk_lcl_nat[num_dofs];
  for(int ii = 0; ii < num_dofs; ++ii)
    chk_lcl_nat[ii] = frst_gbl_dof + ii;
  int chk_lcl_blk[num_dofs];
  for(int ii = 0; ii < num_dofs; ++ii)
    chk_lcl_blk[ii] = frst_gbl_dof + ii * lcl_vrt_cnt;
  // plane aggregation
  // the plane stride is half the gbl_vrt_cnt
  // because we hard-code that there are only two
  // planes and we know the meshes are identical
  int pln_strd = gbl_vrt_cnt >> 1;
  int chk_pln_nat[num_dofs];
  for(int ii = 0; ii < num_dofs; ++ii)
    chk_pln_nat[ii] = frst_gbl_dof + ii;
  int chk_pln_blk[num_dofs];
  for(int ii  = 0; ii < num_dofs; ++ii)
    chk_pln_blk[ii] = frst_pln_dof + (frst_pln_vrt + (ii * pln_strd));
  // global aggregation
  int chk_gbl_nat[num_dofs];
  for(int ii = 0; ii < num_dofs; ++ii)
    chk_gbl_nat[ii] = frst_gbl_dof + ii;
  int chk_gbl_blk[num_dofs];
  for(int ii = 0; ii < num_dofs; ++ii)
    chk_gbl_blk[ii] = (frst_gbl_vrt + (ii * gbl_vrt_cnt));
  // check numbering of all fields
  int * dof_id = new int[num_dofs];
  auto chk_dofs =
    [&](int * chk, int * ids) -> bool
    {
      for(int ii = 0; ii < num_dofs; ++ii)
        if(ids[ii] != chk[ii])
          return false;
      return true;
    };
  for(int fld_id = 0; fld_id < fld_cnt; ++fld_id)
  {
    int * chk = NULL;
    int fld_blks = 0;
    int dof_cnt = 0;
    int agg_scp = 0;
    m3dc1_field_getinfo(&fld_id,&fld_blks,&dof_cnt,&agg_scp);
    // find the first owned vertex since that is the first
    //  node we can reliably recreate the expected values for,
    //  ghost node's ids are determined by their owning part,
    //  making it much more difficult to determine the expected
    //  numbering for those nodes
    int vrt_id = 0;
    int own = 0;
    do
    {
      m3dc1_ent_isowner(&vrt_dim,&vrt_id,&own);
      vrt_id++;
    } while(!own);
    vrt_id--;
    m3dc1_ent_getglobaldofid(&vrt_dim,&vrt_id,&fld_id,&dof_id[0],&dof_cnt);
    if(agg_scp == m3dc1_field::LOCAL_AGGREGATION)
    {
      if(fld_blks == 1)
        chk = &chk_lcl_nat[0];
      else if(fld_blks == 12)
        chk = &chk_lcl_blk[0];
    }
    else if(agg_scp == m3dc1_field::PLANE_AGGREGATION)
    {
      if(fld_blks == 1)
        chk = &chk_pln_nat[0];
      else if(fld_blks == 12)
        chk = &chk_pln_blk[0];
    }
    else if(agg_scp == m3dc1_field::GLOBAL_AGGREGATION)
    {
      if(fld_blks == 1)
        chk = &chk_gbl_nat[0];
      else if(fld_blks == 12)
        chk = &chk_gbl_blk[0];
    }
    bool pass = chk_dofs(chk,&dof_id[0]);
    if(!pass)
    {
      int rnk = PCU_Comm_Self();
      std::cerr << "[" << rnk << "] field " << m3dc1_field_getname(&fld_id)
                << " has incorrect numbering!" << std::endl;
      result++;
    }
  }
  /*
  int sz = PCU_Comm_Peers();
  int rnk = PCU_Comm_Self();
  m3dc1_mesh * m3_msh = m3dc1_mesh::instance();
  apf::Mesh2 * msh = m3_msh->get_mesh();
  fld_id = 0;
  m3dc1_field * fld = m3_msh->get_field(fld_id);
  std::vector<int> blk_ids(1,-1);
  for(int rnd = 0; rnd < sz; ++rnd)
  {
    MPI_Barrier(M3DC1_COMM_WORLD);
    if(rnk == rnd)
    {
      std::cout << rnk << "\n";
      int lcl_vrts = 0;
      int vrt_dim = 0;
      m3dc1_mesh_getnument(&vrt_dim,&lcl_vrts);
      for(int vrt_idx = 0; vrt_idx < lcl_vrts; ++vrt_idx)
      {
        apf::MeshEntity * vtx = apf::getMdsEntity(msh,vrt_dim,vrt_idx);
        int vtx_lid = get_ent_localid(msh,vtx);
        int vtx_gid = get_ent_globalid(msh,vtx);
        get_block_ids(fld,false,0,vtx_lid,&blk_ids[0]);
        bool own_vtx = is_ent_original(vtx);
        std::cout << vtx_gid << " " << (own_vtx ? "own" : "gst") << " ";
        for(int blk = 0; blk < 1; ++blk)
          std::cout << blk_ids[blk] << " ";
        std::cout << std::endl;
      }
    }
  }
  */
  delete [] dof_id;
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return result;
}
