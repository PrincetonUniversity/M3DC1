/******************************************************************************
  (c) 2005-2016 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h" // debugging purpose
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
  int num_blks = 12;
  int num_dofs = 1;
  int fld_id = 0;
  int agg_scps[] = {m3dc1_field::LOCAL_AGGREGATION,
                    m3dc1_field::PLANE_AGGREGATION,
                    m3dc1_field::GLOBAL_AGGREGATION};
  // num_dofs and num_blocks are flipped on this call to get a naturally-ordered field
  m3dc1_field_create(&fld_id,"nat",&num_dofs,&num_blks,&agg_scps[0]);
  ++fld_id;
  const char * agg_str[] = {"lcl","pln","gbl"};
  std::stringstream fnm;
  for(int agg_idx = 0; agg_idx < 3; ++agg_idx)
  {
    fnm << agg_str[agg_idx];
    m3dc1_field_create(&fld_id,
                       fnm.str().c_str(),
                       &num_blks,
                       &num_dofs,
                       &agg_scps[agg_idx]);
    fnm.str("");
    fld_id++;
  }
  int lcl_mat_ids[4] = {0, 1, 2, 3};
  int par_mat_ids[4] = {4, 5, 6, 7};
  int lcl_mat_tp = M3DC1_LOCAL_MAT;
  int par_mat_tp = M3DC1_PARALLEL_MAT;
  int scalar_tp = M3DC1_REAL;
  // this only works because the matrix ids are the same as the field ids for the first 4 fields
  int fld_ids[4] = {0, 1, 2, 3};
  for(int ii = 0; ii < 4; ++ii)
    m3dc1_matrix_create(&lcl_mat_ids[ii],&lcl_mat_tp,&scalar_tp,&fld_ids[ii]);
  for(int ii = 0; ii < 4; ++ii)
    m3dc1_matrix_create(&par_mat_ids[ii],&par_mat_tp,&scalar_tp,&fld_ids[ii]);
  /// assemble all matrices and check nonzero usage
  return result;
}
