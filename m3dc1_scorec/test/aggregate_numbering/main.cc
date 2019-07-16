/******************************************************************************
  (c) 2005-2019 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "m3dc1_scorec.h"
#include <PCU.h>
#include <petsc.h>
#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <unistd.h>
int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc, &argv);
  m3dc1_scorec_init();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  int op = 0, scalar_type = M3DC1_REAL;
  int field_1 = 1;
  int num_values = 1, num_dofs = 12;
  int num_dofs_node = num_values * num_dofs;
  int num_vertex, num_own_vertex, vertex_dim = 0;
  int num_elem, num_own_elem, elem_dim = 2;
  if (argc < 3 && !PCU_Comm_Self())
  {
    std::cout << "Usage: " << argv[0] << "[model_file] [mesh_file] [#planes]" << std::endl;
    return M3DC1_FAILURE;
  }
  int num_plane = atoi(argv[3]);
  if(num_plane > 1 && PCU_Comm_Peers() % num_plane == 0)
    m3dc1_model_setnumplane(&num_plane);
  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]);
  if(num_plane > 1)
  {
    int z = 0;
    m3dc1_mesh_build3d(&z,&z,&z);
  }
  int num_layer = 2;
  m3dc1_ghost_create(&num_layer);
  m3dc1_field_create(&field_1,
                     "field_1",
                     &num_values,
                     &scalar_type,
                     &num_dofs);
  int mat_id = 0;
  int mat_par_tp = 1;  // parallel matrix
  int agg_scp_pln = 1; // aggregate dofs per-plane
  m3dc1_matrix_create(&mat_id,
                      &mat_par_tp,
                      &scalar_type,
                      &field_1,
                      &num_dofs,
                      &agg_scp_pln);
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}
