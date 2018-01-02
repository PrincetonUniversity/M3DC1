/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
/* Test access of fields on a ghosted mesh obtained via the workflow:
   APF -> omega_h -> ghost -> APF. This test assumes that the fields are
   real valued and that the mesh is 2d.
*/

#ifndef NOM3DC1
#define NOM3DC1
#endif

#include "m3dc1_ghost.h"
#include "m3dc1_mesh.h"
#include "m3dc1_model.h"
#include "apfOmega_h.h"
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "gmi_mesh.h"
#include "gmi.h"
#include "gmi_null.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "PCU.h"
#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <unistd.h>

using namespace std;

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc, &argv);
  m3dc1_scorec_init();
  int op = 0, scalar_type = M3DC1_REAL;
  int field_1 = 1, field_2 = 2, field_3 = 3;
  int num_values = 3, num_dofs = 6;
  int num_dofs_node = num_values * num_dofs;
  int num_vertex, num_own_vertex, vertex_dim = 0;
  int num_elem, elem_dim = 2;

  int num_plane = 1;

  if (argc < 2 && !PCU_Comm_Self()) {
    printf("Usage: ./ghost_example model mesh");
    return M3DC1_FAILURE;
  }

  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]);
  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);
  m3dc1_mesh_getnumownent (&vertex_dim, &num_own_vertex);
  m3dc1_mesh_getnument(&elem_dim, &num_elem);

  // Create and fill fields on APF mesh
  m3dc1_field_create(&field_1,
		     "field_1",
		     &num_values,
		     &scalar_type,
		     &num_dofs);
  m3dc1_field_create(&field_2,
		     "field_2",
		     &num_values,
		     &scalar_type,
		     &num_dofs);
  m3dc1_field_create(&field_3,
		     "field_3",
		     &num_values,
		     &scalar_type,
		     &num_dofs);

  // fill field_1
  printf("\n");
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    // 2D mesh, z-component = 0
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &field_1,
			 &num_dofs_node, &dofs.at(0));
  }
  m3dc1_field_printcompnorm(&field_1, "field_1 after set info");

  // fill field_2
  printf("\n");
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &field_2,
			 &num_dofs_node, &dofs.at(0));
  }
  m3dc1_field_printcompnorm(&field_2, "field_2 after set info");
  
  // fill field_3
  printf("\n");
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &field_3,
			 &num_dofs_node, &dofs.at(0));
  }
  m3dc1_field_printcompnorm(&field_3, "field_3 after set info");
 
  // Initialize ghosted mesh
  int nlayers = 1;
  m3dc1_ghost_load(&nlayers);

  // Check fields exist on ghosted mesh
  int exists_1, exists_2, exists_3;
  m3dc1_field_exist(&field_1, &exists_1);
  printf("\nField 1 on ghosted mesh exists: %d", exists_1);

  m3dc1_field_exist(&field_2, &exists_2);
  printf("\nField 2 on ghosted mesh exists: %d", exists_2);

  m3dc1_field_exist(&field_3, &exists_3);
  printf("\nField 3 on ghosted mesh exists: %d\n", exists_3);
  
  // Compute the norms of the fields on the ghosted_mesh
  printf("\n");
  m3dc1_field_printcompnorm(&field_1, "field_1 on ghosted mesh");
  m3dc1_field_printcompnorm(&field_2, "field_2 on ghosted mesh");
  m3dc1_field_printcompnorm(&field_3, "field_3 on ghosted mesh");
 
  // Scale field 1 on the ghosted by a constant and check norm again
  double factor = 2.0;
  printf("\nfield_1 scaled by %1.1f\n", factor);
  m3dc1_field_mult(&field_1, &factor, &scalar_type);
  m3dc1_field_printcompnorm(&field_1, "scaled field_1 on ghosted mesh");

  // Add field_2 and field_3 and store result on field_2
  m3dc1_field_add(&field_2, &field_3);
  printf("\nAdded field_3 to field_2.\n");
  m3dc1_field_printcompnorm(&field_2, "Updated field_2 on ghosted mesh");
  
  // Print ghosted simplices
  apf::MeshEntity* e;
  int ent_id;
  int isghost;
  int simplex_dim = 2;
  int num_simplices;
  int ghost_count = 0;

  m3dc1_mesh_getnument(&simplex_dim, &num_simplices);

  std::cout<< "\nNumber of simplices\n" << num_simplices;
  for(int ent_id=0; ent_id<num_simplices; ++ent_id) {
    isghost = 0;
    m3dc1_ent_isghost(&simplex_dim, &ent_id, &isghost);
    if (isghost != 0)
      ++ghost_count;
  }
  std::cout << "\nNumber of ghosted entities of ";
  std::cout << "dimension " << simplex_dim << ": " << ghost_count;
  std::cout << "\n";
  
  // Test field destruction
  m3dc1_field_delete(&field_1);

  // Test ghost mesh deletion
  m3dc1_ghost_delete();
  
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}





