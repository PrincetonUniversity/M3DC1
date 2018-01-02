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
  int num_elem, num_own_elem, elem_dim = 2;

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
  m3dc1_mesh_getnumownent (&elem_dim, &num_own_elem);

  // Create and fill fields on APF mesh
  m3dc1_field_create(&field_1,
		     "field_1",
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
 
  // Initialize ghosted mesh
  int nlayers = 2;
  m3dc1_ghost_load(&nlayers);
  int simplex_dim = 2;
  apf::Mesh2* mesh = m3dc1_ghost::instance()->mesh;
  int num_simplices = mesh->count(simplex_dim);
  std::cout << "\n[Proc " << PCU_Comm_Self();
  std::cout << "] mesh->count(simplex_dim) = " << num_simplices;
  std::cout << "\n[Proc " << PCU_Comm_Self();
  std::cout << "] m3dc1_mesh_getnument()   = " << num_elem;
  std::cout << "\n[Proc " << PCU_Comm_Self();
  std::cout << "] m3dc1_mesh_getnumownent()   = " << num_own_elem;
  std::cout << "\n";

  // Retrieve global element id
  for(int ielem=0; ielem < num_simplices; ielem++) {
    int gielem;
    m3dc1_ent_getglobalid(&ielem, &simplex_dim, &gielem);
    std::cout << "\n[Proc " << PCU_Comm_Self();
    std::cout << "] Element: local index = " << ielem;
    std::cout << "]  global index = " << gielem;
  }
  
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}




