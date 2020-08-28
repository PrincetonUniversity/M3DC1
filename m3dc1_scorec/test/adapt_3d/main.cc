/****************************************************************************** 

  (c) 2005-2020 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include "apfMesh.h"
#include <iostream>
#include "ma.h"
#include "PCU.h"
#include <mpi.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include <math.h>
#include <stdlib.h>

int main( int argc, char* argv[])
{
  std::cout << "3d Adapt Unit Test" << "\n";
  
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();

  if (argc<4 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes"<<endl;
    return M3DC1_FAILURE;
  }
  
  m3dc1_model_load(argv[1]);
  int  num_plane =  atoi(argv[3]);
 // int num_plane = 8;
  m3dc1_model_setnumplane (&num_plane);
  if (!PCU_Comm_Self()) std::cout<<" num plane "<<num_plane<<endl;
  m3dc1_mesh_load(argv[2]);
  // set/get field dof values
  int num_vertex;
  int vertex_dim = 0;
  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);

  int num_edge;
  int edge_dim = 1;
  m3dc1_mesh_getnument(&edge_dim, &num_edge);

  int num_elem;
  int elem_dim = 2;
  m3dc1_mesh_getnument(&elem_dim, &num_elem);

  // After getting the input model,mesh and number of planes, build 3D mesh
 
  // Set input for m3dc1_mesh_build3d()i
  int num_field = 0;	
  int field_id = 0;
  int dof_per_value = 0;

  // printStats (Mesh m): print global mesh entity counts per dimension
  printStats(m3dc1_mesh::instance()->mesh);			// Print Mesh stats before converting to 3D 
  m3dc1_mesh_build3d(&num_field, &field_id, &dof_per_value);	// Build 3d
  printStats(m3dc1_mesh::instance()->mesh);			// Print Mesh stats after converting to 3D

  // Setup the parameters needed to calculate the node error
  // First we need to find element_error_sum for get_node_error_from_elm()

  int NUM_TERM = 12;

  double edge_error[num_edge][NUM_TERM];
  double elem_error [num_elem][NUM_TERM];
  double elm_error_sum[2][num_elem];
  double elem_error_res[2][num_elem];
  
  
  // node_error output from here will go as input to find_sizefield()
  double node_error[2][num_vertex];
  double final_node_error[num_vertex];  
  int size = 1; 
  for (int j=0; j<num_vertex; ++j)
  {
  	for (int i=0; i<2; ++i)
  	{
  		node_error_3d_mesh (&elm_error_sum[i][j], &size, &node_error[i][j]);
		node_error[i][j] = sqrt(node_error[i][j]);
  	}
	final_node_error[j] = sqrt((node_error[1][j]*node_error[1][j])+(node_error[2][j]*node_error[2][j]));
  	double error_aimed = 0.005;		// Will come from C1 input file (Parameter: adapt_target_error)
  	int max_adapt_node = 10000;		// Will come from C1 input file (Parameter: iadapt_max_node)
  	int adapt_option = 1;			// Will come from C1 input file (Parameter: adapt_control)

  	find_sizefield(num_plane, final_node_error, &error_aimed, &max_adapt_node, &adapt_option );
   }


//  apf::Mesh2*  mesh =m3dc1_mesh::instance()->mesh;
  
//  apf::writeVtkFiles("after",mesh);
//  mesh->destroyNative();
//  apf::destroyMesh(mesh);
  MPI_Finalize();

  return 0;
}
