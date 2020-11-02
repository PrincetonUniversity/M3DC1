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
#include "pumi.h"
#include <mpi.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "petscksp.h"
#include <math.h>
#include <stdlib.h>

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

  if (argc<3 & !pumi_rank())
  {
    cout<<"Usage: ./main  model mesh #planes"<<endl;
    return M3DC1_FAILURE;
  }
  
  int num_plane=pumi_size();
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && pumi_size()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);
  }

  if (m3dc1_model_load(argv[1])) // model loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  // m3dc1_model_print();

  if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  // set/get field dof values
  int num_vertex, nvertex=m3dc1_mesh::instance()->mesh->count(0);
  int num_edge, nedge=m3dc1_mesh::instance()->mesh->count(1);
  int num_elem, nelem=m3dc1_mesh::instance()->mesh->count(2);

  MPI_Allreduce(&nvertex, &num_vertex, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  MPI_Allreduce(&nedge, &num_edge, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  MPI_Allreduce(&nelem, &num_elem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 

  // After getting the input model,mesh and number of planes, build 3D mesh
 
  // Set input for m3dc1_mesh_build3d()
  int num_field = 0;	
  int field_id = 0;
  int dof_per_value = 0;

  // printStats (Mesh m): print global mesh entity counts per dimension
  m3dc1_mesh_build3d(&num_field, &field_id, &dof_per_value);	// Build 3d
  apf::writeVtkFiles("3d",m3dc1_mesh::instance()->mesh); 
  pumi_mesh_print(m3dc1_mesh::instance()->mesh);

  m3dc1_mesh::instance()->remove_wedges();
  pumi_mesh_print(m3dc1_mesh::instance()->mesh);
  apf::writeVtkFiles("2.5d",m3dc1_mesh::instance()->mesh);
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
   }
  double error_aimed = 0.005;             // Will come from C1 input file (Parameter: adapt_target_error)
  int max_adapt_node = 10000;             // Will come from C1 input file (Parameter: iadapt_max_node)
  int adapt_option = 1;                   // Will come from C1 input file (Parameter: adapt_control)

  find_sizefield(final_node_error, &error_aimed, &max_adapt_node, &adapt_option );

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
