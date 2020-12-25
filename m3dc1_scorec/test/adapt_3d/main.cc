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
#include "PCU.h"

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

  // Set input for m3dc1_mesh_build3d()
  int num_field = 0;	
  int field_id = 0;
  int dof_per_value = 0;

  // printStats (Mesh m): print global mesh entity counts per dimension
  m3dc1_mesh_build3d(&num_field, &field_id, &dof_per_value);	// Build 3d
//  apf::writeVtkFiles("3d",m3dc1_mesh::instance()->mesh); 
//  pumi_mesh_print(m3dc1_mesh::instance()->mesh);

//  apf::writeVtkFiles("before-adapt",m3dc1_mesh::instance()->mesh);

  int shouldSnap=1;
  int shouldRunPreZoltan=1;
  int shouldRunPostZoltan=1;
  int shouldRefineLayer=0;
  int maximumIterations=5;
  double goodQuality =0.2;
 
  if (!PCU_Comm_Self()) std::cout << "start adaptation\n";

  m3dc1_mesh_adapt(&shouldSnap, &shouldRunPreZoltan,
      &shouldRunPostZoltan, &shouldRefineLayer, &maximumIterations, &goodQuality);

  if (!PCU_Comm_Self()) std::cout << "adaptation completed\n";

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
