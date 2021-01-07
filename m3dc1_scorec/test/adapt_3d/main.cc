/****************************************************************************** 

  (c) 2005-2021 Scientific Computation Research Center, 
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
#include "petscksp.h"
#include "PCU.h"

void find_isofield (ma::Mesh* m, int* field_id1, int* field_id2, double* dir)
{
  double average = ma::getAverageEdgeLength(m);

  double data = average/2.0;

  int size_data=1;
  for (int i=0; i<m3dc1_mesh::instance()->mesh->count(0); ++i)
  {
    m3dc1_node_setfield(&i, field_id1, &data, &size_data);
    m3dc1_node_setfield(&i, field_id2, &data, &size_data);

    dir[i*3+0] = 1.0;
    dir[i*3+1] = 0.0;
    dir[i*3+2] = 0.0;
  }
}

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

  if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  // Set input for m3dc1_mesh_build3d()
  int zero=0;

  // printStats (Mesh m): print global mesh entity counts per dimension
  if (num_plane>1)
    m3dc1_mesh_build3d(&zero, &zero, &zero);

  int shouldSnap=0;
  if (argc>4) 
    shouldSnap=atoi(argv[4]);
  else
    if (!PCU_Comm_Self()) std::cout <<"Missing Input Args argv[4]: shouldSnap, argv[5]: shouldRunPreZoltan, "
         <<"argv[6]: shouldRunPostZoltan, argv[7]: shouldefineLayer\n"
         <<"Default value 0 is used for all\n"; 

  int shouldRunPreZoltan=0;
  int shouldRunPostZoltan=0;
  int shouldRefineLayer=0;

  if (argc>5) shouldRunPreZoltan=atoi(argv[5]);
  if (argc>6) shouldRunPostZoltan=atoi(argv[6]);
  if (argc>7) shouldRefineLayer=atoi(argv[7]);
    
//  apf::writeVtkFiles("before-adapt",m3dc1_mesh::instance()->mesh);
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;

  int fid_size1=1;
  int fid_size2=2;
  int scalar_type=0;
  int num_value=1;

  m3dc1_field_create (&fid_size1, "size 1", &num_value, &scalar_type, &num_value);
  m3dc1_field_create (&fid_size2, "size 2", &num_value, &scalar_type, &num_value);
   
  double* dir = new double[m->count(0)*3];

  find_isofield (m, &fid_size1, &fid_size2, dir);

  int maximumIterations=5;
  double goodQuality =0.2;
 
  if (!PCU_Comm_Self()) std::cout << "start adaptation\n";

  m3dc1_mesh_adapt(&fid_size1, &fid_size2, dir, &shouldSnap, &shouldRunPreZoltan,
      &shouldRunPostZoltan, &shouldRefineLayer, &maximumIterations, &goodQuality);

  if (!PCU_Comm_Self()) std::cout << "adaptation completed\n";

  delete [] dir;

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
