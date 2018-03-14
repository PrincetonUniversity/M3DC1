/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <apf.h>
#include <PCU.h>
#include <stdlib.h>
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apfMesh.h"

char* meshFile = 0;
char* modelFile = 0;
int num_plane=1;

void getConfig(int argc, char** argv)
{
  if ( argc < 4 ) 
  {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in-model> <in-mesh.smb> <#plane>\n", argv[0]);
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return;
  }
  modelFile = argv[1];
  meshFile = argv[2];
  if (atoi(argv[3])>1) num_plane=atoi(argv[3]);
}

int main(int argc, char** argv)
{
  m3dc1_scorec_init(&argc,&argv);

  getConfig(argc,argv);

  if (num_plane>1) // construct 3D
  {
    if (PCU_Comm_Peers()%num_plane!=0) 
    {
      if ( !PCU_Comm_Self() )
        printf("ERROR: invalid #planes\n");
      m3dc1_scorec_finalize();
      MPI_Finalize();
      return 1;
    }

    m3dc1_model_setnumplane (&num_plane);
  }


  m3dc1_model_load(modelFile);
  m3dc1_mesh_load(meshFile);
//  m3dc1_mesh::instance()->print(0);

  if (num_plane>1)
  {
    int zero=0;
    m3dc1_mesh_build3d(&zero, &zero, &zero);
  }

 // printStats(m3dc1_mesh::instance()->mesh);

//Is there any way to tell how many vertices per plane after partition is applied from just the mesh file without going to actually M3DC1 code? And how many vertices per partition as well, without counting the ghost vertices which belong to other rank?
  m3dc1_mesh::instance()->print(1);

  // finalize
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return 0;
}

