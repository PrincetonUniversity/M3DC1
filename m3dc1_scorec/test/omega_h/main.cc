/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
/* Load an APF mesh using a mesh file and convert to omega_h file;
   ghost omega_h file, and convert back to APF.
*/

#ifndef NOM3DC1
#define NOM3DC1
#endif

#include "m3dc1_ghost.h"
#include "apfOmega_h.h"
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "gmi_mesh.h"
#include "gmi_null.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "PCU.h"
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in.smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  // Load mesh using null model
  gmi_register_null();
  apf::Mesh2* mesh = apf::loadMdsMesh(".null", argv[1]);
  mesh->verify();

  // APF -> omega_h -> APF
  osh_t om = osh::fromAPF(mesh);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  mesh = apf::makeEmptyMdsMesh(gmi_load(".null"), osh_dim(om), false);
  osh_ghost(om, 1);
  osh::toAPF(om, mesh);

  // Clean up and finalize
  osh_free(om);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}



