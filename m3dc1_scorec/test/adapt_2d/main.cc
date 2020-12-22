#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include "apfMesh.h"
#include <gmi_mesh.h>
#include <iostream>
#include "ma.h"
#include "pumi.h"
#include <mpi.h>
#include <lionPrint.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "petscksp.h"
#include <math.h>
#include <stdlib.h>


int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  m3dc1_scorec_init();
  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]);
  lion_set_verbosity(1);
  gmi_register_mesh();
  bool logInterpolation = true;
  
  apf::Mesh2*  mesh = m3dc1_mesh::instance()->mesh;


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
  apf::writeVtkFiles("adapt_before",mesh);


  SetSizeField sf(mesh);
  ma::Input* in = ma::configure(mesh, &sf,0,logInterpolation);

  std::cout << "cansnap() : " << mesh->canSnap() << "\n";
  in->shouldSnap = false;
  in->shouldTransferParametric = false;  
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;  
  in->goodQuality = 0.2;

  ma::adapt(in);
  mesh->verify();
  apf::writeVtkFiles("adapt_after",mesh);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return 0;
}



