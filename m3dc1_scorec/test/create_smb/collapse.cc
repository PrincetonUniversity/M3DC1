#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <parma.h>
#include <pcu_util.h>
#include <gmi_null.h>
#include <cstdlib>
#include <iostream>

namespace {

  const char* modelFile = 0;
  const char* meshFile = 0;
  const char* outFile = 0;
  int partitionFactor = 1;

  struct GroupCode : public Parma_GroupCode {
    apf::Mesh2* mesh;
    void run(int) {
      mesh->writeNative(outFile);
    }
  };

  void getConfig(int argc, char** argv)
  {
    if ( argc < 4 ) {
      if ( !PCU_Comm_Self() )
        printf("Usage: mpirun -np <inPartCount> %s <mesh> <outMesh> <factor>\n"
               "Reduce the part count of mesh from inPartCount to inPartCount/factor.\n",
               argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (argc == 5)
    {
      meshFile = argv[2];
      outFile = argv[3];
      partitionFactor = atoi(argv[4]);
    }
    if (argc == 4)
    {
      meshFile = argv[1];
      outFile = argv[2];
      partitionFactor = atoi(argv[3]);
    }
    PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  double t0 = PCU_Time();
  PCU_Protect();
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_null();
  getConfig(argc,argv);
  GroupCode code;
  gmi_model* nullModel = gmi_load(".null");
  code.mesh = apf::loadMdsMesh(nullModel, meshFile);
  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, partitionFactor, code);
  if (!PCU_Comm_Self())
    std::cout<<"[PUMI INFO] mesh collapsed to "<<PCU_Comm_Peers()/partitionFactor<<" parts in "
             <<PCU_Time()-t0<<"(sec)\n";

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
