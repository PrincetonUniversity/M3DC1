#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_null.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <iostream>

#ifdef PPPL
#include <SimLicense.h>
#ifdef STELLAR
char simLic[128]="/home/PPPL/simmetrix/license/simmetrix.lic";
#else
char simLic[128]="/usr/pppl/Simmetrix/simmodsuite.lic";
#endif
#else
char simLic[128]="/net/common/meshSim/license/license.txt";
#endif

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (argc < 3) {
    if (0==PCU_Comm_Self())
      std::cout<< "usage: " << argv[0] << " <simmetrix mesh (.sms)> <pumi mesh (.smb)>\n";
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }

#ifdef PPPL
  SimLicense_start("geomsim_core,geomsim_adv,meshsim_surface,meshsim_adapt,meshsim_adv",simLic);
#else
  Sim_readLicenseFile(simLic);
#endif
  std::cout<<argv[0]<<": license "<<simLic<<"\n";

  SimPartitionedMesh_start(&argc,&argv);
  gmi_sim_start();

  gmi_register_sim();

  pParMesh sim_mesh = PM_load(argv[1], NULL, NULL);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);

  gmi_register_null();
  gmi_model* mdl = gmi_load(".null");
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
//  apf::printStats(mesh);

  std::cout<<argv[0]<<": Simmetrix mesh converted to "<<argv[2]<<"\n";
  std::cout<<"=== mesh entity count ===\n";
  for (int i=0; i<PCU_Comm_Peers(); ++i)
  {
    if (PCU_Comm_Self()==i)
    {
      std::cout<<"P"<<PCU_Comm_Self()<<": V "<<mesh->count(0)<<", E "
               <<mesh->count(1)<<", F "
               <<mesh->count(2);
      if (mesh->count(3))  std::cout<<", R "<<mesh->count(3);
      std::cout<<"\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  std::cout<<"\n";

  apf::destroyMesh(simApfMesh);

  M_release(sim_mesh);
  mesh->writeNative(argv[2]);
  mesh->destroyNative();
  apf::destroyMesh(mesh);

  gmi_sim_stop();
  SimPartitionedMesh_stop();

#ifdef PPPL
  SimLicense_stop();
#endif

  PCU_Comm_Free();
  MPI_Finalize();
}
