#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <apfZoltan.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

namespace {

const char* meshFile = 0;

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

void getConfig(int argc, char** argv)
{
  if ( argc < 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
}

void print_info (apf::Mesh2* m)
{
  if (!PCU_Comm_Self()) std::cout<<"\n===== mesh size =====\n";

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  int* own_entity_count = new int[4*PCU_Comm_Peers()];

  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=own_entity_count[i]=0;

  apf::MeshEntity* e;
  int self = PCU_Comm_Self();

  for (int d=0; d<4;++d)
  {
    local_entity_count[4*self+d] = m->count(d);
    apf::MeshIterator* it = m->begin(d); 
    while ((e = m->iterate(it)))
    {
      if (m->getOwner(e)==self)
        ++own_entity_count[4*PCU_Comm_Self()+d];
    }
    m->end(it);
  }

  int* global_local_entity_count = new int[4*PCU_Comm_Peers()];
  int* global_own_entity_count = new int[4*PCU_Comm_Peers()];

  MPI_Allreduce(local_entity_count, global_local_entity_count, 4*PCU_Comm_Peers(),
                MPI_INT, MPI_SUM, PCU_Get_Comm());

  MPI_Allreduce(own_entity_count, global_own_entity_count, 4*PCU_Comm_Peers(),
                MPI_INT, MPI_SUM, PCU_Get_Comm());

  if (!PCU_Comm_Self())
  {
    int* global_entity_count = new int[4];
    global_entity_count[0]=global_entity_count[1]=global_entity_count[2]=global_entity_count[3]=0;
    for (int d=0; d<4;++d)
    {
      for (int p=0; p<PCU_Comm_Peers();++p)
        global_entity_count[d] += global_own_entity_count[p*4+d];
    }

    for (int p=0; p<PCU_Comm_Peers(); ++p)
      std::cout<<"(p"<<p<<") # local ent: v "<<global_local_entity_count[p*4]
        <<", e "<<global_local_entity_count[p*4+1]
        <<", f "<<global_local_entity_count[p*4+2]
        <<", r "<<global_local_entity_count[p*4+3]<<"\n";
    std::cout<<"\n";
    for (int p=0; p<PCU_Comm_Peers(); ++p)
      if (global_own_entity_count[p*4])
        std::cout<<"(p"<<p<<") # own ent: v "<<global_own_entity_count[p*4]
          <<", e "<<global_own_entity_count[p*4+1]
          <<", f "<<global_own_entity_count[p*4+2]
          <<", r "<<global_own_entity_count[p*4+3]<<"\n";
    std::cout<<"\n";

    std::cout<<"# global ent: v "<<global_entity_count[0]<<", e "<<global_entity_count[1]
              <<", f "<<global_entity_count[2]<<", r "<<global_entity_count[3]<<"\n\n";

    delete [] global_entity_count;
  }

  delete [] local_entity_count;
  delete [] global_local_entity_count;
  delete [] own_entity_count;
  delete [] global_own_entity_count;
}
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  double t0 = PCU_Time();

  gmi_register_null();
  getConfig(argc,argv);
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = 0;
  m = apf::loadMdsMesh(g, argv[1]);
  apf::verify(m, false);
  //Parma_PrintPtnStats(m, "");
  print_info(m);

  if (!PCU_Comm_Self())
    std::cout<<"check_smb: completed\n";
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
