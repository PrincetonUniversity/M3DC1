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

using std::cout;

namespace {

  const char* modelFile = 0;
  const char* meshFile = 0;
  const char* outFile = 0;
  int partitionFactor = 1;
  void write_topo(apf::Mesh2*, const char*, int);

  struct GroupCode : public Parma_GroupCode {
    apf::Mesh2* mesh;
    void run(int) {
      mesh->writeNative(outFile);
      write_topo(mesh, "topo", 0);
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

#include <assert.h>

void write_topo(apf::Mesh2* m, const char* filename, int start_index)
{
  apf::MeshEntity* e;
  char filename_buff[256];
  apf::MeshIterator* it;

  apf::MeshTag* tag = m->findTag("_global_id_");
  int id, dim=2;
  if (m->count(3)) dim=3;

  sprintf(filename_buff, "%s",filename);
  FILE * fp =fopen(filename_buff, "w");

  fprintf(fp, "%d\t%d\t%d\n", dim, m->count(0), m->count(dim));

  it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    apf::Vector3 xyz;
    assert(m->hasTag(e, tag));
    m->getIntTag(e, tag, &id);
    m->getPoint(e, 0, xyz);
    fprintf(fp, "%d\t%d\t%d\t%lf\t%lf\t%lf\n", id,
            m->getModelType(m->toModel(e)), m->getModelTag(m->toModel(e)),
            xyz[0],xyz[1],xyz[2]);
  }
  m->end(it);

  it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    assert(m->hasTag(e, tag));
    m->getIntTag(e, tag, &id);
    apf::Downward down;
    int num_down = m->getDownward(e,0,down);
    assert (num_down==3 || num_down==6);
    switch (num_down)
    {
      case 3: 
	    fprintf(fp, "%d\t2\t%d\t%d\t%d\t%d\n", id,
	    m->getModelType(m->toModel(e)), // m->getModelTag(m->toModel(e)),
	    getMdsIndex(m, down[0]), getMdsIndex(m, down[1]), getMdsIndex(m, down[2]));
	    break;
      case 6: 
	    fprintf(fp, "%d\t3\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", id,
            m->getModelType(m->toModel(e)),// m->getModelTag(m->toModel(e)),
            getMdsIndex(m, down[0]), getMdsIndex(m, down[1]), getMdsIndex(m, down[2]),
	    getMdsIndex(m, down[3]), getMdsIndex(m, down[4]), getMdsIndex(m, down[5]));
	    break;
      default: break;
    }
  }
  m->end(it);
  fclose(fp);
}

}

void create_globalid(apf::Mesh2* m, int dim)
{ 
  apf::MeshTag* tag = m->findTag("_global_id_");
  if (!tag)
    tag = m->createIntTag("_global_id_",1);

  apf::Sharing* sh = apf::getSharing(m);
  apf::MeshEntity* e;
  apf::MeshEntity* remote_ent;
  int num_own=0;
  
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    if (sh->isOwned(e))
      ++num_own;
  m->end(it);
  
  int num_global=0;
  MPI_Allreduce(&num_own, &num_global, 1, MPI_INT, MPI_SUM, PCU_Get_Comm());
  if (!PCU_Comm_Self())
    std::cout<<__func__<<": #ent of dim "<<dim<<" "<<num_global<<"\n";

  PCU_Exscan_Ints(&num_own,1);
  int initial_id=num_own;
  
  PCU_Comm_Begin();
  it = m->begin(dim);
  while ((e = m->iterate(it)))
  { 
    if (!sh->isOwned(e))
      continue;
    
    m->setIntTag(e, tag, &initial_id);
    apf::CopyArray remotes;
    sh->getCopies(e, remotes);
    
    APF_ITERATE(apf::CopyArray, remotes, rit)
    { 
      PCU_COMM_PACK(rit->peer, rit->entity);
      PCU_Comm_Pack(rit->peer, &initial_id, sizeof(int));
    }
    ++initial_id;
  }
  m->end(it);

  PCU_Comm_Send();
  int global_id;
  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
    { 
      PCU_COMM_UNPACK(remote_ent);
      PCU_Comm_Unpack(&global_id, sizeof(int));
      m->setIntTag(remote_ent, tag, &global_id);
    }
}

void delete_globalid(apf::Mesh2* m)
{
  apf::MeshTag* tag = m->findTag("_global_id_");
  if (!tag) return;

  for (int i=0; i<4; ++i)
    apf::removeTagFromDimension(m, tag, i);
  m->destroyTag(tag);
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

  create_globalid(code.mesh, 0);
  code.mesh->count(3)>0?create_globalid(code.mesh,3):create_globalid(code.mesh,2);

  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, partitionFactor, code);
  if (!PCU_Comm_Self())
    std::cout<<"[PUMI INFO] mesh collapsed to "<<PCU_Comm_Peers()/partitionFactor<<" parts in "
             <<PCU_Time()-t0<<"(sec)\n";

 delete_globalid(code.mesh);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
