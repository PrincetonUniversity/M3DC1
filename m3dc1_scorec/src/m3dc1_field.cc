#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include "PCU.h"
#include "apf.h"
#include "apfField.h"
#include "apfMesh2.h"

//*******************************************************
void synchronize_field(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;       

  int n = countComponents(f);
  double* sender_data = new double[n];
  double* dof_data = new double[n]; 

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(m,e) || (!m->isShared(e)&&!m->isGhosted(e)))
      continue;
      
    getComponents(f, e, 0, dof_data);

    if (m->isShared(e))
    {
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(apf::Copies,remotes,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
      }
    }
    if (m->isGhosted(e))
    {
      apf::Copies ghosts;
      m->getGhosts(e,ghosts);
      APF_ITERATE(apf::Copies,ghosts,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
      } 
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      setComponents(f, r, 0, sender_data);
    }
  delete [] dof_data;
  delete [] sender_data;
}
