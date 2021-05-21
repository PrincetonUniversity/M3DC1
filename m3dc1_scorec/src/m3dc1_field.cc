#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include "PCU.h"
#include "apf.h"
#include "apfField.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include <assert.h>
#include <vector>

void group_complex_dof (apf::Field* field, int option)
{
  //if (!PCU_Comm_Self()) cout<<" regroup complex number field with option "<<option<<endl;
  int num_dof_double = countComponents(field);
  assert(num_dof_double/6%2==0);
  int num_dof = num_dof_double/2;
  std::vector<double> dofs(num_dof_double);
  std::vector<double> newdofs(num_dof_double);

  apf::MeshIterator* ent_it=m3dc1_mesh::instance()->mesh->begin(0);
  apf::MeshEntity* e;
  while ((e = m3dc1_mesh::instance()->mesh->iterate(ent_it)))
  {
    getComponents(field, e, 0, &(dofs[0]));
    for (int j=0; j<num_dof/6; j++)
    {
      if (option)
      {
        for (int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+k)=dofs.at(2*j*6+2*k);
          newdofs.at(2*j*6+6+k)=dofs.at(2*j*6+2*k+1);
        }
      }
      else
      {
        for (int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+2*k)=dofs.at(2*j*6+k);
          newdofs.at(2*j*6+2*k+1)=dofs.at(2*j*6+6+k);
        }
      }
    }
    setComponents(field, e, 0, &(newdofs[0]));
  }
  m3dc1_mesh::instance()->mesh->end(ent_it);
}

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
