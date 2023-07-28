/******************************************************************************

  (c) 2005-2020 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_mesh.h"
#include "m3dc1_model.h"
#include "gmi.h"
#include "apf.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include "PCU.h"
#include <vector>
#include <iostream>
#include <assert.h>
#include <stdlib.h>

using namespace apf;

int get_ent_gclasdim(Mesh2* m, MeshEntity* e)
{
  gmi_ent* clas=(gmi_ent*)m->toModel(e);
  return gmi_dim(m3dc1_model::instance()->model, clas);
}

// **********************************************
bool is_ent_original(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  if (mesh->isGhost(e)) return false;

  return (PCU_Comm_Self() == get_ent_ownpartid(mesh, e));
}

// **********************************************
int get_ent_ownpartid(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  int own_partid;

  if (mesh->hasTag(e, m3dc1_mesh::instance()->own_partid_tag))
    mesh->getIntTag(e, m3dc1_mesh::instance()->own_partid_tag, &own_partid);
  else
    own_partid=mesh->getOwner(e);
  return own_partid;
}

// **********************************************
MeshEntity* get_ent_owncopy(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  if (!(mesh->isShared(e))) // internal ent
    return e;

  int own_partid = get_ent_ownpartid(mesh, e);
  if (own_partid==PCU_Comm_Self()) return e;

  if (mesh->isShared(e)) 
  {
    Copies remotes;
    mesh->getRemotes(e,remotes);
    return remotes[own_partid];
  }
  if (mesh->isGhost(e))
  {
    Copies ghosts;
    mesh->getGhosts(e,ghosts);
    return ghosts[own_partid];
  }
  assert(0); // this should not called
  return NULL;
}

// *********************************************************
// return contiguous local entity ID if available
int get_ent_localid (Mesh2* mesh, MeshEntity* e)
// *********************************************************
{
  if (mesh->hasTag(e,m3dc1_mesh::instance()->local_entid_tag))
  {
    int local_id;
    mesh->getIntTag(e, m3dc1_mesh::instance()->local_entid_tag, &local_id);
    return local_id;
    // this is not true after adaptation
    // assert(local_id== getMdsIndex(mesh, e));
  }
  return getMdsIndex(mesh, e);
}

// *********************************************************
int get_ent_globalid(Mesh2* m, MeshEntity* e)
{
// *********************************************************
  MeshTag* tag = m->findTag("global_id");
  assert(m->hasTag(e, tag));
  int id;
  m->getIntTag(e, tag, &id);
  return id;
}

// exchange adjacency info of part boundary entities
void send_upadj(Mesh2* m, MeshEntity* e, int up_dim, MeshTag* tag)
{
  void* msg_send;
  MeshEntity** s_ent;
  int gid;

  // upward adjacency
  Adjacent adjacent;
  m->getAdjacent(e,up_dim,adjacent);
  int num_adj = adjacent.getSize();
  if (num_adj!=1)
    std::cout<<"[M3DC1 WARNING] (p"<<PCU_Comm_Self()<<") #upward adj of element "
             <<getMdsIndex(m,e)<<" is "<<num_adj<<"\n";

  size_t msg_size=sizeof(MeshEntity*)+num_adj*sizeof(int);
  Copies remotes;
  m->getRemotes(e,remotes);
  APF_ITERATE(Copies,remotes,rit)
  {
    int to = rit->first;
    msg_send = malloc(msg_size);
    s_ent = (MeshEntity**)msg_send; 
    *s_ent = rit->second; 
    int *s_data = (int*)((char*)msg_send+sizeof(MeshEntity*));
    for (int i=0; i<num_adj; ++i)
    {
      m->getIntTag(adjacent[i], tag, &gid);
      s_data[i] = gid;
    }
    PCU_Comm_Write(to, (void*)msg_send, msg_size);
    free(msg_send);
  }
}

void receive_adj(Mesh2* m, MeshTag* adj_pid_tag, MeshTag* adj_gid_tag)
{
  void *msg_recv;
  int pid_from;
  size_t msg_size;
  MeshEntity* e;
  int res=0;

  while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
  {
    e = *((MeshEntity**)msg_recv); 
    int* r_values = (int*)((char*)msg_recv+sizeof(MeshEntity*)); 
    int n = (msg_size-sizeof(MeshEntity*))/sizeof(int);

    m->setIntTag(e, adj_pid_tag, &pid_from);
    m->setIntTag(e, adj_gid_tag, &(r_values[0]));
  } // while
}

// this works only for elements (ent_dim==mesh_dim)
// *********************************************************
int get_ent_global2ndadj (Mesh2* m, int ent_dim, int adj_dim,
                        std::vector<MeshEntity*>& ents,
                        std::vector<int>& num_adj_ent,
                        std::vector<int>& adj_ent_pid,
                        std::vector<int>& adj_ent_gid)
// *********************************************************
{
  int mesh_dim=m->count(3)?3:2;
  assert(ent_dim==mesh_dim && adj_dim==ent_dim);

  //for all mesh entities of ent_dim-1 on part boundary, send adjcency info to remote copies
  PCU_Comm_Begin();
  MeshEntity* e;
  MeshTag* gid_tag = m->findTag("global_id");
  MeshTag* adj_pid_tag = m->createIntTag("remote adj pid", 1);
  MeshTag* adj_gid_tag = m->createIntTag("remote adj gid", 1);

  int down_size;
  MeshEntity* down_e;
  MeshIterator* it = m->begin(ent_dim-1);
  while ((e = m->iterate(it)))
  {
    if (m->isShared(e))
      send_upadj(m, e, ent_dim, gid_tag);
  }
  m->end(it);

  PCU_Comm_Send();
  
  receive_adj(m, adj_pid_tag, adj_gid_tag); 

  int total_num_adj=0, num_adj, pid, gid, myrank=PCU_Comm_Self();
  for (std::vector<MeshEntity*>::iterator vit=ents.begin(); vit!=ents.end(); ++vit)
  {
    e = *vit;
    Adjacent adjacent;
    getBridgeAdjacent(m, e, ent_dim-1, ent_dim, adjacent);
    num_adj = adjacent.getSize();
    for (int i=0; i<num_adj; ++i)
    {
      adj_ent_pid.push_back(myrank);
      m->getIntTag(adjacent[i], gid_tag, &gid);
      adj_ent_gid.push_back(gid);
    }

    Downward downward;
    down_size = m->getDownward(e, ent_dim-1, downward);

    for (int i=0; i<down_size; ++i)
    {
      down_e = downward[i];
      if (m->hasTag(down_e, adj_pid_tag))
      {
        assert (m->isShared(down_e));
        ++num_adj;
        m->getIntTag(down_e, adj_pid_tag, &pid);
        m->getIntTag(down_e, adj_gid_tag, &gid);
        adj_ent_pid.push_back(pid);
        adj_ent_gid.push_back(gid);
      }
    }
    num_adj_ent.push_back(num_adj);

// FIXME: wrong geom-class in 3D results in failure with the verification
#ifdef DEBUG
    int correct_num_adj=(mesh_dim==3)?5:3;

    int bdry_cnt=0;
    for (int i=0; i<down_size; ++i)
    {
      down_e = downward[i];
      if (get_ent_gclasdim(m,down_e)==ent_dim-1)
        ++bdry_cnt;
    }
    if (bdry_cnt+num_adj!=correct_num_adj)
    {
      std::cout<<"("<<PCU_Comm_Self()<<") elm "<<getMdsIndex(m,e)<<" num adj elm "<<num_adj<<"\n";
      for (int i=0; i<down_size; ++i)
      {
        down_e = downward[i];
        std::cout<<"("<<PCU_Comm_Self()<<") elm "<<getMdsIndex(m,e)<<" down "<<i
                 <<" g-dim "<<get_ent_gclasdim(m,down_e)<<"\n";
      }
    }
   // assert(bdry_cnt+num_adj==correct_num_adj);
#endif

    total_num_adj+=num_adj;
  }

  m->destroyTag(adj_pid_tag);
  m->destroyTag(adj_gid_tag);
  return total_num_adj;
}

// *********************************************************
void get_ent_numglobaladj (Mesh2* m, int ent_dim, int adj_dim,
                        std::vector<MeshEntity*>& ent_vec,
                        std::vector<int>& num_adj_ent)
// *********************************************************
{
  std::vector<int> adj_ent_pid;
  std::vector<int> adj_ent_gid;  
  get_ent_global2ndadj (m, ent_dim, adj_dim, ent_vec,
                        num_adj_ent, adj_ent_pid, adj_ent_gid);
}
