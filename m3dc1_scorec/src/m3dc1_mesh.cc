/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_mesh.h"
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "PCU.h"
#include "gmi_base.h"
#include <map>
#include <vector>
#include <set>
#include <string.h>
#include <iostream>
#include <assert.h>

using namespace apf;

//*******************************************************
void compute_globalid(apf::Mesh2* m, int d)
//*******************************************************
{
  apf::MeshTag* tag = m->findTag("global_id");
  if (!tag)  // update existing tag
    tag = m->createIntTag("global_id",1);

  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[d];

  apf::MeshEntity* e;
  PCU_Exscan_Ints(&num_own_ent,1);
  int start=num_own_ent;
  
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(d);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(m,e)) continue;
    m->setIntTag(e, tag, &start);
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&start,sizeof(int));
    }
    if (m->isGhosted(e))
    {
      apf::Copies ghosts;
      m->getGhosts(e, ghosts);
      APF_ITERATE(apf::Copies, ghosts, it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&start,sizeof(int));
      }
    }
    ++start;
  }
  m->end(it);
  PCU_Comm_Send();

  int value;
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&value,sizeof(int));
      m->setIntTag(r, tag, &value);
    }
  }
}

void set_remote(Mesh2* m, MeshEntity* e, int p, MeshEntity* r)
{
  Copies remotes;
  m->getRemotes(e,remotes);
  bool found=false;
  APF_ITERATE(Copies, remotes, it)
    if (it->first==p) 
    {
      found=true;
      break;
    }
  if (found)
  {
    remotes[p] = r;
    m->setRemotes(e,remotes);
  }
  else
    m->addRemote(e, p, r);
}

// m3dc1_mesh
// *********************************************************
m3dc1_mesh::m3dc1_mesh()
// *********************************************************
{
  mesh = NULL;
  field_container=NULL;
  reset();
  local_entid_tag=own_partid_tag=num_global_adj_node_tag=num_own_adj_node_tag=NULL;
}

// *********************************************************
m3dc1_mesh::~m3dc1_mesh()
// *********************************************************
{}

void m3dc1_mesh::clean()
{
 // destroy field AND numbering
  if (field_container)
  {
    for (std::map<FieldID, m3dc1_field*>::iterator f_it=field_container->begin(); f_it!=field_container->end();)
    {
      if (!PCU_Comm_Self()) std::cout<<" destroy field "<<getName(f_it->second->get_field())<<std::endl;
      FieldID id = f_it->first;
      std::map<FieldID, m3dc1_field*>::iterator it_next=++f_it;
      m3dc1_field_delete(&id);
      f_it=it_next;
    }
    //field_container->clear();
  }
  delete field_container; field_container=0;

  // destroy tag data
  for (int d=0; d<4; ++d)
  {
    removeTagFromDimension(mesh, local_entid_tag, d);
    removeTagFromDimension(mesh, own_partid_tag, d);
  }
  removeTagFromDimension(mesh, num_global_adj_node_tag, 0);
  removeTagFromDimension(mesh, num_own_adj_node_tag, 0);

  // destroy tag
  mesh->destroyTag(local_entid_tag);
  mesh->destroyTag(own_partid_tag);
  mesh->destroyTag(num_global_adj_node_tag);
  mesh->destroyTag(num_own_adj_node_tag);
}

m3dc1_mesh* m3dc1_mesh::_instance=NULL;
m3dc1_mesh* m3dc1_mesh::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_mesh();
  return _instance;
}

// *********************************************************
void m3dc1_mesh::reset()
// *********************************************************
{
  for (int i=0; i<4; ++i)
  {
    num_local_ent[i] = 0;
    num_own_ent[i] = 0;
    num_global_ent[i] = 0;
  }
}

// **********************************************
void push_new_entities (Mesh2* mesh, std::map<MeshEntity*, MeshEntity*>& new_entities)
// **********************************************
{
  MeshEntity* e;
  MeshEntity* new_e;
  int value;

  PCU_Comm_Begin();

  for (std::map<MeshEntity*, MeshEntity*>::const_iterator ent_it=new_entities.begin(); ent_it!=new_entities.end();++ent_it)
  {
    e = ent_it->first; 
    new_e = ent_it->second;             // current global pid

    Copies remotes;
    mesh->getRemotes(e, remotes);

    APF_ITERATE(Copies, remotes, it)
    {
       if (it->first==m3dc1_model::instance()->prev_plane_partid) continue;
       PCU_COMM_PACK(it->first, it->second);
       PCU_COMM_PACK(it->first, new_e);
       PCU_Comm_Pack(it->first, &(m3dc1_model::instance()->prev_plane_partid), sizeof(int));
    }
  }
  PCU_Comm_Send();

   // receive phase begins
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      PCU_COMM_UNPACK(new_e);
      PCU_Comm_Unpack(&value,sizeof(int));
      mesh->addRemote(e, value, new_e);
    }
  }
}

  
// **********************************************
void bounce_orig_entities (Mesh2* mesh, std::vector<MeshEntity*>& mesh_ents, int rank_to)
// **********************************************
{
  MeshEntity* e;
  PCU_Comm_Begin();
  int num_remote, index, own_partid, local_rank = PCU_Comm_Self();
  for (std::vector<MeshEntity*>::iterator ent_it=mesh_ents.begin(); ent_it!=mesh_ents.end();++ent_it)
  {
    e = *ent_it;
    own_partid=get_ent_ownpartid(mesh,e);

    Copies remotes;
    mesh->getRemotes(e, remotes);
    num_remote = remotes.size();
    PCU_COMM_PACK(rank_to, remotes[rank_to]);
    PCU_Comm_Pack(rank_to, &own_partid,sizeof(int));  
    PCU_Comm_Pack(rank_to, &num_remote,sizeof(int));

    MeshEntity** remote_copy = new MeshEntity*[num_remote];
    int* remote_pid = new int[num_remote];
    index=0;
    APF_ITERATE(Copies, remotes, it)
    { 
      if (it->first==rank_to)
      {
        remote_pid[index] = local_rank;
        remote_copy[index] = e;
      }
      else
      {
        remote_pid[index] = it->first;
        remote_copy[index] = it->second;
      }
      ++index;
    }  // APF_ITERATE
    PCU_Comm_Pack(rank_to, &(remote_pid[0]), num_remote*sizeof(int));
    PCU_Comm_Pack(rank_to, &(remote_copy[0]), num_remote*sizeof(MeshEntity*));
    delete [] remote_copy;
    delete [] remote_pid;
  } 
  PCU_Comm_Send();

  // receive phase begins
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      MeshEntity* e;
      int own_partid, num_remote;
      PCU_COMM_UNPACK(e);
      PCU_Comm_Unpack(&own_partid,sizeof(int));
      mesh->setIntTag(e, m3dc1_mesh::instance()->own_partid_tag, &own_partid);
      PCU_Comm_Unpack(&num_remote,sizeof(int));
      int* remote_pid = new int [num_remote];
      MeshEntity** remote_copy = new MeshEntity*[num_remote];
      PCU_Comm_Unpack(&remote_pid[0],num_remote*sizeof(int));
      PCU_Comm_Unpack(&remote_copy[0],num_remote*sizeof(MeshEntity*));
      for (int i=0; i<num_remote; ++i)
        mesh->addRemote(e, remote_pid[i], remote_copy[i]);
      delete [] remote_pid;
      delete [] remote_copy;
    }
  }
}

// *********************************************************
void m3dc1_receiveVertices(Mesh2* mesh, std::vector<MeshEntity*>* vertices,
     MeshTag* partbdry_id_tag, 
     std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  int proc_grp_rank = PCU_Comm_Self()/m3dc1_model::instance()->group_size;
  int proc_grp_size = m3dc1_model::instance()->group_size;
  int local_rank = PCU_Comm_Self();
  int num_ent, num_remote, own_partid;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      double* v_coords = new double[num_ent*3];
      double* v_params = new double[num_ent*3];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];

      PCU_Comm_Unpack(&(v_coords[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(v_params[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid;
      if (num_remote) 
      {
        e_rmt_pid = new int[num_remote];
        PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));
      }
      // create vertex
      int rinfo_pos=0;
      Vector3 coord;
      Vector3 param;     
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
        for (int k=0; k<3; ++k)
        {
          coord[k] = v_coords[index*3+k];
          param[k] = v_params[index*3+k];
        }
        new_ent = mesh->createVertex((ModelEntity*)geom_ent, coord, param);
        vertices->push_back(new_ent);

        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        own_partid=proc_grp_rank*proc_grp_size+e_own_partid[index];
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &own_partid);

        if (own_partid==local_rank)
          ++(m3dc1_mesh::instance()->num_own_ent[0]);

        if (e_global_id[index]!=-1)
        {
          partbdry_entities[0][e_global_id[index]] = new_ent;
          mesh->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }

        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          mesh->addRemote(new_ent, proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index
      delete [] v_coords;
      delete [] v_params;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      if (num_remote)
        delete [] e_rmt_pid;
      m3dc1_mesh::instance()->num_local_ent[0] = mesh->count(0);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void m3dc1_receiveEdges(Mesh2* mesh, std::vector<MeshEntity*>* vertices,
        std::vector<MeshEntity*>* edges, 
        MeshTag* partbdry_id_tag, std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  int local_rank=PCU_Comm_Self();
  int proc_grp_rank = local_rank/m3dc1_model::instance()->group_size;
  int proc_grp_size = m3dc1_model::instance()->group_size;

  int num_ent, num_remote, own_partid;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  Downward down_ent; 
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      int* e_down_lid = new int[num_ent*2];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];
      PCU_Comm_Unpack(&(e_down_lid[0]), num_ent*2*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid = new int[num_remote];
      if (num_remote) PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));

      // create edge
      int rinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
        down_ent[0] =  vertices->at(e_down_lid[index*2]);
        down_ent[1] =  vertices->at(e_down_lid[index*2+1]);

        new_ent = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)geom_ent, down_ent);
        edges->push_back(new_ent);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        own_partid=proc_grp_rank*proc_grp_size+e_own_partid[index];
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &own_partid);

        if (own_partid==local_rank)
          ++(m3dc1_mesh::instance()->num_own_ent[1]);

        if (e_global_id[index]!=-1)
        {
          partbdry_entities[1][e_global_id[index]] = new_ent;
          mesh->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }
        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          //ent_bps[new_ent].insert(proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos]);
          mesh->addRemote(new_ent, proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index      
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      delete [] e_rmt_pid;
      m3dc1_mesh::instance()->num_local_ent[1] = mesh->count(1);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())

}

// *********************************************************
void m3dc1_receiveFaces(Mesh2* mesh, 
     std::vector<MeshEntity*>* edges, std::vector<apf::MeshEntity*>* faces)
// *********************************************************
{
  int num_ent, num_down;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  Downward down_ent;
  int local_rank = PCU_Comm_Self();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      PCU_Comm_Unpack(&num_down, sizeof(int));

      int* f_down_num = new int[num_ent];
      int* e_down_lid = new int[num_down];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      PCU_Comm_Unpack(&(f_down_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_down_lid[0]), num_down*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      // create face
      int dinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent =gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
        num_down = f_down_num[index];
        for (int i=0; i<num_down; ++i)
        {
          down_ent[i] = edges->at(e_down_lid[dinfo_pos]);
          ++dinfo_pos;
        }
        if (num_down==3)
          new_ent= mesh->createEntity(apf::Mesh::TRIANGLE, (ModelEntity*)geom_ent, down_ent);
        else
          new_ent= mesh->createEntity(apf::Mesh::QUAD, (ModelEntity*)geom_ent, down_ent);
        faces->push_back(new_ent);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &local_rank);
      } // for index
      delete [] f_down_num;
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      m3dc1_mesh::instance()->num_local_ent[2] = mesh->count(2);
      m3dc1_mesh::instance()->num_own_ent[2] = mesh->count(2);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void m3dc1_stitchLink(Mesh2* mesh, MeshTag* partbdry_id_tag, 
                      std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  PCU_Comm_Begin();
  MeshEntity* e;
  
  int ent_info[3];
  int global_id; 
  for (int dim=0; dim<2; ++dim)
  {
    for (std::map<int, MeshEntity*>::iterator ent_it = partbdry_entities[dim].begin();
         ent_it!= partbdry_entities[dim].end(); ++ent_it)
    {    
      e = ent_it->second;
      Copies remotes;
      mesh->getRemotes(e, remotes);
      mesh->getIntTag(e, partbdry_id_tag, &global_id);
      ent_info[0] = dim;
      ent_info[1] = global_id; 
      ent_info[2] = PCU_Comm_Self(); 
      APF_ITERATE(Copies, remotes, it)
      {  
        PCU_Comm_Pack(it->first, &(ent_info[0]),3*sizeof(int));
        PCU_COMM_PACK(it->first, e);
      }
    }
  } //  for (int dim=0; dim<2; ++dim)

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      int sender_ent_info[3];
      MeshEntity* sender;
      PCU_Comm_Unpack(&(sender_ent_info[0]), 3*sizeof(int));
      PCU_COMM_UNPACK(sender);
      e = partbdry_entities[sender_ent_info[0]][sender_ent_info[1]];
      set_remote(mesh, e, sender_ent_info[2], sender);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void sendEntities(Mesh2* mesh, int dim, MeshTag* partbdry_id_tag)
// *********************************************************
{
  int own_partid, num_ent = m3dc1_mesh::instance()->num_local_ent[dim];
  if (dim>2 || !num_ent) return;

  double* v_coords;
  double* v_params;
  int* e_geom_type = new int[num_ent];
  int* e_geom_tag = new int[num_ent];
  int* f_down_num; // for faces
  int* e_own_partid; // for vertices and edges
  int* e_global_id; // for vertices and edges
  int* e_rmt_num;  // for vertices and edges
  
  MeshEntity* e;
  gmi_ent* geom_ent;

  if (dim==0)
  {
    v_coords = new double[num_ent*3];  // for vertices
    v_params = new double[num_ent*3];
  }
  if (dim==2)
    f_down_num = new int[num_ent]; // for faces
  else // for vertices and edges
  {
    e_own_partid = new int[num_ent]; 
    e_global_id = new int[num_ent]; 
    e_rmt_num = new int[num_ent];  
  }

  int num_down=0, num_remote=0;
  MeshIterator* it = mesh->begin(dim);
  while ((e = mesh->iterate(it)))
  {
    switch (mesh->getType(e))
    {
      case apf::Mesh::TRIANGLE: num_down+=3; break;
      case apf::Mesh::QUAD: num_down+=4; break;
      case apf::Mesh::EDGE: num_down+=2; break;
      default: break;
    }

    if (mesh->isShared(e))
    {
      Copies remotes;
      mesh->getRemotes(e, remotes);
      num_remote+=remotes.size();
    }
  }
  mesh->end(it);

  std::vector<int> e_down_lid;
  int* e_rmt_pid = new int[num_remote];

  int global_id, rinfo_pos=0, index=0, num_down_ent;

  it = mesh->begin(dim);
  Vector3 coord;
  Vector3 param;
  Downward down_ent;
  while ((e = mesh->iterate(it)))
  {
    switch (getDimension(mesh, e))
    {
      case 0: mesh->getPoint(e, 0, coord);
              mesh->getParam(e, param);
              for (int i=0; i<3; ++i)
              {
                v_coords[index*3+i] = coord[i];
                v_params[index*3+i] = param[i];
              }
              break;
       default: { 
                  num_down_ent =  mesh->getDownward(e, dim-1, down_ent); 
                  if (dim==2) 
                    f_down_num[index] = num_down_ent;
                  for (int i=0; i<num_down_ent; ++i)
                    e_down_lid.push_back(get_ent_localid(mesh, down_ent[i]));
                }
    } // switch

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    e_geom_type[index] = gmi_dim(m3dc1_model::instance()->model, geom_ent);
    e_geom_tag[index] = gmi_tag(m3dc1_model::instance()->model, geom_ent);

    if (dim!=2) // for vertices and edges
    {
      own_partid=get_ent_ownpartid(mesh,e);
      e_own_partid[index] = own_partid;
      Copies remotes;
      mesh->getRemotes(e, remotes);
      if (mesh->hasTag(e,partbdry_id_tag)) // part-bdry entity
      {
        mesh->getIntTag(e, partbdry_id_tag, &global_id);
        e_global_id[index] = global_id; 
        e_rmt_num[index] = remotes.size();
        APF_ITERATE(Copies, remotes, it)
        {
          e_rmt_pid[rinfo_pos]=it->first;
          ++rinfo_pos;
        }
      }
      else
      {
        e_global_id[index] = -1; 
        e_rmt_num[index] = 0;
      }
    }
    ++index;
  }
  mesh->end(it);

  int proc=PCU_Comm_Self()+m3dc1_model::instance()->group_size;
  while (proc<PCU_Comm_Peers())
  {  
    PCU_Comm_Pack(proc, &num_ent, sizeof(int));
    switch (dim)
    { 
      case 0: PCU_Comm_Pack(proc, &(v_coords[0]),num_ent*3*sizeof(double));
              PCU_Comm_Pack(proc, &(v_params[0]),num_ent*3*sizeof(double));
              break;
      default: if (dim==2)
               {
                 PCU_Comm_Pack(proc, &num_down, sizeof(int));
                 PCU_Comm_Pack(proc, &(f_down_num[0]),num_ent*sizeof(int));
               }
               PCU_Comm_Pack(proc, &e_down_lid.at(0), num_down*sizeof(int));
    }
    PCU_Comm_Pack(proc, &(e_geom_type[0]),num_ent*sizeof(int));
    PCU_Comm_Pack(proc, &(e_geom_tag[0]),num_ent*sizeof(int));
    if (dim<2)  
    {      
      PCU_Comm_Pack(proc, &(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_global_id[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_rmt_num[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &num_remote, sizeof(int));
      if (num_remote)
        PCU_Comm_Pack(proc, &(e_rmt_pid[0]), num_remote*sizeof(int));
    }
    proc+=m3dc1_model::instance()->group_size;
  }

  if (dim==0)
  {
    delete [] v_coords;
    delete [] v_params;
  }
  delete [] e_geom_type;
  delete [] e_geom_tag;

  if (dim==2)
    delete [] f_down_num;
  else 
  {    
    delete [] e_own_partid;
    delete [] e_global_id;
    delete [] e_rmt_num;
    if (num_remote) delete [] e_rmt_pid;
  }
}
// *********************************************************
void  assign_uniq_partbdry_id(Mesh2* mesh, int dim, MeshTag* partbdry_id_tag)
// *********************************************************
{
  MeshEntity* e;

  int own_partid, num_own_partbdry=0, local_rank=PCU_Comm_Self();

  MeshIterator* it = mesh->begin(dim);
  while ((e = mesh->iterate(it)))
  {
    own_partid = get_ent_ownpartid(mesh, e);
    if (own_partid==local_rank && mesh->isShared(e))
      ++num_own_partbdry;
  }
  mesh->end(it);

  PCU_Exscan_Ints(&num_own_partbdry,1);
  int initial_id=num_own_partbdry;

  PCU_Comm_Begin();
  it = mesh->begin(dim);
  while ((e = mesh->iterate(it)))
  {
    own_partid = get_ent_ownpartid(mesh, e);
    if (own_partid==local_rank && mesh->isShared(e))
    {
      mesh->setIntTag(e, partbdry_id_tag, &initial_id);
      Copies remotes;
      mesh->getRemotes(e, remotes);
      APF_ITERATE(Copies, remotes, it)
      {
        PCU_COMM_PACK(it->first, it->second);
        PCU_Comm_Pack(it->first, &initial_id, sizeof(int));
      }
      ++initial_id;
    }
  }
  mesh->end(it);

  PCU_Comm_Send();
  int global_id;
  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
    {
      MeshEntity* remote_ent;
      PCU_COMM_UNPACK(remote_ent);
      PCU_Comm_Unpack(&global_id, sizeof(int));
      mesh->setIntTag(remote_ent, partbdry_id_tag, &global_id);
    }
}

int get_local_planeid(int p)
{
  return p/m3dc1_model::instance()->group_size;
}

void update_remotes(apf::Mesh2* mesh, MeshEntity* e)
{
  int local_rank=PCU_Comm_Self();

  Copies remotes;
  Copies new_remotes;
  Parts parts;
  mesh->getRemotes(e,remotes);

  APF_ITERATE(Copies,remotes,rit)
  {
    int p=rit->first;
    if (get_local_planeid(p)==m3dc1_model::instance()->local_planeid) // the same local plane
    {
      new_remotes[p]=rit->second; 
      parts.insert(p);
    }
  }
  parts.insert(local_rank);
  mesh->clearRemotes(e);
  mesh->setRemotes(e,new_remotes);
  mesh->setResidence(e, parts);  // set pclassification
}
    
void m3dc1_mesh::remove_wedges()
{       
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshEntity* e;
  int local_rank = PCU_Comm_Self();
  int num_2d_face=mesh->count(3);
 
  // delete regions
  MeshIterator* ent_it = mesh->begin(3);
  while ((e = mesh->iterate(ent_it)))
    mesh->destroy(e);
  mesh->end(ent_it);
  m3dc1_mesh::instance()->num_own_ent[3]=m3dc1_mesh::instance()->num_local_ent[3]=0;
  // delete faces
  int cnt=0;
  ent_it = mesh->begin(2);
  while ((e = mesh->iterate(ent_it)))
  {
    // original entity for 2D plane
    if (cnt<num_2d_face)
      update_remotes(mesh, e);
    else
       mesh->destroy(e);
    ++cnt;
  }
  mesh->end(ent_it);
  m3dc1_mesh::instance()->num_local_ent[2] = mesh->count(2);

  // remote edges and vertices
  for (int dim=1; dim>=0; --dim)
  {
    ent_it = mesh->begin(dim);
    while ((e = mesh->iterate(ent_it)))
    {
      apf::Adjacent adjacent;
      mesh->getAdjacent(e,dim+1,adjacent);
      if (adjacent.getSize()) // original entity for 2D plane
        update_remotes(mesh, e);
      else
        mesh->destroy(e);
    }
    mesh->end(ent_it);
    m3dc1_mesh::instance()->num_local_ent[dim]=mesh->count(dim);
  }

  mesh->acceptChanges(); // update partition model
  changeMdsDimension(mesh, 2);
  initialize(false);

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Wedges between 2D planes removed ***\n";
}



void update_field (int field_id, int ndof_per_value, int num_2d_vtx, MeshEntity** remote_vertices);
void update_field_3d (int field_id, std::map<int, apf::Field*>* new_fields,
                      int num_2d_vtx, 
                      std::vector<MeshEntity*>* vertices, MeshEntity** remote_vertices);

// copy master plane to non-master
// new version which works with both initial & adapted mesh
void copy_master_plane(apf::Mesh2* mesh, MeshTag* local_entid_tag, 
     std::vector<MeshEntity*>* vertices,
     std::vector<MeshEntity*>* edges,
     std::vector<MeshEntity*>* faces,
     MeshEntity** remote_vertices, MeshEntity** remote_edges, MeshEntity** remote_faces, 
     std::vector<MeshEntity*>& btw_plane_edges, std::vector<MeshEntity*>& btw_plane_faces, 
     std::vector<MeshEntity*>& btw_plane_regions)
{
  int local_partid=PCU_Comm_Self();
  apf::MeshIterator* ent_it;
  MeshEntity* e;
  MeshEntity* new_ent;
  gmi_ent* geom_ent; 

  int index, new_id, local_id, local_id_0, local_id_1, local_id_2;
  double local_plane_phi = m3dc1_model::instance()->get_phi(local_partid);
  double next_plane_phi = m3dc1_model::instance()->get_phi(m3dc1_model::instance()->next_plane_partid);
  bool flip_wedge=false;
  if (local_plane_phi>next_plane_phi) // last plane
    flip_wedge=true;

  // create remote copy of vertices on next plane
  gmi_ent* new_geom_ent = NULL;
  int num_local_vtx=mesh->count(0);
  Vector3 cur_coord;
  Vector3 new_coord;

  index=0;
  std::vector<MeshEntity*>::iterator eit;
  for (eit=vertices->begin(); eit!=vertices->end(); ++eit)
  {
    e = *eit;
    mesh->getPoint(e, 0, cur_coord);
    new_coord[0]=cur_coord[0];
    new_coord[1]=cur_coord[1];
    new_coord[2]=next_plane_phi;
    Vector3 param(0,0,0);     
    mesh->getParam(e,param);
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    assert(geom_ent);
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_ent = mesh->createVertex((ModelEntity*)new_geom_ent, new_coord, param);
    new_id = num_local_vtx+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_vertices[index++]=new_ent;
    // update the z-coord of the current coordinate based on the prev plane partid
    cur_coord[2] = local_plane_phi;
    mesh->setPoint(e, 0, cur_coord);
  }

  // create remote copy of edges on next plane
  int num_local_edge=mesh->count(1);
  Downward down_vtx;
  Downward new_down_vtx;
 
  index=0;
  for (eit=edges->begin(); eit!=edges->end(); ++eit)
  {  
    e = *eit;
    mesh->getDownward(e, 0, down_vtx);
    local_id_0 = get_ent_localid(mesh, down_vtx[0]);
    local_id_1 = get_ent_localid(mesh, down_vtx[1]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_down_vtx[0] = remote_vertices[local_id_0];
    new_down_vtx[1] = remote_vertices[local_id_1];
    new_ent = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);

    new_id = num_local_edge+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_edges[index++]=new_ent;
  }

  // create remote copy of faces on next plane
  int num_local_face=mesh->count(2);
  Downward down_edge;
  Downward new_down_edge;

  index=0;
  for (eit=faces->begin(); eit!=faces->end(); ++eit) 
  {  
    e = *eit;
    mesh->getDownward(e, 1, down_edge);
    local_id_0 = get_ent_localid(mesh, down_edge[0]);
    local_id_1 = get_ent_localid(mesh, down_edge[1]);
    local_id_2 = get_ent_localid(mesh, down_edge[2]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_down_edge[0] = remote_edges[local_id_0];
    new_down_edge[1] = remote_edges[local_id_1];
    new_down_edge[2] = remote_edges[local_id_2];

    new_ent = mesh->createEntity(apf::Mesh::TRIANGLE, (ModelEntity*)new_geom_ent, new_down_edge);
    new_id = num_local_face+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_faces[index++]=new_ent;
  }

  // create edges, faces and regions between planes. 
  // if edges are classified on geom edge, a new face is classified on geom_surface_face
  // otherwise, a new face is classified on geom_region

  Downward quad_faces;
  Downward wedge_faces;

  int edge_counter=0, face_counter=0, rgn_counter=0;
  std::vector<MeshEntity*> vertex_edges;
  std::vector<MeshEntity*> edge_faces;

  Downward edgesNextPlane;
  Downward edgesBtwPlane;
  int num_upward;

  for (eit=faces->begin(); eit!=faces->end(); ++eit)
  {
    e = *eit;
    mesh->getDownward(e, 0, down_vtx);
    mesh->getDownward(e, 1, down_edge);
    for (int pos=0; pos<3; ++pos)
    {
      local_id = get_ent_localid(mesh, down_edge[pos]);
      edgesNextPlane[pos] = remote_edges[local_id];
    }

    /**create edges between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the new created face btw planes*/
      edgesBtwPlane[pos]=NULL;
      vertex_edges.clear();
      num_upward = mesh->countUpward(down_vtx[pos]);
      for (int i=0; i<num_upward; ++i)
        vertex_edges.push_back(mesh->getUpward(down_vtx[pos], i));

      for (unsigned int i=0; i<vertex_edges.size();++i)
      {
        // get the local id
        local_id = get_ent_localid(mesh, vertex_edges[i]);
        if (local_id>=num_local_edge) // edge is between planes
        {
          edgesBtwPlane[pos]=vertex_edges[i];
          break; // get out of for loop
        }
      }

      if (edgesBtwPlane[pos]!=NULL) continue;

      // create new edges between vertices[pos] and its remote vertex if not found
      local_id = get_ent_localid(mesh, down_vtx[pos]);

      geom_ent = (gmi_ent*)(mesh->toModel(down_vtx[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
      new_down_vtx[0] = down_vtx[pos];
      new_down_vtx[1] = remote_vertices[local_id],

      edgesBtwPlane[pos] = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);

      new_id =num_local_edge*2+edge_counter;
      mesh->setIntTag(edgesBtwPlane[pos], local_entid_tag, &new_id);
      btw_plane_edges.push_back(edgesBtwPlane[pos]);
      ++edge_counter;
    }// for (int pos=0; pos<3; ++pos)

    /**create quads between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the newly created face btw planes*/
      quad_faces[pos]=NULL;
      edge_faces.clear();

      num_upward = mesh->countUpward(down_edge[pos]);

      for (int i=0; i<num_upward; ++i)
        edge_faces.push_back(mesh->getUpward(down_edge[pos],i));

      for (unsigned int i=0; i<edge_faces.size(); ++i)
      {
        local_id = get_ent_localid(mesh, edge_faces[i]);
        if (local_id>=num_local_face) // face is between planes
        {
          quad_faces[pos]=edge_faces[i];
          break; // get out of for loop
        }
      }
      if (quad_faces[pos]!=NULL) continue;

      /**create new quad between two planes if found*/
      Downward quad_edges;
      quad_edges[0]=down_edge[pos];
      Downward down_edge_vtx;
      mesh->getDownward(down_edge[pos], 0, down_edge_vtx);
      MeshEntity* vtx_1=down_edge_vtx[1];
      for (int i=0; i<3; ++i)
      { 
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_1 || edgeBtw_down[1]==vtx_1)
        {
          quad_edges[1]=edgesBtwPlane[i];
          break; // get out of for loop
        }
      }

      quad_edges[2]=edgesNextPlane[pos];
      MeshEntity* vtx_2 = down_edge_vtx[0];

      for (int i=0; i<3; ++i)
      {
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_2 || edgeBtw_down[1]==vtx_2)
        {
          quad_edges[3]=edgesBtwPlane[i];
          break;// get out of for loop
        }
      }
      geom_ent = (gmi_ent*)(mesh->toModel(down_edge[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);

      quad_faces[pos]= mesh->createEntity(apf::Mesh::QUAD, (ModelEntity*)new_geom_ent, quad_edges);
      btw_plane_faces.push_back(quad_faces[pos]);
      new_id = num_local_face*2+face_counter;
      mesh->setIntTag(quad_faces[pos], local_entid_tag, &new_id);
      ++face_counter;
    } //     for (int pos=0;pos<3; ++pos)

    // create regions per face on local plane
    wedge_faces[0]=e;
    wedge_faces[1]=quad_faces[0];
    wedge_faces[2]=quad_faces[1];
    wedge_faces[3]=quad_faces[2];
    local_id = get_ent_localid(mesh, e);
    wedge_faces[4]=remote_faces[local_id];

    if (flip_wedge) // flip top & bottom of wedge to avoid negative volume in the last plane
    {
      wedge_faces[4]=e;
      wedge_faces[0]=remote_faces[local_id];
    }

    /**create new region between two planes*/
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
//    std::cout<<"[M3D-C1 INFO] (p"<<PCU_Comm_Self()<<") create PRISM with face "<<get_ent_localid(mesh, wedge_faces[0])<<"(t="<< mesh->getType(wedge_faces[0])<<"), "<<get_ent_localid(mesh, wedge_faces[1])<<"(t="<< mesh->getType(wedge_faces[1])<<"), "<<get_ent_localid(mesh, wedge_faces[2])<<"(t="<< mesh->getType(wedge_faces[2])<<"), "<<get_ent_localid(mesh, wedge_faces[3])<<"(t="<< mesh->getType(wedge_faces[3])<<"), "<<get_ent_localid(mesh, wedge_faces[4])<<"(t="<< mesh->getType(wedge_faces[4])<<") (#face="<<mesh->count(2)<<")"<<std::endl;
    new_ent = mesh->createEntity(apf::Mesh::PRISM, (ModelEntity*)new_geom_ent, wedge_faces);
    btw_plane_regions.push_back(new_ent);
    mesh->setIntTag(new_ent, local_entid_tag, &rgn_counter);
    ++rgn_counter;
  }
}

// this routine works only with the initial mesh where local entity ID is consecutive
void copy_master_plane_obsolete(apf::Mesh2* mesh, MeshTag* local_entid_tag, 
     MeshEntity** remote_vertices, MeshEntity** remote_edges, MeshEntity** remote_faces, 
     std::vector<MeshEntity*>& btw_plane_edges, std::vector<MeshEntity*>& btw_plane_faces, 
     std::vector<MeshEntity*>& btw_plane_regions)
{
  int local_partid=PCU_Comm_Self();
  apf::MeshIterator* ent_it;
 
  MeshEntity* e;
  MeshEntity* new_ent;
  gmi_ent* geom_ent; 

  int index, new_id, local_id, local_id_0, local_id_1, local_id_2;
  double local_plane_phi = m3dc1_model::instance()->get_phi(local_partid);
  double next_plane_phi = m3dc1_model::instance()->get_phi(m3dc1_model::instance()->next_plane_partid);
  bool flip_wedge=false;
  if (local_plane_phi>next_plane_phi) // last plane
    flip_wedge=true;

  // create remote copy of vertices on next plane
  gmi_ent* new_geom_ent = NULL;
  int num_local_vtx=mesh->count(0);
  Vector3 cur_coord;
  Vector3 new_coord;

  index=0;
  ent_it = mesh->begin(0);
  while ((e = mesh->iterate(ent_it)))
 // for (index=0;index<num_local_vtx;++index)
  {
    mesh->getPoint(e, 0, cur_coord);
    new_coord[0]=cur_coord[0];
    new_coord[1]=cur_coord[1];
    new_coord[2]=next_plane_phi;
    Vector3 param(0,0,0);     
    mesh->getParam(e,param);
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);
    new_ent = mesh->createVertex((ModelEntity*)new_geom_ent, new_coord, param);
    new_id = num_local_vtx+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_vertices[index]=new_ent;
    // update the z-coord of the current coordinate based on the prev plane partid
    cur_coord[2] = local_plane_phi;
    mesh->setPoint(e, 0, cur_coord);
    ++index;
    if (index==num_local_vtx) break;
  }
  mesh->end(ent_it);

  // create remote copy of edges on next plane
  int num_local_edge=mesh->count(1);
  Downward down_vtx;
  Downward new_down_vtx;
 
  index=0;
  ent_it = mesh->begin(1);
  while ((e = mesh->iterate(ent_it)))
  //for (index=0;index<num_local_edge;++index)
  {  
    mesh->getDownward(e, 0, down_vtx);
    local_id_0 = get_ent_localid(mesh, down_vtx[0]);
    local_id_1 = get_ent_localid(mesh, down_vtx[1]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_down_vtx[0] = remote_vertices[local_id_0];
    new_down_vtx[1] = remote_vertices[local_id_1];
    new_ent = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);

    new_id = num_local_edge+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_edges[index]=new_ent;
    ++index;
    if (index==num_local_edge) break;
  }
  mesh->end(ent_it);

  // create remote copy of faces on next plane
  int num_local_face=mesh->count(2);
  Downward down_edge;
  Downward new_down_edge;
  index=0;
  ent_it = mesh->begin(2);
  while ((e = mesh->iterate(ent_it)))
  {  
    mesh->getDownward(e, 1, down_edge);
    local_id_0 = get_ent_localid(mesh, down_edge[0]);
    local_id_1 = get_ent_localid(mesh, down_edge[1]);
    local_id_2 = get_ent_localid(mesh, down_edge[2]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_down_edge[0] = remote_edges[local_id_0];
    new_down_edge[1] = remote_edges[local_id_1];
    new_down_edge[2] = remote_edges[local_id_2];

    new_ent = mesh->createEntity(apf::Mesh::TRIANGLE, (ModelEntity*)new_geom_ent, new_down_edge);
    new_id = num_local_face+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_faces[index]=new_ent;
    ++index;
    if (index==num_local_face) break;
  }
  mesh->end(ent_it);

  // create edges, faces and regions between planes. 
  // if edges are classified on geom edge, a new face is classified on geom_surface_face
  // otherwise, a new face is classified on geom_region

  Downward quad_faces;
  Downward wedge_faces;

  int edge_counter=0, face_counter=0, rgn_counter=0;
  std::vector<MeshEntity*> vertex_edges;
  std::vector<MeshEntity*> edge_faces;

  Downward edgesNextPlane;
  Downward edgesBtwPlane;
  int num_upward;

  index=0;
  ent_it = mesh->begin(2);
  while ((e = mesh->iterate(ent_it)))
  {  
    mesh->getDownward(e, 0, down_vtx);
    mesh->getDownward(e, 1, down_edge);

    for (int pos=0; pos<3; ++pos)
    {
      local_id = get_ent_localid(mesh, down_edge[pos]);
      edgesNextPlane[pos] = remote_edges[local_id];
    }

    /**create edges between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the new created face btw planes*/
      edgesBtwPlane[pos]=NULL;
      vertex_edges.clear();
      num_upward = mesh->countUpward(down_vtx[pos]);
      for (int i=0; i<num_upward; ++i)
        vertex_edges.push_back(mesh->getUpward(down_vtx[pos], i));

      for (unsigned int i=0; i<vertex_edges.size();++i)
      {
        // get the local id
        local_id = get_ent_localid(mesh, vertex_edges[i]);
        if (local_id>=num_local_edge) // edge is between planes
        {
          edgesBtwPlane[pos]=vertex_edges[i];
          break; // get out of for loop
        }
      }

      if (edgesBtwPlane[pos]!=NULL) continue;

      // create new edges between vertices[pos] and its remote vertex if not found
      local_id = get_ent_localid(mesh, down_vtx[pos]);

      geom_ent = (gmi_ent*)(mesh->toModel(down_vtx[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
      new_down_vtx[0] = down_vtx[pos];
      new_down_vtx[1] = remote_vertices[local_id],

      edgesBtwPlane[pos] = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);

      new_id =num_local_edge*2+edge_counter;
      mesh->setIntTag(edgesBtwPlane[pos], local_entid_tag, &new_id);
      btw_plane_edges.push_back(edgesBtwPlane[pos]);
      ++edge_counter;
    }// for (int pos=0; pos<3; ++pos)

    /**create quads between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the newly created face btw planes*/
      quad_faces[pos]=NULL;
      edge_faces.clear();
      num_upward = mesh->countUpward(down_edge[pos]);

      for (int i=0; i<num_upward; ++i)
        edge_faces.push_back(mesh->getUpward(down_edge[pos],i));

      for (unsigned int i=0; i<edge_faces.size(); ++i)
      {
        local_id = get_ent_localid(mesh, edge_faces[i]);
        if (local_id>=num_local_face) // face is between planes
        {
          quad_faces[pos]=edge_faces[i];
          break; // get out of for loop
        }
      }
      if (quad_faces[pos]!=NULL) continue;

      /**create new quad between two planes if found*/
      Downward quad_edges;
      quad_edges[0]=down_edge[pos];
      Downward down_edge_vtx;
      mesh->getDownward(down_edge[pos], 0, down_edge_vtx);
      MeshEntity* vtx_1=down_edge_vtx[1];
      for (int i=0; i<3; ++i)
      { 
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_1 || edgeBtw_down[1]==vtx_1)
        {
          quad_edges[1]=edgesBtwPlane[i];
          break; // get out of for loop
        }
      }

      quad_edges[2]=edgesNextPlane[pos];
      MeshEntity* vtx_2 = down_edge_vtx[0];

      for (int i=0; i<3; ++i)
      {
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_2 || edgeBtw_down[1]==vtx_2)
        {
          quad_edges[3]=edgesBtwPlane[i];
          break;// get out of for loop
        }
      }
      geom_ent = (gmi_ent*)(mesh->toModel(down_edge[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);

      quad_faces[pos]= mesh->createEntity(apf::Mesh::QUAD, (ModelEntity*)new_geom_ent, quad_edges);
      btw_plane_faces.push_back(quad_faces[pos]);
      new_id = num_local_face*2+face_counter;
      mesh->setIntTag(quad_faces[pos], local_entid_tag, &new_id);
      ++face_counter;
    } //     for (int pos=0;pos<3; ++pos)

    // create regions per face on local plane
    wedge_faces[0]=e;
    wedge_faces[1]=quad_faces[0];
    wedge_faces[2]=quad_faces[1];
    wedge_faces[3]=quad_faces[2];
    local_id = get_ent_localid(mesh, e);
    wedge_faces[4]=remote_faces[local_id];

    if (flip_wedge) // flip top & bottom of wedge to avoid negative volume in the last plane
    {
      wedge_faces[4]=e;
      wedge_faces[0]=remote_faces[local_id];
    }

    /**create new region between two planes*/
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
//    std::cout<<"[M3D-C1 INFO] (p"<<PCU_Comm_Self()<<") create PRISM with face "<<get_ent_localid(mesh, wedge_faces[0])<<"(t="<< mesh->getType(wedge_faces[0])<<"), "<<get_ent_localid(mesh, wedge_faces[1])<<"(t="<< mesh->getType(wedge_faces[1])<<"), "<<get_ent_localid(mesh, wedge_faces[2])<<"(t="<< mesh->getType(wedge_faces[2])<<"), "<<get_ent_localid(mesh, wedge_faces[3])<<"(t="<< mesh->getType(wedge_faces[3])<<"), "<<get_ent_localid(mesh, wedge_faces[4])<<"(t="<< mesh->getType(wedge_faces[4])<<") (#face="<<mesh->count(2)<<")"<<std::endl;
    new_ent = mesh->createEntity(apf::Mesh::PRISM, (ModelEntity*)new_geom_ent, wedge_faces);
    btw_plane_regions.push_back(new_ent);
    mesh->setIntTag(new_ent, local_entid_tag, &rgn_counter);
    ++rgn_counter;
    ++index;
    if (index==num_local_face) break;
  }
  mesh->end(ent_it);
}

void create_localid(apf::Mesh2* mesh, int dim)
{
  MeshIterator* ent_it = mesh->begin(dim);
  apf::MeshEntity* e;

  // assign sequential local ID starting with 0
  int index=0;
  while ((e = mesh->iterate(ent_it)))
  {
    mesh->setIntTag(e, m3dc1_mesh::instance()->local_entid_tag, &index);
    ++index;
  }
  mesh->end(ent_it);
}

void exchange_copies(apf::Mesh2* mesh,
  std::vector<MeshEntity*>* vertices, std::vector<MeshEntity*>* edges,
  std::vector<MeshEntity*>* faces, 
  MeshEntity** remote_vertices, MeshEntity** remote_edges, MeshEntity** remote_faces,
  std::map<MeshEntity*, MeshEntity*>& partbdry_ent_map,
  std::vector<MeshEntity*>& ent_vec)
{
  int num_local_vtx = vertices->size();
  int num_local_edge = edges->size();
  int num_local_face = faces->size();
  MeshEntity* e;

  PCU_Comm_Begin();
  int num_entities = num_local_vtx+num_local_edge+num_local_face;
  MeshEntity** entities = new MeshEntity*[num_entities];
  int pos=0;
  for (int index=0; index<num_local_vtx; ++index)
    entities[pos++]=remote_vertices[index];
  for (int index=0; index<num_local_edge; ++index)
    entities[pos++]=remote_edges[index];

  for (int index=0; index<num_local_face; ++index)
    entities[pos++]=remote_faces[index];

  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &num_entities, sizeof(int));
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &(entities[0]),num_entities*sizeof(MeshEntity*));
  PCU_Comm_Send();
  delete [] entities;

  // receive phase begins
  std::vector<MeshEntity*>::iterator eit;
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      int num_ent;
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      MeshEntity** s_ent = new MeshEntity*[num_ent];
      PCU_Comm_Unpack(&(s_ent[0]), num_ent*sizeof(MeshEntity*));
      pos=0;
      for (eit=vertices->begin(); eit!=vertices->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }

      for (eit=edges->begin(); eit!=edges->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }

      for (eit=faces->begin(); eit!=faces->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        ent_vec.push_back(e);
        ++pos;
      }
      delete [] s_ent;
    }
  }
}

// *********************************************************
void m3dc1_mesh::build3d(int num_field, int* field_id, int* num_dofs_per_value, 
                         int** ge_tag, bool initial_setup,
                         std::map<int, apf::Field*>* new_fields)
// *********************************************************
{
  int local_partid=PCU_Comm_Self();

  MeshIterator* ent_it;
  apf::MeshEntity* e;
  int index; 

  // assign uniq id to part bdry entities
  MeshTag* partbdry_id_tag = mesh->findTag("m3dc1_pbdry_globid");
  if (!partbdry_id_tag)
   partbdry_id_tag = mesh->createIntTag("m3dc1_pbdry_globid", 1);

  for (int dim=0; dim<2; ++dim)
    assign_uniq_partbdry_id(mesh, dim, partbdry_id_tag);

  // copy 2D mesh in process group 0 to other process groups
  std::map<int, MeshEntity*> partbdry_entities[2];

  std::vector<apf::MeshEntity*>* vertices = new std::vector<apf::MeshEntity*>;
  std::vector<apf::MeshEntity*>* edges = new std::vector<apf::MeshEntity*>;
  std::vector<apf::MeshEntity*>* faces = new std::vector<apf::MeshEntity*>;

  std::vector<MeshEntity*>::iterator eit;

  if (!m3dc1_model::instance()->local_planeid)
  {
    ent_it = mesh->begin(0);
    while ((e = mesh->iterate(ent_it)))
      vertices->push_back(e);
    mesh->end(ent_it);

    ent_it = mesh->begin(1);
    while ((e = mesh->iterate(ent_it)))
      edges->push_back(e);
    mesh->end(ent_it);

    ent_it = mesh->begin(2);
    while ((e = mesh->iterate(ent_it)))
      faces->push_back(e);
    mesh->end(ent_it);
  }     
  
  for (int dim=0; dim<3; ++dim)
  {
    PCU_Comm_Begin();
    sendEntities(mesh, dim, partbdry_id_tag);
    PCU_Comm_Send();
    switch (dim)
    {
      case 0: m3dc1_receiveVertices(mesh, vertices, partbdry_id_tag, partbdry_entities); break;
      case 1: m3dc1_receiveEdges(mesh, vertices, edges, partbdry_id_tag, partbdry_entities); break;
      case 2: m3dc1_receiveFaces(mesh, edges, faces); break;
      default: break;
    }   
  }

  m3dc1_stitchLink(mesh, partbdry_id_tag, partbdry_entities);

  for (int dim=0; dim<2; ++dim)
  {
    ent_it = mesh->begin(dim);
    while ((e = mesh->iterate(ent_it)))
    {
      if (!mesh->isShared(e)) continue;
      Copies remotes;
      Parts parts;
      mesh->getRemotes(e, remotes);
      APF_ITERATE(Copies, remotes, it)
        parts.insert(it->first);
      parts.insert(local_partid);
      mesh->setResidence (e, parts); // set pclassification
      mesh->removeTag(e, partbdry_id_tag);
    }
    mesh->end(ent_it);
  }
  mesh->destroyTag(partbdry_id_tag);

  // update global ent counter
  MPI_Allreduce(num_own_ent, num_global_ent, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // construct 3D model
  changeMdsDimension(mesh, 3);

  if (initial_setup)
    m3dc1_model::instance()->create3D();
  else
  {
    gmi_model* model = m3dc1_model::instance()->model;
    gmi_iter* g_it;
    gmi_ent* ge;

    for (int dim=0; dim<=2; ++dim)
    {
      g_it = gmi_begin(model, dim);
      int cnt=0;
      while ((ge = gmi_next(model, g_it)))
        gmi_base_set_tag (model, ge, ge_tag[dim][cnt++]);
      gmi_end(model, g_it);
    }
  }

  int num_local_vtx = mesh->count(0);
  int num_local_edge = mesh->count(1);
  int num_local_face = mesh->count(2);

  MeshEntity** remote_vertices=new MeshEntity*[num_local_vtx];
  MeshEntity** remote_edges=new MeshEntity*[num_local_edge];
  MeshEntity** remote_faces=new MeshEntity*[num_local_face];
  std::vector<MeshEntity*> btw_plane_edges;
  std::vector<MeshEntity*> btw_plane_faces;
  std::vector<MeshEntity*> btw_plane_regions;

  // temporary keep 2 versions not to break existing codes
  if (initial_setup)
    copy_master_plane_obsolete(mesh, local_entid_tag, remote_vertices, remote_edges, 
      remote_faces, btw_plane_edges, btw_plane_faces, btw_plane_regions);
  else
    copy_master_plane(mesh, local_entid_tag, 
      vertices, edges, faces,
      remote_vertices, remote_edges,
      remote_faces, btw_plane_edges, btw_plane_faces, btw_plane_regions);

  // exchange remote copies to set remote copy links
  std::vector<MeshEntity*> ent_vec;
  std::map<MeshEntity*, MeshEntity*> partbdry_ent_map;
 
  PCU_Comm_Begin();
  int num_entities = num_local_vtx+num_local_edge+num_local_face;

  MeshEntity** entities = new MeshEntity*[num_entities];
  // FIXME: use memcpy
/*
  memcpy(&(entities[0]),&(remote_vertices[0]), num_local_vtx*sizeof(int));
  memcpy(&(entities[num_local_vtx]), &(remote_edges[0]), num_local_edge*sizeof(int));
  memcpy(&(entities[num_local_vtx+num_local_edge]), &(remote_faces[0]), num_local_face*sizeof(int));
*/
  int pos=0;
  for (int i=0; i<num_local_vtx; ++i)
    entities[pos++]=remote_vertices[i];
  for (int i=0; i<num_local_edge; ++i)
    entities[pos++]=remote_edges[i];
  for (int i=0; i<num_local_face; ++i)
    entities[pos++]=remote_faces[i];
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &num_entities, sizeof(int));
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &(entities[0]),num_entities*sizeof(MeshEntity*));
  PCU_Comm_Send();
  delete [] entities;

  // receive phase begins
  while (PCU_Comm_Listen())
  {
    int index=0;
    while (!PCU_Comm_Unpacked())
    {
      int num_ent;
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      MeshEntity** s_ent = new MeshEntity*[num_ent];
      PCU_Comm_Unpack(&(s_ent[0]), num_ent*sizeof(MeshEntity*));
      pos=0;
      for (eit=vertices->begin(); eit!=vertices->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }

      for (eit=edges->begin(); eit!=edges->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }

      for (eit=faces->begin(); eit!=faces->end(); ++eit)
      {
        e = *eit;
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        ent_vec.push_back(e);
        ++pos;
      }
      delete [] s_ent;
    }
  }

  push_new_entities(mesh, partbdry_ent_map);

  bounce_orig_entities(mesh, ent_vec, m3dc1_model::instance()->prev_plane_partid);

  // update partition classification
  for (eit=vertices->begin(); eit!=vertices->end(); ++eit)
  {
    e = *eit;
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (index=0; index<num_local_vtx; ++index)
  {
    e = remote_vertices[index];
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (eit=edges->begin(); eit!=edges->end(); ++eit)
  {
    e = *eit;
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first); 
    parts.insert(local_partid);  
    mesh->setResidence(e, parts);
  }

  for (index=0; index<num_local_edge; ++index)
  {
    e = remote_edges[index];
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (eit=faces->begin(); eit!=faces->end(); ++eit)
  {
    e = *eit;
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (index=0; index<num_local_face; ++index)
  {
    e = remote_faces[index];
    Copies remotes;
    Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  // update partition model
  mesh->acceptChanges();

  // update m3dc1_mesh internal data (# local, owned and global entities)
  update_partbdry(remote_vertices, remote_edges, remote_faces, btw_plane_edges, btw_plane_faces, btw_plane_regions);

  // delete existing local numbering
  apf::Numbering* local_n = mesh->findNumbering(mesh->getShape()->getName());
  if (local_n) destroyNumbering(local_n);
  
  // re-create the field and copy field data on master process group to non-master
  for (int i=0; i<num_field; ++i)
    if (!initial_setup && m3dc1_model::instance()->num_plane)
      update_field_3d(field_id[i], new_fields, num_local_vtx, vertices, remote_vertices);
    else
      update_field(field_id[i], num_dofs_per_value[i], num_local_vtx, remote_vertices);

  // remove vectors
  vertices->clear();
  edges->clear();
  faces->clear();
  delete vertices;
  delete edges;
  delete faces;

  // clear temp memory
  delete [] remote_vertices;
  delete [] remote_edges;
  delete [] remote_faces;
  set_node_adj_tag();

  if (!PCU_Comm_Self())
    std::cout<<"\n*** MESH SWITCHED TO 3D ***\n\n";
}

// *********************************************************
void m3dc1_mesh::initialize(bool update_adj)
// *********************************************************
{
  if (!local_entid_tag) local_entid_tag = mesh->createIntTag("m3dc1_local_ent_id", 1);
  if (!own_partid_tag) own_partid_tag = mesh->createIntTag("m3dc1_own_part_id", 1);
  if (!num_global_adj_node_tag) num_global_adj_node_tag = mesh->createIntTag("m3dc1_num_global_adj_node", 1);
  if (!num_own_adj_node_tag) num_own_adj_node_tag = mesh->createIntTag("m3dc1_num_own_adj_node", 1);

  reset();
  MeshEntity* e;

  int counter = 0, own_partid, local_partid=PCU_Comm_Self();

  for (int d=0; d<4; ++d)
  {
    num_local_ent[d] = mesh->count(d);
    counter=0;
    MeshIterator* it = mesh->begin(d);
    while ((e = mesh->iterate(it)))
    {
      mesh->setIntTag(e, local_entid_tag, &counter);
      own_partid = mesh->getOwner(e);
      mesh->setIntTag(e, own_partid_tag, &own_partid);
      if (own_partid==local_partid)
        ++num_own_ent[d];
      ++counter;
    }
    mesh->end(it);
  }

  MPI_Allreduce(num_own_ent, num_global_ent, 4, MPI_INT, MPI_SUM, PCU_Get_Comm());
  if (update_adj) set_node_adj_tag();
}


// *********************************************************
void m3dc1_mesh::update_partbdry(MeshEntity** remote_vertices, MeshEntity** remote_edges, MeshEntity** remote_faces,
    std::vector<MeshEntity*>& btw_plane_edges, std::vector<MeshEntity*>& btw_plane_faces, std::vector<MeshEntity*>& btw_plane_regions)
// *********************************************************
{
  int num_orig_vtx = num_local_ent[0];
  int num_orig_edge = num_local_ent[1];
  int num_orig_face = num_local_ent[2];

  num_local_ent[0] = num_orig_vtx*2;
  num_local_ent[1] = num_orig_edge*2+btw_plane_edges.size();
  num_local_ent[2] = num_orig_face*2+btw_plane_faces.size();
  num_local_ent[3] = btw_plane_regions.size();

  //num_own_ent[0] = num_orig_vtx;
  num_own_ent[1] += btw_plane_edges.size();
  num_own_ent[2] += btw_plane_faces.size();
  num_own_ent[3] = btw_plane_regions.size();


}

// upon  mesh modification, update field wrt memory for dof data and numbering 
// and rebuild the existing dof data in updated field.
// *********************************************************
void update_field_3d (int field_id,
                      std::map<int, apf::Field*>* new_fields, 
                      int num_2d_vtx, 
                      std::vector<MeshEntity*>* vertices, MeshEntity** remote_vertices)
// *********************************************************
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  int num_rank_plane = m3dc1_model::instance()->group_size; 
  MeshIterator* ent_it;
  apf::MeshEntity* e;

  int num_rank = PCU_Comm_Peers();
  int local_rank = PCU_Comm_Self();

  // get the field info and save it for later creation
  char f_name[100];
  int num_values, scalar_type, m3dc1_ndof;   

  apf::Field* field;
 
  m3dc1_field_getinfo (&field_id, f_name, &num_values, &scalar_type, &m3dc1_ndof);
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[field_id]->get_field();
  int to_pid, num_dof = countComponents(f); 
  double* dof_data = new double[num_dof*num_2d_vtx]; 


  PCU_Comm_Begin();
  if (!m3dc1_model::instance()->local_planeid)
  {
    apf::freezeFields(mesh); // turning field data from tag to array
    to_pid = local_rank+num_rank_plane;
    while (to_pid<num_rank)
    {
      PCU_Comm_Pack(to_pid, &num_2d_vtx, sizeof(int));
      field = (*new_fields)[num_rank*field_id+to_pid]; 
      if (!field)
        std::cout<<"[p"<<PCU_Comm_Self()<<"] cannot find new_field "<<num_rank*field_id+to_pid<<"\n";
      else
        std::cout<<"[p"<<PCU_Comm_Self()<<"] sending field "
                 <<field_id<<"(new_fields["<<num_rank*field_id+to_pid<<"]) to p"<<to_pid<<"\n";
      assert(field);
      memcpy(&(dof_data[0]), apf::getArrayData(field), num_dof*num_2d_vtx*sizeof(double));
      PCU_Comm_Pack(to_pid, &(dof_data[0]), num_dof*num_2d_vtx*sizeof(double));
      to_pid += num_rank_plane;
    }
  }

  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_data;

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      recv_dof_data = new double[num_dof*recv_num_ent];
      PCU_Comm_Unpack(&(recv_dof_data[0]), num_dof*recv_num_ent*sizeof(double));
    }
  }

  apf::unfreezeFields(mesh); // turning field data from array to tag

  // delete the field
  m3dc1_field_delete(&field_id);

  // create the field with same attributes
  m3dc1_field_create (&field_id, f_name, &num_values, &scalar_type, &m3dc1_ndof);

  // re-construct the dof data
  f = (*(m3dc1_mesh::instance()->field_container))[field_id]->get_field();
  int new_total_ndof = apf::countComponents(f);
  double* dof_val = new double[new_total_ndof];
  if (!m3dc1_model::instance()->local_planeid)
  {
    for (int vid=0; vid<num_2d_vtx; ++vid)
    {
      e = vertices->at(vid);
      for (int i=0; i<num_dof; ++i)
        dof_val[i] = dof_data[vid*num_dof+i];
      if (new_total_ndof>num_dof)
        for (int i=num_dof; i<new_total_ndof; ++i)
          dof_val[i] = 0.0;
      setComponents(f, e, 0, dof_val);
      setComponents(f, remote_vertices[vid++], 0, dof_val);
    }
  }
  else
  {
    for (int vid=0; vid<recv_num_ent; ++vid)
    {
      e = vertices->at(vid);
      for (int i=0; i<num_dof; ++i)
         dof_val[i] = recv_dof_data[vid*num_dof+i];
      if (new_total_ndof>num_dof)
        for (int i=num_dof; i<new_total_ndof; ++i)
          dof_val[i] = 0.0;
      setComponents(f, e, 0, dof_val);
      setComponents(f, remote_vertices[vid++], 0, dof_val);
    } // vid
    delete [] recv_dof_data;     
  } // while

  delete [] dof_val;
  delete [] dof_data;
}

// upon  mesh modification, update field wrt memory for dof data and numbering 
// and rebuild the existing dof data in updated field.
// *********************************************************
void update_field (int field_id, int ndof_per_value, int num_2d_vtx, MeshEntity** remote_vertices)
// *********************************************************
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshIterator* ent_it;
  apf::MeshEntity* e;

  // get the field info and save it for later creation
  char f_name[100];
  int num_values, scalar_type, old_numdof;   
  m3dc1_field_getinfo (&field_id, f_name, &num_values, &scalar_type, &old_numdof);

  // copy the existing dof data
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[field_id]->get_field();
  int old_total_ndof = apf::countComponents(f);
  double* dof_val; 

  freeze(f); // switch dof data from tag to array
  
  PCU_Comm_Begin();

  if (!m3dc1_model::instance()->local_planeid)
  {
    dof_val = new double[old_total_ndof*num_2d_vtx];
    memcpy(&(dof_val[0]), apf::getArrayData(f), old_total_ndof*num_2d_vtx*sizeof(double));
    int proc=PCU_Comm_Self()+m3dc1_model::instance()->group_size;
    while (proc<PCU_Comm_Peers())
    {       
      PCU_Comm_Pack(proc, &num_2d_vtx, sizeof(int));
      PCU_Comm_Pack(proc, &(dof_val[0]), old_total_ndof*num_2d_vtx*sizeof(double));
      proc+=m3dc1_model::instance()->group_size;
    }
  }
  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_val;

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked())
    { 
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      recv_dof_val = new double[old_total_ndof*recv_num_ent];
      PCU_Comm_Unpack(&(recv_dof_val[0]), old_total_ndof*recv_num_ent*sizeof(double));
    }
  }
  // delete the field
  m3dc1_field_delete(&field_id);

  // create the field with same attributes
  m3dc1_field_create (&field_id, f_name, &num_values, &scalar_type, &ndof_per_value);
  
  // re-construct the dof data
  f = (*(m3dc1_mesh::instance()->field_container))[field_id]->get_field();
  int new_total_ndof = apf::countComponents(f);
  double* dof_data = new double[new_total_ndof];
  if (!m3dc1_model::instance()->local_planeid)
  {
    int index=0;
    ent_it = mesh->begin(0);
    while ((e = mesh->iterate(ent_it)))
    // for (int index=0; index<num_2d_vtx; ++index)
    {
      for (int i=0; i<old_total_ndof; ++i)
        dof_data[i] = dof_val[index*old_total_ndof+i];
      if (new_total_ndof>old_total_ndof)
        for (int i=old_total_ndof; i<new_total_ndof; ++i)
          dof_data[i] = 0.0;
      setComponents(f, e, 0, dof_data);
      setComponents(f, remote_vertices[index++], 0, dof_data);
    }
    mesh->end(ent_it);
    delete [] dof_val;
  }
  else
  {
    int index=0;
    ent_it = mesh->begin(0);
    while ((e = mesh->iterate(ent_it)))
    //for (int index=0; index<recv_num_ent; ++index)
    {
      for (int i=0; i<old_total_ndof; ++i)
         dof_data[i] = recv_dof_val[index*old_total_ndof+i];
      if (new_total_ndof>old_total_ndof)
        for (int i=old_total_ndof; i<new_total_ndof; ++i)
          dof_data[i] = 0.0;
      setComponents(f, e, 0, dof_data);
      setComponents(f, remote_vertices[index++], 0, dof_data);
    } // index
    mesh->end(ent_it);
    delete [] recv_dof_val;     
  } // while

  delete [] dof_data;
}

 struct entMsg
  {
    int pid;
    MeshEntity* ent;
    entMsg( int pid_p=0, MeshEntity* ent_p=NULL)
    {
      pid=pid_p;
      ent=ent_p;
    }
  };
  struct classcomp
  {
    bool operator() (const entMsg& lhs, const entMsg& rhs) const
    {
      if (lhs.ent==rhs.ent) return lhs.pid<rhs.pid;
      else return lhs.ent<rhs.ent;
    }
  };

// FIXME: this crashes when called with ghosted mesh
// **********************************************
void m3dc1_mesh::set_node_adj_tag()
// **********************************************
{
  int value;
  int brgType = (num_local_ent[3])?3:2;

  apf::MeshEntity* e;
  apf::MeshIterator* it = mesh->begin(0);
  PCU_Comm_Begin();
  while ((e = mesh->iterate(it)))
  {
    int num_adj_node=0;
    Adjacent elements;
    getBridgeAdjacent(mesh, e, brgType, 0, elements);
    int num_adj = elements.getSize();

    for (int i=0; i<num_adj; i++)
    {
      if (is_ent_original(mesh, elements[i]))
        ++num_adj_node;
    }
    mesh->setIntTag(e, num_own_adj_node_tag, &num_adj_node);

    if (!mesh->isShared(e)) continue;
    // first pass msg size to owner
    int own_partid = get_ent_ownpartid(mesh, e);
    MeshEntity* own_copy = get_ent_owncopy(mesh, e);

    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_copy);
    PCU_Comm_Pack(own_partid, &num_adj,sizeof(int));
  }
  mesh->end(it);

  PCU_Comm_Send();

  std::map<apf::MeshEntity*, std::map<int, int> > count_map;
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      PCU_Comm_Unpack(&value,sizeof(int));
      count_map[e][PCU_Comm_Sender()]=value;
    }
  }
  
  // pass entities to ownner
  std::map<apf::MeshEntity*, std::set<entMsg, classcomp> > count_map2;
  it = mesh->begin(0);
  PCU_Comm_Begin();
  while ((e = mesh->iterate(it)))
  {
    // pass entities to ownner

    std::vector<entMsg> msgs;
    Adjacent elements;
    getBridgeAdjacent(mesh, e, brgType, 0, elements);

    MeshEntity* ownerEnt=get_ent_owncopy(mesh, e);
    int own_partid = get_ent_ownpartid(mesh, e);
    for (int i=0; i<elements.getSize(); i++)
    {
      MeshEntity* ownerEnt2=get_ent_owncopy(mesh, elements[i]);
      int owner=get_ent_ownpartid(mesh, elements[i]);
      msgs.push_back(entMsg(owner, ownerEnt2));
      if (own_partid==PCU_Comm_Self()) 
      {
        count_map2[e].insert(*msgs.rbegin());
      }
    }

    if (own_partid!=PCU_Comm_Self())
    {
      PCU_COMM_PACK(own_partid, ownerEnt);
      PCU_Comm_Pack(own_partid, &msgs.at(0),sizeof(entMsg)*msgs.size());
    }
  }
  mesh->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      int sizeData = count_map[e][PCU_Comm_Sender()];
      std::vector<entMsg> data(sizeData);
      PCU_Comm_Unpack(&data.at(0),sizeof(entMsg)*sizeData);
      for (int i=0; i<data.size(); i++)
      {
        count_map2[e].insert(data.at(i));
      }
    }
  }

  for (std::map<apf::MeshEntity*, std::set<entMsg,classcomp> >::iterator mit=count_map2.begin(); mit!=count_map2.end(); ++mit)
  {
    e = mit->first;
    int num_global_adj =count_map2[e].size();
    mesh->setIntTag(mit->first, num_global_adj_node_tag, &num_global_adj);
  }
}

// **********************************************
//void m3dc1_mesh::set_node_adj_tag()
// **********************************************
/*
{
  int num_adj_node, value, own_partid;
  int brgType = (num_local_ent[3])?3:2;

  MeshEntity* e;
  MeshEntity* own_copy;

  std::map<MeshEntity*, std::set<MeshEntity*> > copy_map;

  PCU_Comm_Begin();
  MeshIterator* it = mesh->begin(0);
  while ((e = mesh->iterate(it)))
  {
    num_adj_node=0;
    Adjacent elements;
    getBridgeAdjacent(mesh, e, brgType, 0, elements);
    int num_adj = elements.getSize();
    MeshEntity** adj_copy=new MeshEntity*[num_adj];

    for (int i=0; i<num_adj; ++i)
    {
      if (is_ent_original(mesh, elements[i])) 
         ++num_adj_node;         
      adj_copy[i] = get_ent_owncopy(mesh, elements[i]);
      if (is_ent_original(mesh, e)) 
        copy_map[e].insert(adj_copy[i]);  // copy_map[e] @ owner_part contains owner copies of adj entities
    }
    mesh->setIntTag(e, num_own_adj_node_tag, &num_adj_node);

    if (!is_ent_original(mesh, e))
    {
      own_partid = get_ent_ownpartid(mesh, e);
      own_copy = get_ent_owncopy(mesh, e);
      PCU_COMM_PACK(own_partid, own_copy);
      PCU_Comm_Pack(own_partid, &num_adj, sizeof(int));
      PCU_Comm_Pack(own_partid, &adj_copy[0], sizeof(MeshEntity*)*num_adj);
    }
    delete [] adj_copy;
  }
  mesh->end(it);

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      int num_adj;
      MeshEntity* own_e;
      PCU_COMM_UNPACK(own_e);
      PCU_Comm_Unpack(&num_adj,sizeof(int));
      MeshEntity** adj_ent = new MeshEntity*[num_adj];
      PCU_Comm_Unpack(&adj_ent[0],sizeof(MeshEntity*)*num_adj);
      for (int i=0; i<num_adj; i++)
        copy_map[own_e].insert(adj_ent[i]); 
      delete [] adj_ent;
    }
  }

  int num_global_adj; 
  for (std::map<MeshEntity*, std::set<MeshEntity*> >::iterator mit=copy_map.begin(); mit!=copy_map.end(); ++mit)
  { 
    num_global_adj = copy_map[mit->first].size();
    mesh->setIntTag(mit->first, num_global_adj_node_tag, &num_global_adj);
  }
}
*/

// *********************************************************
void m3dc1_mesh::print(int LINE)
// *********************************************************
{
  for (int i=0; i<=3; ++i)
    assert(mesh->count(i) == num_local_ent[i]);
         
  if (LINE<0) LINE=0;
  if (!PCU_Comm_Self())
    std::cout<<"\n[M3D-C1 Mesh Count] "<<LINE<<" #ent (VEFR): \n";

  for (int i=0; i<PCU_Comm_Peers(); ++i)
  {
    if (PCU_Comm_Self()==i)
    {
      std::cout<<"  plane "<<m3dc1_model::instance()->local_planeid<<" proc "<<PCU_Comm_Self()
             <<": #local: ("<<num_local_ent[0]<<", "<< num_local_ent[1]<<", "<<num_local_ent[2]<<", "<<num_local_ent[3]
             <<"), #own: ("<<num_own_ent[0]<<", "<< num_own_ent[1]<<", "<<num_own_ent[2]<<", "<<num_own_ent[3]
             <<")\n";
    }  
    MPI_Barrier(PCU_Get_Comm());
  }
  if (PCU_Comm_Self()==PCU_Comm_Peers()-1)
    std::cout<<"  #global: ("<<num_global_ent[0]<<", "<< num_global_ent[1]<<", "<<num_global_ent[2]<<", "<<num_global_ent[3]<<")\n";

/*  // verify the own ent counter 
  pPartEntIter vtx_iter;
  MeshEntity* ent;
  int counter=0;
  PUMI_PartEntIter_Init(mesh->getPart(0), 0, PUMI_ALLTOPO, vtx_iter);
  while (PUMI_PartEntIter_GetNext(vtx_iter, ent)==PUMI_SUCCESS)
  {
    if (get_ent_ownpartid (mesh->getPart(0), ent) == PCU_Comm_Self())
      ++counter;
  }
  assert(counter==num_own_ent[0]);
  PUMI_PartEntIter_Del(vtx_iter);

  pPartEntIter edge_iter;
  counter=0;
  PUMI_PartEntIter_Init(mesh->getPart(0), 1, PUMI_ALLTOPO, edge_iter);
  while (PUMI_PartEntIter_GetNext(edge_iter, ent)==PUMI_SUCCESS)
  {
    if (get_ent_ownpartid (mesh->getPart(0), ent) == PCU_Comm_Self())
      ++counter;
  }
  assert(counter==num_own_ent[1]);
  PUMI_PartEntIter_Del(edge_iter);

  pPartEntIter face_iter;
  counter=0;
  PUMI_PartEntIter_Init(mesh->getPart(0), 2, PUMI_ALLTOPO, face_iter);
  while (PUMI_PartEntIter_GetNext(face_iter, ent)==PUMI_SUCCESS)
  {
    if (get_ent_ownpartid (mesh->getPart(0), ent) == PCU_Comm_Self())
      ++counter;
  }
  assert(counter==num_own_ent[2]);
  PUMI_PartEntIter_Del(face_iter);
*/
}
