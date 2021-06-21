/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_mesh.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <PCU.h>
#include "apfMDS.h"
#include "ReducedQuinticImplicit.h"
#include "m3dc1_slnTransfer.h"
#include "apfShape.h" // getLagrange
#include "pumi.h"
#include "ma.h"
#include "gmi_base.h"
#include "apfMesh.h"

using namespace apf;

int get_local_planeid(int p)
{
  return p/m3dc1_model::instance()->group_size;
}

void update_remotes(apf::Mesh2* mesh, MeshEntity* e)
{
  int myrank=PCU_Comm_Self();

  apf::Copies remotes;
  apf::Copies new_remotes;
  apf::Parts parts;
  mesh->getRemotes(e,remotes);

  APF_ITERATE(apf::Copies,remotes,rit)
  {
    int p=rit->first;
    if (get_local_planeid(p)==m3dc1_model::instance()->local_planeid) // the same local plane
    {
      new_remotes[p]=rit->second; 
      parts.insert(p);
    }
  }
  parts.insert(myrank);
  mesh->clearRemotes(e);
  mesh->setRemotes(e,new_remotes);
  mesh->setResidence(e, parts);  // set pclassification
}

void m3dc1_mesh::remove_wedges()
{       
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshEntity* e;
  int myrank = PCU_Comm_Self();
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
  set_mcount();

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Wedges between 2D planes removed ***\n";
}

void m3dc1_mesh::create_wedges()
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshEntity* e;
  int local_partid = PCU_Comm_Self();

  // construct 3D mesh
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
  MeshEntity** remote_vertices=new MeshEntity*[num_local_vtx];
  apf::Vector3 cur_coord;
  apf::Vector3 new_coord;

  index=0;
  MeshIterator* it = mesh->begin(0);
  while ((e = mesh->iterate(it)))
  {
    mesh->getPoint(e, 0, cur_coord);
    new_coord[0]=cur_coord[0];
    new_coord[1]=cur_coord[1];
    new_coord[2]=next_plane_phi;
    apf::Vector3 param(0,0,0);     
    mesh->getParam(e,param);
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);
    new_ent = mesh->createVertex((ModelEntity*)new_geom_ent, new_coord, param);
    new_id = num_local_vtx+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_vertices[index++]=new_ent;
    // update the z-coord of the current coordinate based on the prev plane partid
    cur_coord[2] = local_plane_phi;
    mesh->setPoint(e, 0, cur_coord);
  }
  mesh->end(it);

  // create remote copy of edges on next plane
  int num_local_edge=mesh->count(1);
  Downward down_vtx;
  Downward new_down_vtx;
  MeshEntity** remote_edges=new MeshEntity*[num_local_edge];
  index=0;
   
  it = mesh->begin(1);
  while ((e = mesh->iterate(it)))
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
    remote_edges[index++]=new_ent;
  }
  mesh->end(it);

  // create remote copy of faces on next plane
  int num_local_face=mesh->count(2);
  Downward down_edge;
  Downward new_down_edge;
  MeshEntity** remote_faces=new MeshEntity*[num_local_face];
  index=0;
  it = mesh->begin(2);

  while ((e = mesh->iterate(it)))
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
    remote_faces[index++]=new_ent;
  }
  mesh->end(it);

  // create edges, faces and regions between planes. 
  // if edges are classified on geom edge, a new face is classified on geom_surface_face
  // otherwise, a new face is classified on geom_region

  Downward quad_faces;
  Downward wedge_faces;

  int edge_counter=0, face_counter=0, rgn_counter=0;
  std::vector<MeshEntity*> btw_plane_edges;
  std::vector<MeshEntity*> btw_plane_faces;
  std::vector<MeshEntity*> btw_plane_regions;
  std::vector<MeshEntity*> vertex_edges;
  std::vector<MeshEntity*> edge_faces;

  Downward edgesNextPlane;
  Downward edgesBtwPlane;
  int num_upward;

  it = mesh->begin(2);
  while ((e = mesh->iterate(it)))
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
  }
  mesh->end(it);

  // exchange remote copies to set remote copy links
  PCU_Comm_Begin();
  int num_entities = num_local_vtx+num_local_edge+num_local_face;
  MeshEntity** entities = new MeshEntity*[num_entities];
  int pos=0;
  for (int index=0; index<num_local_vtx; ++index)
  {
    entities[pos]=remote_vertices[index];
    ++pos;
  }
  for (int index=0; index<num_local_edge; ++index)
  {
    entities[pos]=remote_edges[index];
    ++pos;
  }

  for (int index=0; index<num_local_face; ++index)
  {
    entities[pos]=remote_faces[index];
    ++pos;
  }
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &num_entities, sizeof(int));
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &(entities[0]),num_entities*sizeof(MeshEntity*));
  PCU_Comm_Send();
  delete [] entities;

  // receive phase begins
  std::vector<MeshEntity*> ent_vec;
  std::map<MeshEntity*, MeshEntity*> partbdry_ent_map;
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      int num_ent;
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      MeshEntity** s_ent = new MeshEntity*[num_ent];
      PCU_Comm_Unpack(&(s_ent[0]), num_ent*sizeof(MeshEntity*));
      pos=0;
      for (int index=0; index<num_local_vtx; ++index)
      {
        e = getMdsEntity(mesh, 0, index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }
      for (int index=0; index<num_local_edge; ++index)
      {
        e = getMdsEntity(mesh, 1, index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }
      for (int index=0; index<num_local_face; ++index)
      {
        e = getMdsEntity(mesh, 2, index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        ent_vec.push_back(e);
        ++pos;
      }
      delete [] s_ent;
    }
  }

  push_new_entities(mesh, partbdry_ent_map);
  bounce_orig_entities(mesh, ent_vec, m3dc1_model::instance()->prev_plane_partid,remote_vertices,remote_edges,remote_faces);

  // update partition classification
  for (int index=0; index<num_local_vtx; ++index)
  {
    e = getMdsEntity(mesh, 0, index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_vtx; ++index)
  {
    e = remote_vertices[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_edge; ++index)
  {
    e = getMdsEntity(mesh, 1, index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first); 
    parts.insert(local_partid);  
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_edge; ++index)
  {
    e = remote_edges[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_face; ++index)
  {
    e = getMdsEntity(mesh, 2, index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_face; ++index)
  {
    e = remote_faces[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  // update partition model
  mesh->acceptChanges();

  // update m3dc1_mesh internal data (# local, owned and global entities)
  update_partbdry(remote_vertices, remote_edges, remote_faces, 
                  btw_plane_edges, btw_plane_faces, btw_plane_regions);

  // delete existing local numbering
  apf::Numbering* local_n = mesh->findNumbering(mesh->getShape()->getName());
  if (local_n) destroyNumbering(local_n);
  
  // re-create the field and copy field data on master process group to non-master
  std::map<FieldID, m3dc1_field*>::iterator fit=m3dc1_mesh::instance()->field_container->begin();
  while(fit!=m3dc1_mesh::instance()->field_container->end())
  {
    update_field (fit->second->get_id(), fit->second->get_dof_per_value(), 
                  num_local_vtx, remote_vertices);
  }

  // clear temp memory
  delete [] remote_vertices;
  delete [] remote_edges;
  delete [] remote_faces;
  set_node_adj_tag();

#ifdef DEBUG
  printStats(mesh);
#endif

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Wedges between 2D planes created ***\n";
}

void compute_size_and_frame_fields(apf::Mesh2* m, double* size_1, double* size_2, 
     double* angle, apf::Field* sizefield, apf::Field* framefield)
{
  MeshEntity* e;
  MeshIterator* ent_it = m->begin(0);
  int i=0;
  while ((e = m->iterate(ent_it)))
  {
    double h1 = size_1[i];
    double h2 = size_2[i];

    double angle_1[3];
    angle_1[0] = angle[(i*3)];
    angle_1[1] = angle[(i*3)+1];
    angle_1[2] = angle[(i*3)+2];

 // Calculate the second unit vector
/*  double a, b;
    double frac_1, frac_2;
    frac_1 = (angle_1[0])*(angle_1[0]);
    frac_2 = (angle_1[0])*(angle_1[0]) + (angle_1[1])*(angle_1[1]);


    b = sqrt (frac_1/frac_2);
    a = -(angle_1[1]*b)/angle_1[0];

    double mag = sqrt (a*a + b*b);
    double dir_2[3];
    dir_2[0] = a /mag;
    dir_2[1] = b /mag;
    dir_2[2] = 0.0;
*/
    double dir_2[3];
    dir_2[0] = -angle_1[1];
    dir_2[1] =  angle_1[0];
    dir_2[2] =  0.0;

    ma::Vector h(h1, h2, h2);

    ma::Matrix r;
    r[0][0]=angle_1[0];
    r[0][1]=angle_1[1];
    r[0][2]=0.0;

    r[1][0]= dir_2[0];
    r[1][1]= dir_2[1];
    r[1][2]=0.0;

    r[2][0]=0;
    r[2][1]=0;
    r[2][2]=1.;

    apf::setVector(sizefield, e, 0, h);
    apf::setMatrix(framefield, e, 0, r);
//  apf::setMatrix(framefield, e, 0, apf::transpose(r));	// For Shock Test Case
  }
  m->end(ent_it);

  // sync the fields to make sure verts on part boundaries end up with the same size and frame
  apf::synchronize(sizefield);
  apf::synchronize(framefield);
}

// *********************************************************
void copy_field (apf::Mesh2* mesh, apf::Field* f, double* master_data_array)
// *********************************************************
{
  int num_2d_vtx=mesh->count(0);

  // copy the existing dof data
  int total_ndof = apf::countComponents(f);

  double* dof_val = new double[total_ndof*num_2d_vtx];

  freeze(f); // switch dof data from tag to array

  PCU_Comm_Begin();

  int proc=PCU_Comm_Self()%m3dc1_model::instance()->group_size;
  if (m3dc1_model::instance()->local_planeid)
  {
    memcpy(&(dof_val[0]), apf::getArrayData(f), total_ndof*num_2d_vtx*sizeof(double));
    PCU_Comm_Pack(proc, &num_2d_vtx, sizeof(int));
    PCU_Comm_Pack(proc, &(dof_val[0]), total_ndof*num_2d_vtx*sizeof(double));
  }
  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_val;

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      int index=total_ndof*recv_num_ent*from;
      PCU_Comm_Unpack(&(master_data_array[index]), total_ndof*recv_num_ent*sizeof(double));
    }
  }
  delete [] dof_val;
}


// *********************************************************
void copy_field (apf::Mesh2* mesh, apf::Field* f, 
                 map<int, apf::Field*>& new_fields, int field_id)
// *********************************************************
{
  int num_2d_vtx=mesh->count(0);
  int num_rank = pumi_size();
  // copy the existing dof data
  int total_ndof = apf::countComponents(f);

  double* dof_data = new double[total_ndof*num_2d_vtx];

  freeze(f); // switch dof data from tag to array

  PCU_Comm_Begin();

  int proc=PCU_Comm_Self()%m3dc1_model::instance()->group_size;
  if (m3dc1_model::instance()->local_planeid)
  {
    memcpy(&(dof_data[0]), apf::getArrayData(f), total_ndof*num_2d_vtx*sizeof(double));
    PCU_Comm_Pack(proc, &num_2d_vtx, sizeof(int));
    PCU_Comm_Pack(proc, &(dof_data[0]), total_ndof*num_2d_vtx*sizeof(double));
  }
  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_data = new double[total_ndof*num_2d_vtx];

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      if (!new_fields[field_id*num_rank+from])
        std::cout<<"p"<<PCU_Comm_Self()<<"] cannot find new_fields["<<field_id*num_rank+from<<"\n";
      dof_data = apf::getArrayData(new_fields[field_id*num_rank+from]);
      PCU_Comm_Unpack(&(dof_data[0]), total_ndof*recv_num_ent*sizeof(double));
#ifdef DEBUG
      std::cout<<"[p"<<PCU_Comm_Self()<<"] saving field "
               <<field_id<<"@p"<<from<<" into new_fields["<<field_id*num_rank+from<<"] \n";
#endif
    }
  }
}

// run mesh adaptation in multiple planes in 3D
void adapt_mesh (int field_id_h1, int field_id_h2, double* dir,
    int shouldSnap, int shouldRunPreZoltan ,int shouldRunPostZoltan,
    int shouldRefineLayer, int maximumIterations, double goodQuality)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  apf::Field* f_h1 = (*m3dc1_mesh::instance()->field_container)[field_id_h1]->get_field();
  synchronize_field(f_h1);
  int num_dof = countComponents(f_h1);
  if (!isFrozen(f_h1)) freeze(f_h1);
  double* data_h1;

  if (num_dof==1)
    data_h1 = apf::getArrayData(f_h1);
  else 
  {
    data_h1 = new double[mesh->count(0)];
    double* temp_data = apf::getArrayData(f_h1);
    for (int i=0; i<mesh->count(0); ++i)
      data_h1[i] = temp_data[i*num_dof];
  }

  apf::Field* f_h2 = (*m3dc1_mesh::instance()->field_container)[field_id_h2]->get_field();
  synchronize_field(f_h2);
  if (!isFrozen(f_h2)) freeze(f_h2);

  double* data_h2;
  if (num_dof==1)
    data_h2 = apf::getArrayData(f_h2);
  else
  {
    data_h2= new double[mesh->count(0)];
    double* temp_data = apf::getArrayData(f_h2);
    for (int i=0; i<mesh->count(0); ++i)
      data_h2[i] = temp_data[i*num_dof];
  }

  apf::Field* size_field = apf::createField(mesh, "size_field", apf::VECTOR, apf::getLagrange(1));
  apf::Field* frame_field = apf::createField(mesh, "frame_field", apf::MATRIX, apf::getLagrange(1));

  compute_size_and_frame_fields(mesh, data_h1, data_h2, dir, size_field, frame_field);
  	 
   m3dc1_field_delete (&field_id_h1);
   m3dc1_field_delete (&field_id_h2);

  if (num_dof>1)
  {
    delete [] data_h1;
    delete [] data_h2;
  }

  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*>::iterator it=m3dc1_mesh::instance()->field_container->begin();

  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 1);
    fields.push_back(field);
    it++;
  }
  
  apf::unfreezeFields(mesh); // turning field data from array to tag

  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": "<<fields.size() <<" fields will be transfered\n";

  // delete all the matrix
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    apf::destroyNumbering(n);
  }
	
  ReducedQuinticImplicit shape;
  ReducedQuinticTransfer slnTransfer(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh, size_field, frame_field, &slnTransfer);
	
  in->shouldSnap = 0; // FIXME: crash if *shouldSnap==1;
  in->shouldTransferParametric = 0;
  in->shouldRunPreZoltan = shouldRunPreZoltan;
  in->shouldRunPostZoltan = shouldRunPostZoltan;
  in->shouldRunMidParma = 0;
  in->shouldRunPostParma = 0;
  in->shouldRefineLayer = shouldRefineLayer;
  in->maximumIterations=maximumIterations;
  in->goodQuality = goodQuality;

  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": snap "<<shouldSnap
  	  <<", runPreZoltan "<<shouldRunPreZoltan<<", runPostZoltan "<<shouldRunPostZoltan<<"\n";

#ifdef DEBUG
  apf::writeVtkFiles("before-adapt", mesh);
  mesh->writeNative("mesh.smb");
#endif

  apf::MeshEntity* e;

  gmi_model* model = m3dc1_model::instance()->model;
  int** ge_tag = new int*[3];
  gmi_iter* g_it;
  gmi_ent* ge;
  int cnt;

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    // save the current model tag and revert it back to the original
    for (int i = 0; i < 3; ++i)
      ge_tag[i] = new int[model->n[i]];

    for (int dim=0; dim<=2; ++dim)
    {
      g_it = gmi_begin(model, dim);
      cnt=0;
      while ((ge = gmi_next(model, g_it)))
      {
        ge_tag[dim][cnt++]= gmi_tag(model,ge);
        gmi_base_set_tag (model, ge, cnt);
      }
      gmi_end(model, g_it);
    }

    // remove 3D entities
    m3dc1_mesh::instance()->remove_wedges();
    apf::printStats(mesh);

    apf::writeVtkFiles("after-wedge-removal", mesh);
    mesh->writeNative("after-wedge-removal.smb");
  }

  ma::adapt(in);

  mesh->removeField(size_field);
  mesh->removeField(frame_field);
  apf::destroyField(size_field);
  apf::destroyField(frame_field);

  if (!PCU_Comm_Self()) std::cout<<"After adaptation: ";
  m3dc1_mesh::instance()->print(__LINE__);

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    m3dc1_mesh::instance()->set_mcount();
    // re-construct wedges
    m3dc1_mesh::instance()->create_wedges();

    for(int i = 0; i < 3; ++i)
      delete [] ge_tag[i];
    delete [] ge_tag;
  } // 3d 

  // FIXME: crash in 2D if no pre/post zoltan 
  reorderMdsMesh(mesh);

  // FIXME: crash in 3D 
  apf::writeVtkFiles("after-adapt", mesh);

  m3dc1_mesh::instance()->initialize();

  compute_globalid(mesh, 0);
  compute_globalid(mesh, mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    synchronize_field(field);
    it++;
  }

  apf::freezeFields(mesh); // turning field data from tag to array
}
