/******************************************************************************
  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef M3DC1_MESH_H
#define M3DC1_MESH_H
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apf.h"
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "m3dc1_field.h"
#include <PCU.h>
#include <pumi.h>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
void compute_globalid(apf::Mesh2* m, int d);
bool is_ent_original(apf::Mesh2* mesh, apf::MeshEntity* e);
// instead, call pumi_ment_getOwnPID(pMeshEnt e)
//int get_ent_ownpartid(apf::Mesh2* mesh, apf::MeshEntity* ent);
apf::MeshEntity* get_ent_owncopy(apf::Mesh2* mesh, apf::MeshEntity* ent);
int get_ent_localid (apf::Mesh2* mesh, apf::MeshEntity* ent);
int get_ent_globalid (apf::Mesh2* mesh, apf::MeshEntity* ent);
void verify_field(pMesh m, pField f);
class m3dc1_mesh
{
public:
  m3dc1_mesh();
  ~m3dc1_mesh();
  static m3dc1_mesh* instance();
  // functions
  void reset();
  void clean();
  void build3d(int num_field, int* field_id, int* num_dofs_per_value); // old: setup3DMesh(pPart mesh, pGeomMdl model,int ifXYZ) in PlaneManager.h
  void initialize(); // to be called after initial mesh loading. old: updatemeshinfo_
  void update_partbdry(apf::MeshEntity** remote_vertices, apf::MeshEntity** remote_edges,
              apf::MeshEntity** remote_faces, std::vector<apf::MeshEntity*>& btw_plane_edges,
              std::vector<apf::MeshEntity*>& btw_plane_faces, std::vector<apf::MeshEntity*>& btw_plane_regions);
  void print(int);

  void add_field(FieldID fid, m3dc1_field * fld)
  {
    auto loc_new = field_container.insert(std::make_pair(fid,fld));
    assert(loc_new.second && "[M3D-C1 Error] Field with specified ID already exists!");
  }
  m3dc1_field * get_field(FieldID fid)
  {
    auto fld_it = field_container.find(fid);
    assert(fld_it != field_container.end() && "[M3D-C1 Error] Field with specified ID does not exist!");
    return fld_it->second;
  }
  void destroy_field(FieldID fid)
  {
    auto fld_it = field_container.find(fid);
    assert(fld_it != field_container.end() && "[M3D-C1 Error] Field with specified id does not exist!");
    apf::destroyField(fld_it->second->get_field());
    field_container.erase(fld_it);
  }
  bool field_exists(FieldID fid)
  {
    return field_container.find(fid) != field_container.end();
  }
  void verify_fields()
  {
    for(auto fld = field_container.begin(); fld != field_container.end(); ++fld)
    {
      if(!PCU_Comm_Self())
        std::cout << "[M3D-C1 INFO] verifying field " << fld->second->get_name() << std::endl;
      verify_field(m3dc1_mesh::instance()->mesh, fld->second->get_field());
    }
  }
  template <class O>
  void retrieve_fields(O out);

  void load_mesh(const char * fn)
  {
    mesh = pumi_mesh_load(pumi::instance()->model,
                          fn,
                          m3dc1_model::instance()->group_size);
  }
  void create_mesh()
  {
    mesh = pumi_mesh_create(pumi::instance()->model,
                            2,
                            false);
  }
  apf::Mesh2 * get_mesh() { return mesh; }

  // data
  apf::MeshEntity*** ments;

  // local counter for fast info retrieval
  int num_local_ent[4];
  int num_global_ent[4];
  int num_own_ent[4];

  // tag for local entity id
  apf::MeshTag* local_entid_tag;

  // tag for owned partid attached to the part bdry entities
  //  apf::MeshTag* own_partid_tag;

  // tags for second order adjanceny info
  apf::MeshTag* num_global_adj_node_tag;
  apf::MeshTag* num_own_adj_node_tag;
  //void set_node_adj_tag2();
private:
  void set_node_adj_tag();
  apf::Mesh2 * mesh;
  // field container
  std::map<FieldID, m3dc1_field*> field_container;
  static m3dc1_mesh* _instance;
};
template <class O>
void m3dc1_mesh::retrieve_fields(O out)
{
  for(auto fld_it = field_container.begin(); fld_it != field_container.end(); ++fld_it)
    *out++ = fld_it->second;
}
#endif
