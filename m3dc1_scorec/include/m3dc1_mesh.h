/******************************************************************************
  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef M3DC1_MESH_H
#define M3DC1_MESH_H
#include "apfMDS.h"
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
// use pumi_ment_isOwned with m3dc1_mesh::instance()->ownership
bool is_ent_original(apf::MeshEntity* e);

// use pumi_ment_getOwnPID with m3dc1_mesh::instance()->ownership
int get_ent_ownpartid(apf::MeshEntity* ent);

// use pumi_ment_getOwnEnt with m3dc1_mesh::instance()->ownership
apf::MeshEntity* get_ent_owncopy(apf::MeshEntity* ent);

int get_ent_localid (apf::Mesh2* mesh, apf::MeshEntity* ent);
int get_ent_globalid (apf::Mesh2* mesh, apf::MeshEntity* ent);
void verify_field(pMesh m, pField f);

class m3dc1_mesh
{
private:
  void set_node_adj_tag();
  apf::Mesh2 * mesh;
  // field container
  std::map<FieldID, m3dc1_field*> field_container;
  static m3dc1_mesh* _instance;
  // local counter to avoid apf/pumi acually counting the ents every time
  //  need to update after adapted
  int num_local_ent[4];
  int num_global_ent[4];
  int num_own_ent[4];
  // tag for local entity id
  apf::MeshTag * local_entid_tag;
  // tags for second order adjanceny info
  apf::MeshTag * num_global_adj_node_tag;
  apf::MeshTag * num_own_adj_node_tag;
public:
  // user-defined ownership
  pOwnership ownership;
  apf::MeshTag* own_partid_tag;

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

  int get_local_count(int dim)
  {
    return num_local_ent[dim];
  }

  int get_global_count(int dim)
  {
    return num_global_ent[dim];
  }

  int get_own_count(int dim)
  {
    return num_own_ent[dim];
  }

  // TODO : get rid of these setter functions, let the mesh know when the
  //   underlying mesh has been modified and update the cached entity counts
  void update_local_count(int dim)
  {
    num_local_ent[dim] = mesh->count(dim);
  }

  // collective on the current PCU_COMM_WORLD
  void update_global_count(int dim)
  {
    update_local_count(dim);
    MPI_Allreduce(&num_own_ent[dim],&num_global_ent[dim],1,MPI_INT,MPI_SUM,PCU_Get_Comm());
  }

  void update_own_count(int dim)
  {
    num_own_ent[dim] = apf::countOwned(mesh,dim);
  }

  apf::MeshTag * local_id_tag() { return local_entid_tag; }
  // tag for owned partid attached to the part bdry entities
  apf::MeshTag * own_bridge_adj_tag() { return num_own_adj_node_tag; }
  apf::MeshTag * global_bridge_adj_tag() { return num_global_adj_node_tag; }
};
template <class O>
void m3dc1_mesh::retrieve_fields(O out)
{
  for(auto fld_it = field_container.begin(); fld_it != field_container.end(); ++fld_it)
    *out++ = fld_it->second;
}
#endif
