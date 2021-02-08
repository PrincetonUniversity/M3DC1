/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_MESH_H
#define M3DC1_MESH_H
#include "map"
#include <set>
#include "utility"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apf.h"
#include "m3dc1_scorec.h"
#include "m3dc1_field.h"

void compute_globalid(apf::Mesh2* m, int d);

bool is_ent_original(apf::Mesh2* mesh, apf::MeshEntity* e);
int get_ent_ownpartid(apf::Mesh2* mesh, apf::MeshEntity* ent);
apf::MeshEntity* get_ent_owncopy(apf::Mesh2* mesh, apf::MeshEntity* ent);
int get_ent_localid (apf::Mesh2* mesh, apf::MeshEntity* ent);
int get_ent_globalid (apf::Mesh2* mesh, apf::MeshEntity* ent);

int get_ent_global2ndadj (apf::Mesh2*, int, int, std::vector<apf::MeshEntity*>&,
    std::vector<int>&, std::vector<int>&, std::vector<int>&);
void get_ent_numglobaladj(apf::Mesh2*, int, int, std::vector<apf::MeshEntity*>&, std::vector<int>&);

class m3dc1_mesh
{
public:
  m3dc1_mesh();
  ~m3dc1_mesh();
  static m3dc1_mesh* instance();
  // functions
  void reset();
  void clean();
  void remove_wedges();
  void build3d(int num_field, int* field_id, int* num_dofs_per_value); // old: setup3DMesh(pPart mesh, pGeomMdl model,int ifXYZ) in PlaneManager.h
  void initialize(); // to be called after initial mesh loading. old: updatemeshinfo_
  void update_partbdry(apf::MeshEntity** remote_vertices, apf::MeshEntity** remote_edges, 
              apf::MeshEntity** remote_faces, std::vector<apf::MeshEntity*>& btw_plane_edges, 
              std::vector<apf::MeshEntity*>& btw_plane_faces, std::vector<apf::MeshEntity*>& btw_plane_regions);

  void print(int);

  // data
  apf::Mesh2* mesh;
  apf::MeshEntity*** ments;

  // local counter for fast info retrieval
  int num_local_ent[4];
  int num_global_ent[4];
  int num_own_ent[4];

  // field container 
  std::map<FieldID, m3dc1_field*>* field_container;

  // tag for local entity id
  apf::MeshTag* local_entid_tag;

  // tag for owned partid attached to the part bdry entities
  apf::MeshTag* own_partid_tag; 

  // tags for second order adjanceny info
  apf::MeshTag* num_global_adj_node_tag;
  apf::MeshTag* num_own_adj_node_tag;
private:
  void set_node_adj_tag();
  static m3dc1_mesh* _instance;
};
#endif
