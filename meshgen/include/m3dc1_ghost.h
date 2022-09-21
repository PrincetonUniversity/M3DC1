/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_GHOST_H
#define M3DC1_GHOST_H
#include "map"
#include <set>
#include "utility"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apf.h"
#include "m3dc1_scorec.h"
#include "m3dc1_field.h"

class m3dc1_ghost
{
public:
  m3dc1_ghost();
  ~m3dc1_ghost();
  static m3dc1_ghost* instance();
  static void destroy();
  // functions
  void reset();
  void clean();
  void initialize(); // to be called after creating ghost mesh.

  // data
  apf::Mesh2* mesh;
  int nlayers;
  int ordering_opt;

  // local counter for fast info retrieval
  int num_local_ent[4];
  int num_global_ent[4];
  int num_own_ent[4];

  // field container 
  std::map<FieldID, m3dc1_field*>* field_container;

  std::vector<bool>* org_node_flag;

private:
  static m3dc1_ghost* _instance;
};
#endif


