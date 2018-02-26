/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_FIELD_H
#define M3DC1_FIELD_H

#include "apf.h"
#include "apfMesh2.h"

class m3dc1_field
{
public:
  m3dc1_field (int i, const char* str, int n, int t, int ndof);
  ~m3dc1_field();

  apf::Field* get_field(int vid);

  int get_id() { return id; }
  std::string get_name() { return name; }
  int get_num_value() { return num_value; }
  int get_value_type() { return value_type; }
  int get_dof_per_value() {return num_dof;}
private:
  int id;
  std::string name;
  apf::Field** fields; // name and #dofs are available from apf::Field
  int num_value;
  int value_type;
  int num_dof;
};

void get_ent_localdofid(m3dc1_field* mf, int ent_lid, int* dof_id, int* dof_cnt);
void get_ent_globaldofid(m3dc1_field* mf, int ent_gid, int* dof_id, int* dof_cnt);
void get_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data);
void set_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data);

// given a field ID, create a field out of file
void load_field(apf::Mesh2* m, int field_id, const char* filename);
void write_field(apf::Mesh2* m, m3dc1_field* mf, const char* filename, int start_index);
#endif
