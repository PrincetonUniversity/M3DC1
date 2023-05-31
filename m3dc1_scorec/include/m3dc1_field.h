/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_FIELD_H
#define M3DC1_FIELD_H

#include "apf.h"

void group_complex_dof (apf::Field* field, int option);

class m3dc1_field
{
public:
  m3dc1_field (int i, apf::Field* f, int n, int t, int ndof): id(i), field(f), num_value(n), value_type(t),dof_per_value(ndof) {transfer = false;}
  ~m3dc1_field() {}
  apf::Field* get_field() { return field; }
  void set_field(apf::Field* f) { field = f; }
  int get_id() { return id; }
  int get_num_value() { return num_value; }
  int get_value_type() { return value_type; }
  int get_dof_per_value() {return dof_per_value;}
  bool should_transfer() {return transfer;}
  void mark_for_solutiontransfer() {transfer = true;}
private:
  int id;
  apf::Field* field; // name and #dofs are available from apf::Field
  int num_value;
  int value_type;
  int dof_per_value;
  bool transfer;
};

void synchronize_field(apf::Field* f);
#endif
