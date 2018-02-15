/****************************************************************************** 

  (c) 2005-2018 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include <stdio.h>

m3dc1_field::m3dc1_field (int ID, const char* str, int nv, int t, int ndof)
: id(ID), name(str), num_value(nv), value_type(t), num_dof(ndof)
{
  name=str;
  typedef apf::Field* p_field;
  fields = new p_field[nv];
  int n_components;
  char* field_name=new char[32];
  for (int i=0; i<nv; ++i)
  {
    n_components = (t+1)*num_dof;
    sprintf(field_name,"%s_%d",str,i);
    fields[i]=apf::createPackedField(m3dc1_mesh::instance()->mesh, field_name, n_components);
    apf::freeze(fields[i]); // switch dof data from tag to array
  }
}

m3dc1_field::~m3dc1_field()
{
  for (int i=0; i<num_value; ++i)
    destroyField(fields[i]);

  delete [] fields;
}

apf::Field* m3dc1_field::get_field(int vid)
{ 
  return fields[vid]; 
}

