/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_ghost.h"
#include "m3dc1_mesh.h"
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "m3dc1_field.h"
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "PCU.h"
#include <vector>
#include <set>
#include <string.h>
#include <iostream>

using namespace apf;

// m3dc1_ghost
// *********************************************************
m3dc1_ghost::m3dc1_ghost()
// *********************************************************
{
  mesh = NULL;
  nlayers = 0;
  field_container=NULL;
  reset();
}

// *********************************************************
m3dc1_ghost::~m3dc1_ghost()
// *********************************************************
{
  if (!mesh)
    return;
  clean();
  mesh->destroyNative();
  apf::destroyMesh(mesh);
}

void m3dc1_ghost::clean()
{
 // destroy field
  if (field_container)
  {
    for (std::map<FieldID, m3dc1_field*>::iterator f_it=field_container->begin(); f_it!=field_container->end();)
    {
      FieldID id = f_it->first;
      std::map<FieldID, m3dc1_field*>::iterator it_next=++f_it;
      m3dc1_field_delete(&id);
      f_it=it_next;
    }
  }
  delete field_container; 
  field_container=0;
}

m3dc1_ghost* m3dc1_ghost::_instance=NULL;

m3dc1_ghost* m3dc1_ghost::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_ghost();
  return _instance;
}

void m3dc1_ghost::destroy()
{
  delete _instance;
  _instance = NULL;
}

// *********************************************************
void m3dc1_ghost::reset()
// *********************************************************
{
  for (int i=0; i<4; ++i)
  {
    num_local_ent[i] = 0;
    num_own_ent[i] = 0;
    num_global_ent[i] = 0;
  }
}

// *********************************************************
void m3dc1_ghost::initialize()
// *********************************************************
{
  reset();
  for (int d=0; d<4; ++d)
  {
    num_local_ent[d] = countEntitiesOfType(mesh, d);
    num_own_ent[d] = m3dc1_mesh::instance()->num_own_ent[d];
    num_global_ent[d] = m3dc1_mesh::instance()->num_global_ent[d];
  }

  // Copy the fields on the ghosted mesh into the field container
  if (!m3dc1_ghost::instance()->field_container)
    m3dc1_ghost::instance()->field_container = new std::map<FieldID, m3dc1_field*>;
  for (std::map<FieldID, m3dc1_field*>::iterator it =
	 m3dc1_mesh::instance()->field_container->begin();
       it !=  m3dc1_mesh::instance()->field_container->end();
       ++it) 
  {
    int field_id = it->first;
    int num_values = it->second->get_num_value();
    int scalar_type = it->second->get_value_type();
    int num_dofs_per_value = it->second->get_dof_per_value();
    apf::Field* old_field = it->second->get_field();
    apf::Field* new_field = mesh->findField(apf::getName(old_field));
    
    m3dc1_ghost::instance()->field_container->insert(
        std::map<FieldID, m3dc1_field*>::value_type(field_id,
          new m3dc1_field(field_id,
                          new_field,
                          num_values,
                          scalar_type,
                          num_dofs_per_value)));
    apf::freeze(new_field);
  }
}

