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
#include <sstream>
#include <iostream>
#include <PCU.h>
#include <assert.h>
#include "apfMDS.h"

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

void get_ent_localdofid(m3dc1_field* mf, int ent_lid, int* dof_id, int* dof_cnt)
{
  int num_value = mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int num_dof = num_value*dof_per_value;
  int num_local_node = m3dc1_mesh::instance()->mesh->count(0);

  int i=0;

  for (int nv=0; nv<num_value; ++nv)
    for (int nd=0; nd<dof_per_value; ++nd)
      dof_id[i++] = nv*num_local_node*dof_per_value+ent_lid*dof_per_value+nd;
  *dof_cnt=num_dof;
}

void get_ent_globaldofid(m3dc1_field* mf, int ent_gid, int* dof_id, int* dof_cnt)
{
  int num_value = mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int num_dof = num_value*dof_per_value;
  int num_global_node = m3dc1_mesh::instance()->num_global_ent[0];

  int i=0;

  for (int nv=0; nv<num_value; ++nv)
    for (int nd=0; nd<dof_per_value; ++nd)
      dof_id[i++] = nv*num_global_node*dof_per_value+ent_gid*dof_per_value+nd;
  *dof_cnt=num_dof;
}


//*******************************************************
void get_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
//*******************************************************
{
  int nv = mf->get_num_value();
  apf::Field* f = mf->get_field(0);
  int ndof=countComponents(f);

  for (int i=0; i<nv; ++i)
  {
    f = mf->get_field(i);
    getComponents(f, e, 0, &(dof_data[ndof*i]));
  }
}

//*******************************************************
void set_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
//*******************************************************
{
  int nv = mf->get_num_value();
  apf::Field* f = mf->get_field(0);
  int ndof=countComponents(f);

  for (int i=0; i<nv; ++i)
  {
    f = mf->get_field(i);
    setComponents(f, e, 0, &(dof_data[ndof*i]));
  }
}


//=========================================================================
void load_field(apf::Mesh2* m, int field_id, const char* filename)
{
  std::string in(filename);
  std::stringstream s;
  s <<in<< "-"<<PCU_Comm_Self();
  std::string partFile = s.str();
  FILE * fp =fopen(partFile.c_str(), "r");

  if (!fp) 
    std::cout<<"("<<PCU_Comm_Self()<<") [M3D-C1 ERROR] fail to load file \""<<partFile<<"\"\n";

  apf::MeshEntity* e;

  int gid, lid, did, ndof, nv, nd, vt, start_index; 
  double dof;
  char field_name[32];
  fscanf(fp, "%s %d %d %d %d\n", field_name, &nv, &vt, &nd, &start_index);
  std::cout<< field_name <<" "<<nv<<" "<<vt<<" "<<nd<<" "<<start_index<<"\n";

  if (m3dc1_mesh::instance()->field_container)
    assert(m3dc1_mesh::instance()->field_container->count(field_id)==0);

  m3dc1_field_create (&field_id, field_name, &nv, &vt, &nd);
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[field_id];
  assert(mf->get_num_value()==nv && mf->get_dof_per_value()==nd && mf->get_value_type()==vt);

  int num_dof=mf->get_num_value()*mf->get_dof_per_value();
#ifdef PETSC_USE_COMPLEX
  num_dof*=2;
#endif
  double* dof_data = new double[num_dof];
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    fscanf(fp, "%d %d %d\n", &gid, &lid, &ndof);
    assert(get_ent_globalid(m,e)+start_index==gid && getMdsIndex(m, e)==lid);
    for (int i=0; i<num_dof; ++i)
      dof_data[i]=0.0;
    for (int i=0; i<ndof; ++i)
    {
      fscanf(fp, "%d %lf",&did, &dof);
      dof_data[did]=dof;
    }
    set_ent_dofdata(mf, e, dof_data);
  } // while
  m->end(it);
  fclose(fp);  
  delete [] dof_data;
}


//=========================================================================
void write_field(apf::Mesh2* m, m3dc1_field* mf, const char* filename, int start_index)
{
  std::string in(filename);
  std::stringstream s;
  s << in<< "-"<<PCU_Comm_Self();
  std::string partFile = s.str();
  FILE * fp =fopen(partFile.c_str(), "w");
  assert(fp);

  apf::MeshEntity* e;

  int num_dof=mf->get_num_value()*mf->get_dof_per_value();
#ifdef PETSC_USE_COMPLEX
  num_dof*=2;
#endif
  double* dof_data = new double[num_dof];

  fprintf(fp, "%s %d %d %d %d\n", mf->get_name().c_str(), mf->get_num_value(), 
          mf->get_value_type(), mf->get_dof_per_value(), start_index);
 
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    get_ent_dofdata(mf, e, dof_data);
    int ndof=0;
    for (int i=0; i<num_dof; ++i)
    {     
      if (!m3dc1_double_isequal(dof_data[i], 0.0))
        ++ndof;
    }
    
    fprintf(fp, "%d %d %d\n", get_ent_globalid(m,e)+start_index, getMdsIndex(m, e), ndof);   
    for (int i=0; i<num_dof; ++i)
    {     
      if (!m3dc1_double_isequal(dof_data[i], 0.0))
        fprintf(fp, "%d %lf\n", i, dof_data[i]);
    }
  } // while
  m->end(it);
  fclose(fp);  
  delete [] dof_data;
}

