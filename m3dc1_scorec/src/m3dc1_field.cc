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
#include "m3dc1_numbering.h"

m3dc1_field::m3dc1_field (int ID, const char* str, int nv, int t, int ndof)
  : id(ID)
  , name(str)
  , fld(NULL)
  , num_value(nv)
  , value_type(t)
  , num_dof(ndof)
  , num(NULL)
{
  int n_components = nv * ndof * (t+1);
  fld = apf::createPackedField(m3dc1_mesh::instance()->mesh, str, n_components);
  num = apf::createNumbering(fld);
  apf::freeze(fld);
  aggregateNumbering(MPI_COMM_SELF,num,nv,ndof);
  // freeze the numbering
  /*
  char * field_name = new char[32];
  for (int i=0; i<nv; ++i)
  {
    n_components = (t+1)*num_dof;
    sprintf(field_name,"%s_%d",str,i);
    fields[i]=apf::createPackedField(m3dc1_mesh::instance()->mesh, field_name, n_components);
    apf::freeze(fields[i]); // switch dof data from tag to array
  }
  delete [] field_name;
  */
}

m3dc1_field::~m3dc1_field()
{
  /*
  for (int i=0; i<num_value; ++i)
    destroyField(fields[i]);
  delete [] fields;
  */
  destroyField(fld);
  fld = NULL;
}

void get_ent_localdofid(m3dc1_field * mf, int ent_lid, int* dof_id, int* dof_cnt)
{
  int blks_per_nd = mf->get_num_value();
  int dofs_per_blk = mf->get_dof_per_value();
  int dofs_per_nd = dofs_per_blk * blks_per_nd;
  int num_local_node = m3dc1_mesh::instance()->mesh->count(0);
  int ii = 0;
  for (int blk = 0; blk < blks_per_nd; ++blk)
    for (int dof = 0; dof < dofs_per_blk; ++dof)
      dof_id[ii++] = ((ent_lid + (blk * num_local_node)) * dofs_per_blk) + dof;
  *dof_cnt=dofs_per_nd;
}

// this assumes the entity is a vert
void get_ent_globaldofid(m3dc1_field * mf, int ent_lid, int* dof_id, int* dof_cnt)
{
  //assumes only verts hold nodes/dofs
  apf::Mesh2 * msh = static_cast<apf::Mesh2*>(apf::getMesh(mf->get_field()));
  apf::MeshEntity * vrt = apf::getMdsEntity(msh,0,ent_lid);
  apf::Field * fld = mf->get_field();
  apf::Numbering * num = mf->get_global_numbering();
  *dof_cnt = apf::countComponents(fld);
  for(int cmp = 0; cmp < *dof_cnt; ++cmp)
    dof_id[cmp] = apf::getNumber(num,vrt,0,cmp);
}


//*******************************************************
void get_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
//*******************************************************
{
  apf::Field * f = mf->get_field();
  getComponents(f, e, 0, &(dof_data[0]));
}

//*******************************************************
void set_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
//*******************************************************
{
  apf::Field* f = mf->get_field();
  setComponents(f, e, 0, &(dof_data[0]));
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

