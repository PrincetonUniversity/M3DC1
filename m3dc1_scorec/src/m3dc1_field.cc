/******************************************************************************

  (c) 2005-2018 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include "m3dc1_numbering.h"
#include <apf.h>
#include <apfField.h>
#include <apfMDS.h>
#include <PCU.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
MPI_Comm getAggregationComm(int agg_scp)
{
  return agg_scp == m3dc1_field::LOCAL_AGGREGATION ? MPI_COMM_SELF :
    agg_scp == m3dc1_field::PLANE_AGGREGATION ? m3dc1_model::instance()->getPlaneComm() :
    agg_scp == m3dc1_field::GLOBAL_AGGREGATION ? M3DC1_COMM_WORLD : MPI_COMM_NULL;
}

m3dc1_field::m3dc1_field(int ID,
                         const char * str,
                         int num_blks,
                         int num_dofs,
                         aggregation_scope scp)
  : id(ID)
  , name(str)
  , fld(NULL)
  , gbl_num(NULL)
  , lcl_num(NULL)
  , agg_scp(scp)
  , blks_per_nd(num_blks)
  , dofs_per_blk(num_dofs)
  , dofs_per_nd(num_blks*num_dofs)
{
  int cmps = dofs_per_nd;
#ifdef PETSC_USE_COMPLEX
  cmps *= 2;
#endif
  apf::Mesh2 * msh = m3dc1_mesh::instance()->get_mesh();
  fld = apf::createPackedField(msh, str, cmps);
  std::string gnm(str);
  gnm += "gbl_num";
  gbl_num = apf::createNumbering(msh,gnm.c_str(),fld->getShape(),cmps);
  std::string lnm(str);
  lnm += "lcl_num";
  lcl_num = apf::createNumbering(msh,lnm.c_str(),fld->getShape(),cmps);
  apf::freeze(fld);
  apf::zeroField(fld);
  MPI_Comm cm = getAggregationComm(agg_scp);
  aggregateNumbering(cm,gbl_num,blks_per_nd,dofs_per_blk,M3DC1_COMM_WORLD);
  aggregateNumbering(MPI_COMM_SELF,lcl_num,blks_per_nd,dofs_per_blk,MPI_COMM_SELF);
}

m3dc1_field::~m3dc1_field()
{
  destroyField(fld);
  fld = NULL;
}

void get_ent_xxxdofid(bool lcl, m3dc1_field * mf, int ent_lid, int * dof_id, int * dof_cnt)
{
  if(lcl)
    get_ent_localdofid(mf,ent_lid,dof_id,dof_cnt);
  else
    get_ent_globaldofid(mf,ent_lid,dof_id,dof_cnt);
}

void get_ent_localdofid(m3dc1_field * mf, int ent_lid, int* dof_id, int* dof_cnt)
{
  apf::Mesh2 * msh = static_cast<apf::Mesh2*>(apf::getMesh(mf->getCoreField()));
  apf::MeshEntity * vrt = apf::getMdsEntity(msh,0,ent_lid);
  apf::Field * fld = mf->getCoreField();
  apf::Numbering * num = mf->getLocalNumbering();
  *dof_cnt = apf::countComponents(fld);
  for(int cmp = 0; cmp < *dof_cnt; ++cmp)
    dof_id[cmp] = apf::getNumber(num,vrt,0,cmp);
}

// this assumes the entity is a vert
void get_ent_globaldofid(m3dc1_field * mf, int ent_lid, int* dof_id, int* dof_cnt)
{
  //assumes only verts hold nodes/dofs
  apf::Mesh2 * msh = static_cast<apf::Mesh2*>(apf::getMesh(mf->getCoreField()));
  apf::MeshEntity * vrt = apf::getMdsEntity(msh,0,ent_lid);
  apf::Field * fld = mf->getCoreField();
  apf::Numbering * num = mf->getGlobalNumbering();
  *dof_cnt = apf::countComponents(fld);
  for(int cmp = 0; cmp < *dof_cnt; ++cmp)
    dof_id[cmp] = apf::getNumber(num,vrt,0,cmp);
}

void get_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
{
  apf::Field * f = mf->getCoreField();
  getComponents(f, e, 0, &(dof_data[0]));
}

void set_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data)
{
  apf::Field* f = mf->getCoreField();
  setComponents(f, e, 0, &(dof_data[0]));
}

// TODO : why pass in an apf mesh if we're going to use the m3dc1_mesh anyway?
void load_field(apf::Mesh2 * m, int fid, const char* filename)
{
  std::string in(filename);
  std::stringstream s;
  s <<in<< "-"<<PCU_Comm_Self();
  std::string partFile = s.str();
  FILE * fp =fopen(partFile.c_str(), "r");

  if (!fp)
    std::cout<<"("<<PCU_Comm_Self()<<") [M3D-C1 ERROR] fail to load file \""<<partFile<<"\"\n";

  int gid, lid, did, ndof, nv, nd, start_index, agg_scp;
  double dof;
  char field_name[32];
  fscanf(fp, "%s %d %d %d %d\n", field_name, &nv, &nd, &start_index, &agg_scp);
  std::cout<< field_name <<" "<<nv<<" "<<nd<<" "<<start_index<<" " << agg_scp << "\n";

  assert(!m3dc1_mesh::instance()->field_exists(fid));
  m3dc1_field_create(&fid, field_name, &nv, &nd, &agg_scp);
  m3dc1_field * mf = m3dc1_mesh::instance()->get_field(fid);

  int num_dof = mf->getDofsPerNode();
#ifdef PETSC_USE_COMPLEX
  num_dof *= 2;
#endif
  double * dof_data = new double[num_dof];
  apf::MeshEntity * e = NULL;
  apf::MeshIterator * it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    fscanf(fp, "%d %d %d\n", &gid, &lid, &ndof);
    assert(get_ent_globalid(m,e)+start_index==gid && getMdsIndex(m, e)==lid);
    memset(&dof_data[0],0,num_dof*sizeof(double));
    for (int i = 0; i < ndof; ++i)
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

void write_field(apf::Mesh2* m, m3dc1_field* mf, const char* filename, int start_index)
{
  std::string in(filename);
  std::stringstream s;
  s << in << "-" << PCU_Comm_Self();
  std::string partFile = s.str();
  FILE * fp = fopen(partFile.c_str(), "w");
  assert(fp);
  int num_dof = mf->getDofsPerNode();
#ifdef PETSC_USE_COMPLEX
  num_dof *= 2;
#endif
  double * dof_data = new double[num_dof];

  fprintf(fp, "%s %d %d %d %d\n",
          mf->getName().c_str(),
          mf->getBlocksPerNode(),
          mf->getDofsPerBlock(),
          start_index,
          mf->getAggregationScope());

  apf::MeshEntity * e = NULL;
  apf::MeshIterator * it = m->begin(0);
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

