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
#include "apfNumbering.h"

// todo: template to work with complex
class m3dc1_field
{
public:
  enum aggregation_scope
  {
    LOCAL_AGGREGATION = 0,
    PLANE_AGGREGATION = 1,
    GLOBAL_AGGREGATION = 2
  };
public:
  m3dc1_field(int i,
              const char * str,
              int blks_per_nd,
              int dofs_per_blk,
              aggregation_scope scp);
  ~m3dc1_field();
  int getId() { return id; }
  const std::string & getName() { return name; }
  apf::Field * getCoreField() { return fld; }
  apf::Numbering * getGlobalNumbering() { return gbl_num; }
  apf::Numbering * getLocalNumbering() { return lcl_num; }
  int getBlocksPerNode() { return blks_per_nd; }
  int getDofsPerBlock() { return dofs_per_blk; }
  int getDofsPerNode() { return blks_per_nd * dofs_per_blk; }
  int getAggregationScope() { return agg_scp; }
  void zero() { apf::zeroField(fld); }
private:
  int id;
  std::string name;
  apf::Field * fld;
  apf::Numbering * gbl_num;
  apf::Numbering * lcl_num;
  int agg_scp;
  int blks_per_nd;
  int dofs_per_blk;
  int dofs_per_nd;
};

void get_ent_xxxdofid(bool lcl, m3dc1_field * mf, int ent_lid, int * dof_id, int * dof_cnt);
void get_ent_localdofid(m3dc1_field* mf, int ent_lid, int* dof_id, int* dof_cnt);
void get_ent_globaldofid(m3dc1_field* mf, int ent_lid, int* dof_id, int* dof_cnt);
void get_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data);
void set_ent_dofdata(m3dc1_field* mf, apf::MeshEntity* e, double* dof_data);

// given a field ID, create a field out of file
void load_field(apf::Mesh2* m, int field_id, const char* filename);
void write_field(apf::Mesh2* m, m3dc1_field* mf, const char* filename, int start_index);
#endif
