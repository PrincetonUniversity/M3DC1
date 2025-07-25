/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include "m3dc1_mesh.h"
#include "m3dc1_field.h"
#include <mpi.h>
#include <PCU.h>
#include "pcu_util.h"
#include "gmi_null.h" // FIXME: should be deleted later on since it's added temporarily for null model
#include <gmi_analytic.h>
#include <map>
#include <cstring>
#include <iomanip> // setprecision
#include <fstream> // file input
#include "apfMDS.h"
#include "Expression.h"
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "pumi.h"
#include "apfShape.h" 
#include <spr.h>
// #include <stdio.h>
// #include <stdlib.h>
#ifdef M3DC1_TRILINOS
#include "m3dc1_ls.h"
#endif
#include <alloca.h>
#include "lionPrint.h"

const int dofNode = C1TRIDOFNODE;

#ifdef DEBUG
static const char* get_field_name_from_id(const FieldID id)
{
   m3dc1_field* mf = (*m3dc1_mesh::instance()->field_container)[id];
   return getName(mf->get_field());
}
#endif


/* static functions used for spr-adapt */
static apf::Field* get_field_at_index(apf::Mesh2* m, apf::Field* inField, int index, int numDofs)
{
  int numComps = apf::countComponents(inField);
  int numFields = numComps/numDofs;
  PCU_ALWAYS_ASSERT(index <= numFields);
  PCU_ALWAYS_ASSERT(index > 0);

  apf::Field* targetField = apf::createPackedField(m, "target_field", numDofs);

  apf::NewArray<double> allDofs(numComps);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    apf::getComponents(inField, v, 0, &(allDofs[0]));
    apf::setComponents(targetField, v, 0, &(allDofs[(index-1)*numDofs]));
  }
  m->end(it);
  return targetField;
}

static apf::Field* get_scalerComponent_of_field(apf::Mesh2* m, apf::Field* in, int comp)
{
  int comps = countComponents(in);
  PCU_ALWAYS_ASSERT(comp < comps);
  apf::Field* out = apf::createFieldOn(m, "in_field_comp", apf::SCALAR);

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);

  double dofs[FIXSIZEBUFF];
  while ( (e = m->iterate(it)) )
  {
    getComponents(in, e, 0, &dofs[0]);
    double dof_comp = dofs[comp];
    setComponents(out, e, 0, &dof_comp);
  }
  m->end(it);
  return out;
}

static apf::Field* get_vectorComponent_of_field(apf::Mesh2* m, apf::Field* in)
{
  apf::Field* out = apf::createFieldOn(m, "in_field_comp", apf::VECTOR);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);

  double dofs[FIXSIZEBUFF];
  while ( (e = m->iterate(it)) )
  {
    getComponents(in, e, 0, &dofs[0]);
    double comp_dR = dofs[1];
    double comp_dZ = dofs[2];
    double comp_dR2 = dofs[3];
    double comp_dZ2 = dofs[5];
    double vField[3] = {comp_dR, comp_dZ, 0.0};
    setComponents(out, e, 0, &vField[0]);
  }
  m->end(it);
  return out;
}

static apf::Field* get_ip_field(apf::Mesh2* m, apf::Field* in)
{
  ReducedQuinticImplicit shape;
  int numComps = apf::countComponents(in);
  assert(numComps == dofNode);
  int dim = m->getDimension();
  assert(dim == 2);
  int order = 2;
  apf::Field* ip = apf::createIPField(m, "ip_field", apf::VECTOR, order);

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);

  while ( (e = m->iterate(it)) )
  {
    /* setup the ReducedQuintic Related Info */
    apf::MeshEntity* dvs[3];
    int nd = m->getDownward(e, 0, dvs);
    double coords[3][2];
    for (int i = 0; i < 3; i++) {
      apf::Vector3 p;
      m->getPoint(dvs[i], 0, p);
      coords[i][0] = p[0];
      coords[i][1] = p[1];
    }
    shape.setCoord(coords);

    apf::NewArray<double> values(3*dofNode);
    apf::NewArray<double> dofAtXi(dofNode);
    for (int i = 0; i < 3; i++)
      apf::getComponents(in, dvs[i], 0, &(values[dofNode*i]));

    shape.setDofs(&(values[0]));


    apf::MeshElement* me = apf::createMeshElement(m, e);
    for (int i = 0; i < apf::countIntPoints(me, order); i++) {
      apf::Vector3 xi; /* parametric coords of the point in e at which we are evaluating the field */
      apf::getIntPoint(me, order, i, xi);
      apf::Vector3 p; /* physical coords of the point in e at which we are evaluating the field */
      apf::mapLocalToGlobal(me, xi, p);
      double pArray[3];
      p.toArray(pArray);
      shape.eval_g(pArray, &(dofAtXi[0]));
      apf::Vector3 grad(dofAtXi[1], dofAtXi[2], 0.0);
      apf::setVector(ip, e, i, grad);
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);
  return ip;
}

static void process_size_field(apf::Mesh2* m, apf::Field* in_size, int ts,
    double max_size, int refine_level, int coarsen_level)
{
  /* compute both average and min current size at each vertex */
  apf::Field* sum_field = apf::createFieldOn(m, "sum_field", apf::SCALAR);
  apf::Field* min_field = apf::createFieldOn(m, "min_field", apf::SCALAR);
  apf::Field* cnt_field = apf::createFieldOn(m, "cnt_field", apf::SCALAR);

  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    double current_min = 1.e32;
    double current_sum = 0.;
    double current_cnt = 0.;
    for (int i = 0; i < m->countUpward(v); i++) {
      double edge_length = apf::measure(m, m->getUpward(v, i));
      if (edge_length < current_min)
      	current_min = edge_length;
      current_sum += edge_length;
      current_cnt += 1.;
    }
    apf::setScalar(sum_field, v, 0, current_sum);
    apf::setScalar(min_field, v, 0, current_min);
    apf::setScalar(cnt_field, v, 0, current_cnt);
  }
  m->end(it);
  /* accumulate sum and cnt fields and compute the averages */
  apf::accumulate(sum_field);
  apf::accumulate(cnt_field);
  it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    double total_sum = apf::getScalar(sum_field, v, 0);
    double total_cnt = apf::getScalar(cnt_field, v, 0);
    apf::setScalar(sum_field, v, 0, total_sum/total_cnt);
  }
  m->end(it);

  /* update the min field using share reduction with min op */
  apf::sharedReduction(min_field, 0, false, apf::ReductionMin<double>());

  /* compute min/max of avg_size over the whole mesh */
  double mesh_min = 1.e16;
  double mesh_max = -1.e16;
  it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    double s = apf::getScalar(sum_field, v, 0);
    if (s > mesh_max)
      mesh_max = s;
    if (s < mesh_min)
      mesh_min = s;
  }
  m->end(it);

  PCU_Min_Doubles(&mesh_min, 1);
  PCU_Max_Doubles(&mesh_max, 1);

  if (!PCU_Comm_Self()) {
    printf("min/max of current avg size at time step %d: %f/%f\n", ts, mesh_min, mesh_max);
    printf("user requested max_size at time step %d: %f\n", ts, max_size);
  }

  int bdim = m->getDimension() - 1;

  double refine_factor = 1.;
  for (int i = 0; i < refine_level; i++)
    refine_factor *= 2.;

  double coarsen_factor = 1.;
  for (int i = 0; i < coarsen_level; i++)
    coarsen_factor *= 2.;

  it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    double asked_size = apf::getScalar(in_size, v, 0);
    double avg_size = apf::getScalar(sum_field, v, 0);
    double min_size = apf::getScalar(min_field, v, 0);
    int mtype = m->getModelType(m->toModel(v));
    if (mtype == bdim) {
      asked_size = avg_size;
    }
    else {
      if (asked_size < min_size / refine_factor) /* cap refinement by refine_factor (:=2^refine_level) */
      	asked_size = min_size / refine_factor;
      if (asked_size > min_size * coarsen_factor) /* cap coarsening by coarsen_factor (:=2^coarsen_level) */
      	asked_size = min_size * coarsen_factor;
    }
    /* cap the biggest size to user specified max_size,
 *        thus never allowing the mesh to get coarser that max_size */
    if (asked_size > max_size)
      asked_size = max_size;
    apf::setScalar(in_size, v, 0, asked_size);
  }
  m->end(it);
  apf::synchronize(in_size);

  /* compute min/max of avg_size over the whole mesh */
  double asked_min = 1.e16;
  double asked_max = -1.e16;
  it = m->begin(0);
  while ( (v = m->iterate(it)) )
  {
    double s = apf::getScalar(in_size, v, 0);
    if (s > asked_max)
      asked_max = s;
    if (s < asked_min)
      asked_min = s;
  }
  m->end(it);

  PCU_Min_Doubles(&asked_min, 1);
  PCU_Max_Doubles(&asked_max, 1);


  if (!PCU_Comm_Self())
    printf("min/max of asked size at time step %d: %f/%f\n", ts, asked_min, asked_max);

  /* clean up */
  m->removeField(sum_field);
  m->removeField(min_field);
  m->removeField(cnt_field);
  apf::destroyField(sum_field);
  apf::destroyField(min_field);
  apf::destroyField(cnt_field);
}

/* for 3D
 *   typedefs */
typedef std::vector<apf::Field*> MultiField; // fields on different planes

/* the naming convention of multi-plane fields is "name_pnnn" where
 *    "name" is the original name of the field
 *       "nnn" is the plane number, zero-padded to have length equal to zero_pad_length */
const int zero_pad_length = 3; /* this would allow number of planes to be between 0 and 999 */
const int max_num_plane = 999; /* this should be equal (10^zero_pad_length - 1) */

static const char* get_plane_name_format(int pad = zero_pad_length)
{
  PCU_ALWAYS_ASSERT(pad >= 3 && pad < 7);
  static const char* format[7] =
  {
    "", /* 0 */
    "", /* 1 */
    "", /* 2 */
    "%s_p%03d", /* 3 */
    "%s_p%04d", /* 3 */
    "%s_p%05d", /* 3 */
    "%s_p%06d", /* 3 */
  };
  return format[pad];
}

static void get_original_field_name(const char* plane_name, char* original_name)
{
  std::string plane_name_str(plane_name);
  std::size_t s = plane_name_str.size();
  s -= (zero_pad_length+2);
  plane_name_str.resize(s);

  std::strcpy(original_name, plane_name_str.c_str());
}

/* helper functions */
static bool is_in_plane(apf::MeshEntity* e)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);

  double tol = 2. * M3DC1_PI / m3dc1_model::instance()->num_plane / 1.e6;
  apf::Vector3 p[2];
  apf::MeshEntity* v[2];
  m->getDownward(e, 0, v);
  for (int i = 0; i < 2; i++) {
    m->getPoint(v[i], 0, p[i]);
  }
  return std::fabs(p[0][2] - p[1][2]) < tol;
}

static void destroyElement(apf::Mesh2* m, apf::MeshEntity* e)
{
  int dim = apf::getDimension(m,e);
  if (dim < m->getDimension())
  { /* destruction is a no-op if this entity still supports
       higher-order ones */
    if (m->hasUp(e))
      return;
  }
  if (dim < m->getDimension())
  {
    int etype = m->getType(e);
    if (etype == apf::Mesh::TRIANGLE) return;
    if (etype == apf::Mesh::EDGE)
    {
      if (is_in_plane(e)) return;
    }
    /* if (etype == apf::Mesh::VERTEX) */
  }
  apf::Downward down;
  int nd = 0;
  if (dim > 0)
    nd = m->getDownward(e,dim-1,down);
  m->destroy(e);
  /* destruction applies recursively to the closure of the entity */
  if (dim > 0)
    for (int i=0; i < nd; ++i)
      destroyElement(m,down[i]);
}

static void remove_all_wedges()
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;

  apf::MeshEntity* e;
  apf::MeshIterator* it;

  it = m->begin(3);
  while ( (e = m->iterate(it)) )
  {
    destroyElement(m, e);
  }
  m->acceptChanges();
  changeMdsDimension(m, 2);
}



static void transfer_field_data_on_plane(apf::Field* in, apf::Field* out, int p)
{
  int np = m3dc1_model::instance()->num_plane;
  int local_planeid = m3dc1_model::instance()->local_planeid;
  PCU_ALWAYS_ASSERT_VERBOSE(p >= 0, "p must be strictly positive!");
  PCU_ALWAYS_ASSERT_VERBOSE(p < np, "p must be strictly less than number of planes!");
  int nc = apf::countComponents(in);
  apf::NewArray<double> dofs(nc);
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ( (e = m->iterate(it)) )
  {
    /* if not on the plane continue */
    if (p != local_planeid) continue;
    /* if not owned continue
 *        if (!m->isOwned(e)) continue; */
    apf::getComponents(in, e, 0, &dofs[0]);
    apf::setComponents(out, e, 0, &dofs[0]);
  }
  m->end(it);
  synchronize_field(out);
}

static void move_field_data_down(apf::Field* f, int startp)
{
  int np = m3dc1_model::instance()->num_plane;
  int local_planeid = m3dc1_model::instance()->local_planeid;
  PCU_ALWAYS_ASSERT_VERBOSE(startp > 0, "startp must be strictly positive!");
  PCU_ALWAYS_ASSERT_VERBOSE(startp < np, "startp must be strictly less than number of planes!");
  int nc = apf::countComponents(f);
  apf::NewArray<double> zeros(nc);
  for (int i = 0; i < nc; i++)
    zeros[i] = 0.;


  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  while ( (e = m->iterate(it)) )
  {
    /* continue if the edge is in any of the poloidal planes */
    if (is_in_plane(e)) continue;
    /* continue if the edge is not on local plane */
    if (startp-1 != local_planeid) continue;
    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    apf::Vector3 ps[2];
    for (int i = 0; i < 2; i++)
      m->getPoint(vs[i], 0, ps[i]);
    /* make sure the second vert in vs is on lower plane (i.e. has lower z component) */
    if (ps[1][2] > ps[0][2])
    {
      std::swap(ps[0], ps[1]);
      std::swap(vs[0], vs[1]);
    }
    apf::NewArray<double> dofs(nc);
    apf::getComponents(f, vs[0], 0, &dofs[0]);
    apf::setComponents(f, vs[1], 0, &dofs[0]);
    apf::setComponents(f, vs[0], 0, &zeros[0]);
  }
  m->end(it);
  synchronize_field(f);
  /* PCU_Barrier(); */
}

static void transfer_field_to_main(apf::Field* f, MultiField& pfields)
{
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int np = m3dc1_model::instance()->num_plane;
  int nc = apf::countComponents(f);

  PCU_ALWAYS_ASSERT(np <= max_num_plane);

  for (int i = 0; i < np; i++) {
    char pfieldname[128];
    sprintf(pfieldname, get_plane_name_format(), apf::getName(f), i);
    apf::Field* tmp = createPackedField(m, pfieldname, nc, apf::getShape(f));
    apf::zeroField(tmp);
    pfields.push_back(tmp);
  }


  for (int i = 0; i < (int)pfields.size(); i++)
    transfer_field_data_on_plane(f, pfields[i], i);


  for (int i = 1; i < (int)pfields.size(); i++)
    for (int j = 0; j < i; j++)
      move_field_data_down(pfields[i], i-j);
}

void transfer_field_from_main(MultiField& pfields)
{
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int np = m3dc1_model::instance()->num_plane;
  PCU_ALWAYS_ASSERT(np == (int)pfields.size());
  int lpid = m3dc1_model::instance()->local_planeid;
  int nc = apf::countComponents(pfields[0]);

  /* get the original field corresponding to multiplane pfields */
  char original_name[128];
  get_original_field_name(apf::getName(pfields[0]), original_name);
  apf::Field* f = m->findField(original_name);
  PCU_ALWAYS_ASSERT_VERBOSE(f, "was not able to find the original field!");


  /* first transfer everything to the last plane */
  it = m->begin(1);
  while ( (e = m->iterate(it)) )
  {
    if (lpid != np-1) continue;
    if (is_in_plane(e)) continue;
    apf::MeshEntity* v[2];
    apf::Vector3 p[2];
    m->getDownward(e, 0, v);
    for (int i = 0; i < 2; i++)
      m->getPoint(v[i], 0, p[i]);
    if (p[0][2] > p[1][2])
    {
      std::swap(v[0], v[1]);
      std::swap(p[0], p[1]);
    }
    apf::NewArray<double> dofs(nc);
    for (int j = 0; j < (int)pfields.size(); j++) {
      apf::getComponents(pfields[j], v[0], 0, &dofs[0]);
      apf::setComponents(pfields[j], v[1], 0, &dofs[0]);
    }
  }
  m->end(it);

  for (int i = 0; i < (int)pfields.size(); i++)
    synchronize_field(pfields[i]);

  for (int i = 1; i < np-1 ; i++)
    for (int j = np-1; j > i; j--)
      move_field_data_down(pfields[i], j);

  for (int i = 0; i < (int)pfields.size(); i++)
    transfer_field_data_on_plane(pfields[i], f, i);
}

static apf::Field* compute_multiplane_size_field(const MultiField& pfields, double ar)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int np = m3dc1_model::instance()->num_plane;
  PCU_ALWAYS_ASSERT(np == (int)pfields.size());

  int nc = apf::countComponents(pfields[0]);

  apf::MeshEntity* e;
  apf::MeshIterator* it;

  MultiField sizefields; // size fields computed for each plane

  for (int i = 0; i < np; i++) {
    apf::Field* ip = get_ip_field(m, pfields[i]);
    apf::Field* size = spr::getSPRSizeField(ip, ar);
    char szname[128];
    sprintf(szname, "sz_%s", apf::getName(pfields[i]));
    apf::Field* sz = apf::createField(m, szname, apf::SCALAR, apf::getShape(size));
    apf::copyData(sz, size);
    sizefields.push_back(sz);
    m->removeField(size);
    apf::destroyField(size);
    m->removeField(ip);
    apf::destroyField(ip);
    /* if (i != np-1) { */
    /*   m->removeField(ip); */
    /*   apf::destroyField(ip); */
    /* } */
  }

  apf::Field* mastersize = apf::createField(m, "size", apf::SCALAR, apf::getShape(sizefields[0]));

  it = m->begin(0);
  while ( (e = m->iterate(it)) )
  {
    double minsize = 1.e16;
    for (int i = 0; i < np; i++) {
      double s = apf::getScalar(sizefields[i], e, 0);
      if (s < minsize)
        minsize = s;
    }
    apf::setScalar(mastersize, e, 0, minsize);
  }
  m->end(it);
  /* apf::writeVtkFiles("03_mesh_with_all_sizes", m); */
  for (int i = 0; i < np; i++) {
    m->removeField(sizefields[i]);
    apf::destroyField(sizefields[i]);
  }
  synchronize_field(mastersize);
  /* apf::writeVtkFiles("04_mesh_with_master_size", m); */
  return mastersize;
}

static void zero_fields_on_non_master(const MultiField& pfields)
{
  int np = m3dc1_model::instance()->num_plane;
  PCU_ALWAYS_ASSERT(np == (int)pfields.size());

  int nc = apf::countComponents(pfields[0]);
  int local_planeid = m3dc1_model::instance()->local_planeid;

  double tol = 2. * M3DC1_PI / np / 1.e6;

  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;

  /* get the original field corresponding to multiplane pfields */
  char original_name[128];
  get_original_field_name(apf::getName(pfields[0]), original_name);
  apf::Field* f = m->findField(original_name);
  PCU_ALWAYS_ASSERT_VERBOSE(f, "was not able to find the original field!");


  apf::NewArray<double> zeros(nc);
  for (int i = 0; i < nc; i++)
    zeros[i] = 0.;

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ( (e = m->iterate(it)) )
  {
    apf::setComponents(f, e, 0, &zeros[0]);

    if (local_planeid == 0) continue;

    for (int j = 0; j < (int)pfields.size(); j++)
      apf::setComponents(pfields[j], e, 0, &zeros[0]);
  }
  m->end(it);

  it = m->begin(1);
  while ( (e = m->iterate(it)) )
  {
    if (local_planeid != 0) continue;
    if (is_in_plane(e)) continue;
    apf::MeshEntity* v[2];
    apf::Vector3 p[2];
    m->getDownward(e, 0, v);
    for (int i = 0; i < 2; i++) {
      m->getPoint(v[i], 0, p[i]);
    }

    if (std::fabs(p[0][2]) > tol)
    {
      for (int j = 0; j < (int)pfields.size(); j++)
	apf::setComponents(pfields[j], v[0], 0, &zeros[0]);
    }
    else
    {
      PCU_ALWAYS_ASSERT(std::fabs(p[1][2]) > tol);
      for (int j = 0; j < (int)pfields.size(); j++)
	apf::setComponents(pfields[j], v[1], 0, &zeros[0]);
    }
  }
  m->end(it);
  for (int j = 0; j < np; j++)
    apf::synchronize(pfields[j]);
}

/* the following should be called only on master plane 0 */
static void zero_fields_on_master(apf::Field* f) 
{
  int np = m3dc1_model::instance()->num_plane;
  int nc = apf::countComponents(f);
  int local_planeid = m3dc1_model::instance()->local_planeid;
  PCU_ALWAYS_ASSERT_VERBOSE(local_planeid==0, "function can only be called on master palne!");

  double tol = 2. * M3DC1_PI / np / 1.e6;

  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;

  apf::NewArray<double> zeros(nc);
  for (int i = 0; i < nc; i++)
    zeros[i] = 0.;

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ( (e = m->iterate(it)) )
    apf::setComponents(f, e, 0, &zeros[0]);
  m->end(it);
  for (int j = 0; j < np; j++)
    apf::synchronize(f);
}
/* end of static functions used for spr-adapt */

#ifdef DEBUG
int begin_numVert;
#endif
double begin_mem, begin_time;

// helper routines
void group_complex_dof (apf::Field* field, int option);

void print_ent_error(const char* func, int dim, int id)
{
  std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<": "
          <<func<<" entity dim "<<dim<<" id "<<id<<" doesn't exist\n";
}

bool m3dc1_double_isequal(double A, double B)
{
  double maxDiff = 1e-5;
  double maxRelDiff = 1e-5;
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

//*******************************************************
int m3dc1_scorec_init()
//*******************************************************
{ 
  pumi_start();
  begin_time=MPI_Wtime();
  return M3DC1_SUCCESS; 
}

//*******************************************************
void m3dc1_scorec_verbosity(int* l)
//*******************************************************
{
  lion_set_verbosity(*l);
}

//*******************************************************
int m3dc1_scorec_finalize()
//*******************************************************
{ 
  pumi_mesh_deleteGlobalID(m3dc1_mesh::instance()->mesh);  // delete global id
  m3dc1_mesh::instance()->clean(); // delete tag, field and internal data 

#ifdef DEBUG
  int local_numVert=m3dc1_mesh::instance()->mesh->count(0);
  int global_numVert=0;
  MPI_Allreduce(&local_numVert, &global_numVert, 
                4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (begin_numVert != global_numVert)
    m3dc1_mesh::instance()->mesh->writeNative("adapted-mesh.smb");
#endif

  pumi_mesh_delete(m3dc1_mesh::instance()->mesh);

  if (!pumi_rank()) 
    std::cout<<"\n* [M3D-C1 INFO] run time: "<<MPI_Wtime()-begin_time<<" (sec)\n";
  pumi_finalize();
  return M3DC1_SUCCESS; 
}

/** plane functions */

//*******************************************************
int m3dc1_plane_setnum(int* num_plane)
//*******************************************************
{
  if (m3dc1_mesh::instance()->mesh->getDimension()==3)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: #plane should be set before constructing 3D mesh\n";
    return M3DC1_FAILURE; 
  }

  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) 
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: invalid #planes - mod(#ranks, #planes) should be 0\n";
    return M3DC1_FAILURE; 
  }

  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_plane_getnum(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_plane_getid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphirange(double* min_val, double* max_val)
//*******************************************************
{
  if (m3dc1_mesh::instance()->mesh->getDimension()==3)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: phi should be set before constructing 3D mesh\n";
    return M3DC1_FAILURE; 
  }

  m3dc1_model::instance()->set_phi(*min_val, *max_val);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphi(int* planeid, double* phi)
//*******************************************************
{
  if (m3dc1_mesh::instance()->mesh->getDimension()==3)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: phi should be set before constructing 3D mesh\n";
    return M3DC1_FAILURE; 
  }

  m3dc1_model::instance()->set_phi(*planeid, *phi);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_getphi(int* planeid, double* phi)
//*******************************************************
{
  m3dc1_model::instance()->get_phi(*planeid, phi);
  return M3DC1_SUCCESS; 
}

/** model functions */
//*******************************************************
int m3dc1_model_getplaneid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getmincoord(double* x_min, double* y_min)
//*******************************************************
{
  *x_min = m3dc1_model::instance()->boundingBox[0];
  *y_min = m3dc1_model::instance()->boundingBox[1];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getmaxcoord(double* x_max, double* y_max)
//*******************************************************
{
  *x_max = m3dc1_model::instance()->boundingBox[2];
  *y_max = m3dc1_model::instance()->boundingBox[3];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_load(char* /* in */ model_file)
//*******************************************************
{ 
  FILE *test_in = fopen (model_file,"r");
  if (!test_in)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: model file \""<<model_file<<"\" doesn't exist\n";
    return M3DC1_FAILURE; 
  }
  else
    fclose(test_in);

  pGeom g;
  std::string model_name = model_file; 
  if(model_name.substr(model_name.find_last_of(".") + 1) == "dmg")
  { // modelType = 2 (.dmg model)
    g = pumi_geom_load(model_file, "mesh");
    m3dc1_model::instance()->model = g->getGmi();   
    gmi_register_mesh();
    m3dc1_model::instance()->model = gmi_load(model_file);
    m3dc1_model::instance()->modelType = 2;
  }
  else
  { // modelType = 1 (non .dmg models)
    g = pumi_geom_load(model_file, "analytic");
    m3dc1_model::instance()->model = g->getGmi();
    m3dc1_model::instance()->load_analytic_model(model_file); 
    m3dc1_model::instance()->modelType = 1;
  }
  pumi_geom_freeze(g);

  if (m3dc1_model::instance()->modelType == 1)
    m3dc1_model::instance()->caculateBoundingBox();

  // save the num of geo ent on the oringal plan
  m3dc1_model::instance()->numEntOrig[0]=m3dc1_model::instance()->model->n[0];
  m3dc1_model::instance()->numEntOrig[1]=m3dc1_model::instance()->model->n[1];
  m3dc1_model::instance()->numEntOrig[2]=m3dc1_model::instance()->model->n[2];
  return M3DC1_SUCCESS;
}
//*******************************************************
void m3dc1_model_getmodeltype(int* modeltype)
//*******************************************************
{
  // 1 for .txt model, 2 for .dmg model
  *modeltype = m3dc1_model::instance()->modelType;
}

//*******************************************************
int m3dc1_modelinfo_load(char* /* in */ modelInfo_file)
//*******************************************************
{
  if (m3dc1_model::instance()->modelType == 1)
    return M3DC1_SUCCESS;  // only move forward for .dmg models

  assert(m3dc1_model::instance()->modelType == 2);
  if (std::strcmp(modelInfo_file, "dummyInfo") == 0)
  {
    std::cout << "Make sure to provide model_info file since model is .dmg\n";
    exit(1);  // Need to have model info file
  }

  std::cout << "Model Info file = " << modelInfo_file << "\n";
  
  // Check if the given model info file exists or not:
  std::ifstream ifs(modelInfo_file); 
  if (!ifs) 
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: model file \""<<modelInfo_file<<"\" doesn't exist\n";
    return M3DC1_FAILURE; 
  }

  // Load the content of the model info file
  // First read the bounding box
  double r_min, z_min, r_max, z_max;
  ifs >> r_min >> z_min >> r_max >> z_max;
  m3dc1_model::instance()->boundingBox[0]= r_min;
  m3dc1_model::instance()->boundingBox[1]= z_min;
  m3dc1_model::instance()->boundingBox[2]= r_max;
  m3dc1_model::instance()->boundingBox[3]= z_max;
  
  // Read innermost loop
  int numEdges = 0;
  ifs >> numEdges;
  for (int i = 0; i < numEdges; i++)
  {
    int edgeId;
    ifs >> edgeId;
    m3dc1_model::instance()->innerLoop.push_back(edgeId);
  } 
  
  // Read outermost loop
  ifs >> numEdges;
  for (int i = 0; i < numEdges; i++)
  {
    int edgeId;
    ifs >> edgeId;
    m3dc1_model::instance()->outerLoop.push_back(edgeId);
  } 
}

// ****************************************************
void m3dc1_model_getgeometricloop(int* modelEdges, int* numberOfEdges, int* loopType)
// ****************************************************
{
  std::vector <int> geomEdges;
  if (*loopType == 0)
    geomEdges = m3dc1_model::instance()->innerLoop;
  else if (*loopType == 1)
    geomEdges = m3dc1_model::instance()->outerLoop;
  else  // might add further types in future
  {
    std::cout << "Incorrect Loop Type. Set 0 for inner loop, and 1 for outer loop\n";
    exit(1);
  }

  // Print a warning if the requested loop does not exist.
  if (geomEdges.size() == 0)
    std::cout << "Warning: The requested loop has zero edges - Loop does not exist in the model\n";

  for (int i = 0; i <  geomEdges.size(); i++)
    modelEdges[i] = geomEdges[i];

  *numberOfEdges = geomEdges.size();
}

//*******************************************************
int m3dc1_model_print()
//*******************************************************
{
  if (PCU_Comm_Self() || m3dc1_model::instance()->local_planeid) 
    return M3DC1_SUCCESS;

  gmi_iter* gf_it = gmi_begin(m3dc1_model::instance()->model, 2);
  gmi_ent* ge;

  while ((ge = gmi_next(m3dc1_model::instance()->model, gf_it))) 
  {
    gmi_set* gf_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())    std::cout<<"model face id "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - # adj model edges: "<<gf_edges->n<<"\n";
    if (gf_edges->n)
    {
    if (!PCU_Comm_Self())      std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gf_edges->n; ++i)
    if (!PCU_Comm_Self())        std::cout<<gmi_tag(m3dc1_model::instance()->model,  gf_edges->e[i])<<" "; 
    if (!PCU_Comm_Self())      std::cout<<"\n";
    }
    gmi_free_set(gf_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gf_it);

// verify geom edge
  gmi_iter* ge_it = gmi_begin(m3dc1_model::instance()->model, 1);
  while ((ge = gmi_next(m3dc1_model::instance()->model, ge_it))) 
  {
    M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,ge);
    if (!pn)   
    {
    if (!PCU_Comm_Self())       std::cout<<"["<<PCU_Comm_Self()<<"] model edge "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - failed with gmi_analytic_data retrieval\n";
    }
  }
  gmi_end(m3dc1_model::instance()->model, ge_it);

// verify gv_edges
  gmi_iter* gv_it = gmi_begin(m3dc1_model::instance()->model, 0);
  while ((ge = gmi_next(m3dc1_model::instance()->model, gv_it))) 
  {
    gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())
      std::cout<<"model vertex id "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - # adj model edges: "<<gv_edges->n<<"\n";
    if (gv_edges->n)
    {
    if (!PCU_Comm_Self())
      std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gv_edges->n; ++i)
            if (!PCU_Comm_Self()) std::cout<<gmi_tag(m3dc1_model::instance()->model,  gv_edges->e[i])<<" "; 
          if (!PCU_Comm_Self()) std::cout<<"\n";
    }
    gmi_free_set(gv_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gv_it);

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_setnumplane(int* num_plane)
//*******************************************************
{
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) return M3DC1_FAILURE;
  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getnumplane(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_model_settopo()
//*******************************************************
{
  m3dc1_model::instance()->snapping = true;
}

/** mesh functions */
#include <parma.h>

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    m->setDoubleTag(e, tag, &w);
  m->end(it);
}

apf::MeshTag* setWeights(apf::Mesh* m) {
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}

void clearTags(apf::Mesh* m, apf::MeshTag* t) {
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}

void m3dc1_mesh_load_3d(char* mesh_file, int* num_plane)
{
  // switch COMM to GLOBAL COMM
  MPI_Comm groupComm = PCU_Get_Comm();
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);

  m3dc1_model::instance()->create3D();

  m3dc1_mesh::instance()->mesh = apf::loadMdsMesh(m3dc1_model::instance()->model,
            mesh_file);

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  // delete numbering and tags loaded from file
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    destroyNumbering(n);
  }

  apf::DynamicArray<apf::MeshTag*> tags;
  mesh->getTags(tags);
  for (int i=0; i<tags.getSize(); i++)
  {
    if (mesh->findTag("norm_curv")==tags[i]) continue;
    for (int idim=0; idim<4; idim++)
      apf::removeTagFromDimension(mesh, tags[i], idim);
    mesh->destroyTag(tags[i]);
  }
  m3dc1_mesh::instance()->initialize();
  if (m3dc1_model::instance()->num_plane==1) // 2D problem
  {
    compute_globalid(m3dc1_mesh::instance()->mesh, 0);
    compute_globalid(m3dc1_mesh::instance()->mesh, 2);
  }
  m3dc1_mesh::instance()->set_node_adj_tag();
}

#include <sstream>
//*******************************************************
int m3dc1_mesh_load(char* mesh_file)
//*******************************************************
{ 
  std::string in(mesh_file);
  size_t p = in.rfind('.');
  std::string base = in.substr(0,p);
  std::string ext = in.substr(p+1);
  std::stringstream s;
  s << base << pumi_rank()<<"." << ext;
  std::string partFile = s.str(); 
  FILE* test_in = fopen (partFile.c_str(),"r");

  if (!test_in)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: mesh file \""<<mesh_file<<"\" doesn't exist\n";
    return M3DC1_FAILURE; 
  }
  else
    fclose(test_in);

  if (m3dc1_model::instance()->local_planeid == 0) // master plane
  {
    m3dc1_mesh::instance()->mesh = pumi_mesh_load(pumi::instance()->model, 
                                    mesh_file, m3dc1_model::instance()->group_size);

    /* vertex load balancing */
    //Parma_PrintPtnStats(m3dc1_mesh::instance()->mesh, "initial");

    // clean-up tag, field and numbering loaded from file
    apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
    while(mesh->countFields())
    {
      apf::Field* f = mesh->getField(0);
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": field "<<getName(f)<<" deleted\n";
      destroyField(f);
    }

    while(mesh->countNumberings())
    {
      apf::Numbering* n = mesh->getNumbering(0);
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
      destroyNumbering(n);
    }

    apf::DynamicArray<apf::MeshTag*> tags;
    mesh->getTags(tags);
    for (int i=0; i<tags.getSize(); i++)
    {
      if (mesh->findTag("norm_curv")==tags[i]) continue;
      for (int idim=0; idim<4; idim++)
        apf::removeTagFromDimension(mesh, tags[i], idim);
      mesh->destroyTag(tags[i]);
    }
  } 
  else // non-master plane
    m3dc1_mesh::instance()->mesh = pumi_mesh_create(pumi::instance()->model, 2, false);

  m3dc1_mesh::instance()->initialize();
  
  if (m3dc1_model::instance()->num_plane==1) // 2D problem
  {
    compute_globalid(m3dc1_mesh::instance()->mesh, 0);
    compute_globalid(m3dc1_mesh::instance()->mesh, 2);
  }

#ifdef DEBUG
  pumi_mesh_verify(m3dc1_mesh::instance()->mesh, false);
  begin_numVert=m3dc1_mesh::instance()->mesh->count(0);
  MPI_Allreduce(&local_numVert, &begin_numVert,
                4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_mesh_build3d (int* num_field, int* field_id,  
                        int* num_dofs_per_value)
//*******************************************************
{ 
  // switch COMM to GLOBAL COMM
  MPI_Comm groupComm = PCU_Get_Comm();
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);

  // initialize phi value and construct 3d
  int num_plane = m3dc1_model::instance()->num_plane;
  if (num_plane==1)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: #plane must be greater than 1\n";
    return M3DC1_FAILURE; 
  }

  if (!m3dc1_model::instance()->phi)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 INFO] "<<__func__
             <<": setting phi range (0.0, "<<2.0*M3DC1_PI/num_plane*(num_plane-1)<<")\n";
    m3dc1_model::instance()->set_phi(0.0, 2.0*M3DC1_PI/num_plane*(num_plane-1));
  }

  m3dc1_model::instance()->setupCommGroupsPlane();
  // returns error if num_plane>1
  pumi_mesh_deleteGlobalID(m3dc1_mesh::instance()->mesh);
  m3dc1_mesh::instance()->build3d(*num_field, field_id, num_dofs_per_value);

  // update global ID
  compute_globalid(m3dc1_mesh::instance()->mesh, 0);
  compute_globalid(m3dc1_mesh::instance()->mesh, 3);

#ifdef DEBUG
  pumi_mesh_verify(m3dc1_mesh::instance()->mesh, false);
  int local_numVert=m3dc1_mesh::instance()->mesh->count(0);
  MPI_Allreduce(&local_numVert, &begin_numVert, 
                4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  return M3DC1_SUCCESS; 
}

void m3dc1_dir_export(double* dir, int ts, int dir_size)
{
  char fname[64];
  sprintf(fname, "dirs-ts%d-%d", ts, PCU_Comm_Self());

  FILE * fp =fopen(fname, "w");

  fprintf(fp, "%d ", dir_size);
  for (int i=0; i<dir_size; ++i)
    fprintf(fp, "%lf\n", dir[i]);
  fclose(fp);
}

void m3dc1_dir_import(double* dir, int ts)
{
  char fname[64];
  sprintf(fname, "dirs-ts%d-%d", ts, PCU_Comm_Self());

  FILE * fp =fopen(fname, "r");
  assert(fp);
  int dir_size;
  fscanf(fp, "%d ", &dir_size);
  assert(dir_size==m3dc1_mesh::instance()->mesh->count(0)*3);

  for (int i=0; i<dir_size; ++i)
    fscanf(fp, "%lf", &(dir[i]));
  fclose(fp);
}

/* new mesh adaptation */
/* Input Parameters
 * field_id_h1, field_id_h2: removed before adaptation so it won't be available after adaptation
 * dir: direction per node. The length of dir should be #nodes * 3 
*/
void m3dc1_mesh_adapt(int* field_id_h1, int* field_id_h2, double* dir)
{
#ifdef DEBUG
  static int ts=1;
  if (!PCU_Comm_Self()) 
    std::cout<<"field_id_h1 "<<*field_id_h1
             <<", field_id_h2 "<<*field_id_h2<<"\n";
  m3dc1_field_export (field_id_h1, &ts);
  m3dc1_field_export (field_id_h2, &ts);
  m3dc1_dir_export(dir, ts, m3dc1_mesh::instance()->mesh->count(0)*3);
  ts++;
  // export fields and dirs
#endif
  adapt_mesh (*field_id_h1, *field_id_h2, dir);
#ifdef DEBUG
  printStats(m3dc1_mesh::instance()->mesh);
  pumi_mesh_verify(m3dc1_mesh::instance()->mesh, false);
#endif
}

/* ghosting functions */
//*******************************************************
int m3dc1_ghost_create (int* num_layer )
//*******************************************************
{
  pMesh m = m3dc1_mesh::instance()->mesh;

  // create layer
  pumi_ghost_createLayer (m, 0, pumi_mesh_getDim(m), *num_layer, 1);

  if (m3dc1_mesh::instance()->field_container)
  {
    std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
    while(it!=m3dc1_mesh::instance()->field_container->end())
    {
#ifdef DEBUG
      int fieldId= it->first, isnan;
      m3dc1_field_isnan(&fieldId, &isnan);
      assert(isnan==0);
#endif
      synchronize_field(it->second->get_field());
      ++it;
    }
  }
  if (!pumi_rank()) std::cout<<"[M3D-C1 INFO] "<<* num_layer<<" ghost layer(s) created\n";

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ghost_delete()
//*******************************************************
{
  pumi_ghost_delete (m3dc1_mesh::instance()->mesh);

  if (!pumi_rank()) std::cout<<"[M3D-C1 INFO] ghost layer(s) deleted\n";

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  *num_ent = m3dc1_mesh::instance()->mesh->count(*ent_dim);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  *num_ent = m3dc1_mesh::instance()->num_own_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  *num_ent = m3dc1_mesh::instance()->num_global_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumghostent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_mesh::instance()->mesh->count(*ent_dim) - 
             m3dc1_mesh::instance()->num_local_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_search(int* initial_simplex,
		      double* final_position,
		      int* final_simplex)
//*******************************************************
{
  bool located = false;
  apf::MeshEntity* e = NULL;
  apf::MeshEntity* simplex = NULL;
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::Adjacent adjacent;
  int edge_curr_index, edge_prev_index= -1;  
  int simplex_dim = m->getDimension();
  int vertex_dim = 0, edge_dim = 1;
  int bridge_dim = simplex_dim - 1;
  int count = 0, max_count = m->count(simplex_dim);
  double tol = 1e-15;
  
  simplex = apf::getMdsEntity(m, simplex_dim, *initial_simplex);
  while (not located) 
  {
    int simplex_index = apf::getMdsIndex(m, simplex);
    apf::Downward vertices;
    apf::Vector3 v_coord[3];
    apf::Matrix<3, 3> A;
    int nv = m->getDownward(simplex, vertex_dim, vertices);
    for (int j = 0; j < nv; ++j) {
      m->getPoint(vertices[j], 0, v_coord[j]);
      // Note: second argument is 0 for linear meshes
    }
    // Compute (linear) barycentric coordinates of final position
    // with respect to current simplex
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 3; ++k)
	A[j][k] = v_coord[k][j];
    for (int j = 0; j < 3; ++j)
      A[2][j] = 1;
    apf::Matrix<3, 3> Ainv = apf::invert(A);
    apf::Vector3 b_coords; // = Ainv * final_position
    for (int j = 0; j < 3; ++j) {
      b_coords[j] = (Ainv[j][0] * final_position[0] +
		     Ainv[j][1] * final_position[1] +
		     Ainv[j][2]);
    }
    // If all positive for current simplex, exit.
    if (((b_coords[0] >= -tol) && (b_coords[0] <= (1 + tol))) &&
	((b_coords[1] >= -tol) && (b_coords[1] <= (1 + tol))) &&
	((b_coords[2] >= -tol) && (b_coords[2] <= (1 + tol)))) {
      located = true;
      *final_simplex = simplex_index;
      break;
    }
    // Otherwise, check which coordinates are negative
    bool b_negative[3] = {false, false, false};
    int bneg_index[2] = {0, 0};
    int bneg_count = 0;
    for (int j = 0; j < 3; ++j)
      if (b_coords[j] < 0) 
      {
	b_negative[j] = true;
	bneg_index[bneg_count] = j;
	++bneg_count;
      }
    // Obtain the index of most negative coordinate
    assert(bneg_count > 0 && bneg_count < 3);

    // Ensure bneg_index[0] is the index of vertex whose corresponding
    // barycentric coordinate is most negative; ties automatically resolved
    // as a result.
    if (bneg_count == 2)
      if (fabs(b_coords[bneg_index[0]]) <
	  fabs(b_coords[bneg_index[1]])) {
	int tmp_index = bneg_index[0];
	bneg_index[0] = bneg_index[1];
	bneg_index[1] = tmp_index;
      }

    // Determine edge opposite to this vertex
    apf::MeshEntity* edge_vertices[2];
    int edge_count = 0;
    int edge_type = 1;    // Mesh entity type; see APF documentation.
    for (int j = 0; j < 3; ++j)
      if (j != bneg_index[0]) {
	edge_vertices[edge_count] = vertices[j];
	++edge_count;
      }
    
    // Find neighboring simplex sharing this edge
    e = apf::findElement(m, edge_type, edge_vertices);
    edge_curr_index = apf::getMdsIndex(m, e);

    // If current edge choice is same as previous edge, pick edge
    // opposite to second least (actual, not absolute, valued) barycentric
    // coordinate.
    if (edge_curr_index == edge_prev_index) {
      edge_count = 0;
      for (int j = 0; j < 3; ++j)
	if (j != bneg_index[1]) {
	  edge_vertices[edge_count] = vertices[j];
	  ++edge_count;
	}
      e = apf::findElement(m, edge_type, edge_vertices);
      edge_curr_index = apf::getMdsIndex(m, e);
    }     

    apf::getBridgeAdjacent(m, e,
			   bridge_dim, simplex_dim,
			   adjacent);
    if (adjacent.getSize() == 2) {
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
	int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
	if (new_simplex_index != simplex_index)
	  simplex = adjacent[j];
      }
    }
    else {
      apf::Downward edges;
      int ne = m->getDownward(simplex, edge_dim, edges);
      for (size_t j = 0; j < ne; ++j) {
	int edge_tmp_index = apf::getMdsIndex(m, edges[j]);
	if ((edge_tmp_index != edge_curr_index) &&
	    (edge_tmp_index != edge_prev_index)) {
	  e = edges[j];
	  break;
	}
      }
      edge_curr_index = apf::getMdsIndex(m, e);    
      apf::getBridgeAdjacent(m, e,
			     bridge_dim, simplex_dim,
			     adjacent);
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
	int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
	if (new_simplex_index != simplex_index)
	  simplex = adjacent[j];
      }
    }
    
    // Keep track of edge via which we entered the current simplex
    edge_prev_index = edge_curr_index;
    ++count;

    if (count == max_count){
      // std::cout << "\nError: hit maximum number of simplices to look for.";
      *final_simplex = -2;
      return M3DC1_FAILURE;

    }
  }
  return M3DC1_SUCCESS;
}

// VERIFY FIELDS
void send_dof(pMesh m, pMeshEnt e, pField f)
{
  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;

  double dof_data[FIXSIZEBUFF];
  getComponents(f, e, 0, dof_data);
  int n=countComponents(f);

    msg_size=sizeof(pMeshEnt)+n*sizeof(double);
    Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(Copies,remotes,rit)
    {
      int to = rit->first;
      msg_send = malloc(msg_size);
      s_ent = (pMeshEnt*)msg_send; 
      *s_ent = rit->second; 
      double *s_data = (double*)((char*)msg_send+sizeof(pMeshEnt));
      for (int pos=0; pos<n; ++pos)
        s_data[pos]=dof_data[pos];
      PCU_Comm_Write(to, (void*)msg_send, msg_size);
      free(msg_send);
    }

    if (m->isGhosted(e))
    {
      Copies g;
      m->getGhosts(e, g);
      APF_ITERATE(Copies, g, it)
      {
        int to = it->first;
        msg_send = malloc(msg_size);
        s_ent = (pMeshEnt*)msg_send; 
        *s_ent = it->second; 
        double *s_data = (double*)((char*)msg_send+sizeof(pMeshEnt)+sizeof(int));
        for (int pos=0; pos<n; ++pos)
          s_data[pos]=dof_data[pos];
        PCU_Comm_Write(to, (void*)msg_send, msg_size);
        free(msg_send);
      }
    } //if (m->isGhosted(e))
}

int receive_dof(pField f)
{
  void *msg_recv;
  int pid_from;
  size_t msg_size;
  pMeshEnt e;
  int res=0;

  while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
  {
    e = *((pMeshEnt*)msg_recv); 
    int n=countComponents(f);
    double* r_values = (double*)((char*)msg_recv+sizeof(pMeshEnt)); 
#ifdef DEBUG
    int num_data = (msg_size-sizeof(pMeshEnt))/sizeof(double);
    assert(n==num_data);
#endif 
    double dof_data[FIXSIZEBUFF];
    getComponents(f, e, 0, dof_data);  

    for (int i=0; i<n; ++i)
    {
      if (dof_data[i]!=r_values[i])
      {
        res=1;
        break;
      }
    }
  } // while
  return res;
}

void verify_field(pMesh m, pField f)
{
  PCU_Comm_Begin();
  pMeshIter it = m->begin(0);
  pMeshEnt e;
  while ((e = m->iterate(it)))
  {
    if (is_ent_original(m, e) && (m->isShared(e) || m->isGhosted(e)))
      send_dof(m, e, f);
  }
  m->end(it);
  PCU_Comm_Send();
  int mismatch = receive_dof(f); 

  int global_mismatch = PCU_Max_Int(mismatch);
  if (global_mismatch)
  {
    if (!PCU_Comm_Self())
      std::cout<<": failed\n";
  }
  else
  {
    if (!PCU_Comm_Self())
      std::cout<<": passed\n";
  }
}

void m3dc1_field_verify()
{
  if (!m3dc1_mesh::instance()->field_container)
    return;

  for (std::map<FieldID, m3dc1_field*>::iterator f_it=m3dc1_mesh::instance()->field_container->begin();
     f_it!=m3dc1_mesh::instance()->field_container->end();++f_it)
  {
    if (!pumi_rank()) std::cout<<"[M3D-C1 INFO] verifying field "<<getName(f_it->second->get_field());
    FieldID field_id = f_it->first;
    int isnan=0;  
    m3dc1_field_isnan(&field_id, &isnan);
    assert(isnan==0);
    verify_field(m3dc1_mesh::instance()->mesh, f_it->second->get_field()); 
  }
}

/* mesh entity functions */
//*******************************************************
void m3dc1_ent_getlocalid (int* /* in */ ent_dim, int* /* out */ ent_ids,
            int* /* in */ allocated_size, int* /* out */ num_ent)
//*******************************************************
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  if (*allocated_size<mesh->count(*ent_dim))
  {
    std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<" "<<__func__
               <<" failed: not enough array size for entity id's (allocated: "
               <<*allocated_size<<", needed: "<<mesh->count(*ent_dim)<<"\n";
    return;
  }

  apf::MeshEntity* e;
  apf::MeshIterator* it = mesh->begin(*ent_dim);
  int index=0;
  while ((e = mesh->iterate(it)))
    ent_ids[index++] = getMdsIndex(mesh, e);
  mesh->end(it);

  *num_ent = index;
}

//*******************************************************
int m3dc1_ent_getglobalid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ global_ent_id)
//*******************************************************
{
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);
  *global_ent_id = get_ent_globalid(m3dc1_mesh::instance()->mesh, e);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id, 
            int* /* out */ geom_class_dim, int* /* out */ geom_class_id)
//*******************************************************
{ 
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(ent);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(ent));
  *geom_class_dim = gmi_dim(m3dc1_model::instance()->model,gent);
  *geom_class_id = gmi_tag(m3dc1_model::instance()->model,gent);
  // if 3D mesh, need to return the classification on the original plane
  if ( m3dc1_mesh::instance()->mesh->getDimension() ==3 )
  {
    int numEntOrig[3];
    int numPlane = m3dc1_model::instance()->num_plane;
    memcpy( numEntOrig, m3dc1_model::instance()->numEntOrig, sizeof(numEntOrig));
    *geom_class_id-=1;
    switch (*geom_class_dim)
    {
      case 3: *geom_class_id%=numEntOrig[2]; break;
      case 2: *geom_class_id%=(numEntOrig[1]+numEntOrig[2]); break;
      case 1:  *geom_class_id%=(numEntOrig[0]+numEntOrig[1]); break;
      case 0: *geom_class_id%=(numEntOrig[0]);
    }
    *geom_class_id+=1;
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                      int* /* in */ adj_dim, int* /* out */ adj_ent, 
                      int* /* in */ adj_ent_allocated_size, int* /* out */ adj_ent_size)
//*******************************************************
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e = getMdsEntity(mesh, *ent_dim, *ent_id);
#ifdef DEBUG
  if (!e)
    print_ent_error(__func__, *ent_dim, *ent_id);
#endif
  assert(e);

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    mesh->getAdjacent(e,*adj_dim,adjacent);
    *adj_ent_size = adjacent.getSize();
    if (*adj_ent_allocated_size<*adj_ent_size)
    {
      std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<" "<<__func__
               <<" failed: not enough array size for adjacent entities (allocated: "
               <<*adj_ent_allocated_size<<", needed: "<<*adj_ent_size<<"\n";
      return M3DC1_FAILURE;
    }
    for (int i=0; i<*adj_ent_size; ++i)
      adj_ent[i] = getMdsIndex(mesh, adjacent[i]);
  }
  else if (*adj_dim<*ent_dim)
  {
    apf::Downward downward;
    *adj_ent_size = mesh->getDownward(e, *adj_dim, downward);
    if (*adj_ent_allocated_size<*adj_ent_size)
    {
      std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<" "<<__func__
               <<" failed: not enough array size for adjacent entities (allocated: "
               <<*adj_ent_allocated_size<<", needed: "<<*adj_ent_size<<"\n";
      return M3DC1_FAILURE;
    }
    for (int i=0; i<*adj_ent_size; ++i)
      adj_ent[i] = getMdsIndex(mesh, downward[i]);
    //adjust the order to work with m3dc1
    if (mesh->getDimension()==3 && *ent_dim==3 &&*adj_dim==0 &&adj_ent[0]>adj_ent[3])
    {
      int buff[3];
      memcpy(buff, adj_ent, 3*sizeof(int));
      memcpy(adj_ent, adj_ent+3, 3*sizeof(int));
      memcpy(adj_ent+3, buff, 3*sizeof(int));
    }
  }
  else // element's 2nd order adjacency
  {
    if (mesh->count(3))
      assert(*ent_dim==3);
    else
      assert(*ent_dim==2);

    apf::Adjacent adjacent;
    apf::getBridgeAdjacent(mesh, e, *ent_dim-1, *adj_dim, adjacent);
    *adj_ent_size = adjacent.getSize();
    if (*adj_ent_allocated_size<*adj_ent_size)
    {
      std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<" "<<__func__
               <<" failed: not enough array size for adjacent entities (allocated: "
               <<*adj_ent_allocated_size<<", needed: "<<*adj_ent_size<<"\n";
      return M3DC1_FAILURE;
    }
    for (int i=0; i<*adj_ent_size; ++i)
      adj_ent[i] = getMdsIndex(mesh, adjacent[i]);
  }
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                         int* /* in */ adj_dim, int* /* out */ num_adj_ent)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e || *adj_dim==*ent_dim)
    return M3DC1_FAILURE;

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
  }
  else if (*adj_dim<*ent_dim) 
  {
    apf::Downward downward;
    *num_adj_ent = m3dc1_mesh::instance()->mesh->getDownward(e, *adj_dim, downward);
  }
  return M3DC1_SUCCESS; 
}

//*******************************************************
void m3dc1_ent_getglobaladj (int* /* in */ ent_dim, 
                      int* /* in */ ent_ids, int* /* in */ num_ent,
                      int* /* in */ adj_dim,
                      int* /* out */ num_adj_ent, int* /* out */ adj_ent_pids, int* /* out */ adj_ent_gids, 
                      int* /* in */ adj_ent_allocated_size, int* /* out */ adj_ent_size)
//*******************************************************
{
  if (*adj_dim<*ent_dim)
  {
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 ERROR] "<<__func__
               <<" failed: adjacency dimension ("<<*adj_dim
               <<") should be greater than or equal to entity dimention ("
               <<").\n\t Use m3dc1_ent_getadj for downward adjacency\n";
    return;
  }
  else if (*adj_dim>*ent_dim)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" not supported yet for upward adjacency\n";

  }

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  int mesh_dim=mesh->count(3)?3:2;
  assert(*ent_dim==mesh_dim);

  apf::MeshEntity* e;
  std::vector<apf::MeshEntity*> ent_vec;
  std::vector<int> adj_gid_vec;
  std::vector<int> adj_pid_vec;
  std::vector<int> num_adj_vec;
  for (int i=0; i<*num_ent; ++i)
    ent_vec.push_back(getMdsEntity(mesh, *ent_dim, ent_ids[i]));

  *adj_ent_size = get_ent_global2ndadj(mesh, *ent_dim, *adj_dim, ent_vec, num_adj_vec, adj_pid_vec, adj_gid_vec);
  
  if (*adj_ent_allocated_size<*adj_ent_size)
  {
      std::cout<<"[M3D-C1 ERROR] p"<<PCU_Comm_Self()<<" "<<__func__
               <<" failed: not enough array size for adjacent entities (allocated: "
               <<*adj_ent_allocated_size<<", needed: "<<*adj_ent_size<<"\n";
      return;
  }

  memcpy(num_adj_ent, &(num_adj_vec[0]), num_adj_vec.size()*sizeof(int));
  memcpy(adj_ent_pids, &(adj_pid_vec[0]), adj_pid_vec.size()*sizeof(int));
  memcpy(adj_ent_gids, &(adj_gid_vec[0]), adj_gid_vec.size()*sizeof(int));
}

// allocated size of num_adj_ent should be greater than or equal to the element size
//*******************************************************
void m3dc1_ent_getnumglobaladj (int* /* in */ ent_dim, 
                      int* /* in */ ent_ids, int* /* in */ num_ent,
                      int* /* in */ adj_dim, int* /* out */ num_adj_ent)
//*******************************************************
{
  if (*adj_dim<*ent_dim)
  {
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 ERROR] "<<__func__
               <<" failed: adjacency dimension ("<<*adj_dim
               <<") should be greater than or equal to entity dimention ("
               <<").\n\t Use m3dc1_ent_getnumadj for downward adjacency\n";
    return;
  }
  else if (*adj_dim>*ent_dim)
  {
    if (!PCU_Comm_Self()) 
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" not supported yet for upward adjacency\n";
      
  }

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh; 
  std::vector<apf::MeshEntity*> ent_vec;
  std::vector<int> num_adj_vec;
  for (int i=0; i<*num_ent; ++i)
    ent_vec.push_back(getMdsEntity(mesh, *ent_dim, ent_ids[i]));
  get_ent_numglobaladj(mesh, *ent_dim, *adj_dim, ent_vec, num_adj_vec);
  memcpy(num_adj_ent, &(num_adj_vec[0]), *num_ent*sizeof(int));
}

//*******************************************************
int m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ owning_partid)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);
  *owning_partid = get_ent_ownpartid(m3dc1_mesh::instance()->mesh, e);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_ismine (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ ismine)
//*******************************************************
{
  *ent_id -= 1; //index change from Fortran to C
 
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  assert(e);

  if (is_ent_original(m,e)) 
     *ismine = 1;   // 
  else
     *ismine = 0; 
  return M3DC1_SUCCESS;
}

int m3dc1_ent_isghost(int* /* in */ ent_dim, int* /* in */ ent_id, int* isghost)
{
  pMeshEnt e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (pumi_ment_isGhost(e))
    *isghost=1;
  else
    *isghost=0;
  return M3DC1_SUCCESS;
}

//*******************************************************
// return volume, area or length of the entity
void m3dc1_ent_measure(int* /* in */ ent_dim, int* /* in */ ent_id,
                          double* /* out */ value)
//*******************************************************
{
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);
  *value = apf::measure(m3dc1_mesh::instance()->mesh, e);
}

// node-specific functions
//*******************************************************
void m3dc1_node_setfield (int* /* in */ node_id, int* /* in */ field_id,
                          double* /* in */ data, int* /* in */ size_data)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(e);
  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  if (*size_data != countComponents(f))
  {
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 ERROR] "<<__func__
               <<" failed: #data mismatch for field "<<getName(f)<<" (given: "
               <<*size_data <<", needed: "<<countComponents(f)<<"\n";
  }

  assert(*size_data == countComponents(f));
  apf::setComponents(f, e, 0, data);
}

//*******************************************************
void m3dc1_node_getfield (int* /* in */ node_id, int* /* in */ field_id,
                          double* /* inout */ data, int* /* in */ allocated_data)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(e);
  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  if (*allocated_data < countComponents(f))
  {
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 ERROR] "<<__func__
               <<" failed: not enough array size for field "<<getName(f) <<" (allocated: "
               <<*allocated_data <<", needed: "<<countComponents(f)<<"\n";
      return;
  }
  apf::getComponents(f, e, 0, data);
}

//*******************************************************
int m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(e);
  apf::Vector3 xyz;
  m3dc1_mesh::instance()->mesh->getPoint(e, 0, xyz);
  for (int i=0; i<3; ++i)
    coord[i] = xyz[i]; 
  return M3DC1_SUCCESS;
}

//*******************************************************
void get_gv_edges(gmi_ent* gvertex, std::vector<gmi_ent*>& gedges)
//*******************************************************
{
  gedges.clear();
  gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, gvertex, 1);
  for (int i=0; i<gv_edges->n; ++i)
    gedges.push_back(gv_edges->e[i]);
  gmi_free_set(gv_edges);
  assert(gedges.size()>=1);
}

//*******************************************************
int m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyzt)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);
  xyzt[2]=0.0;

  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if (gType !=  1 && gType !=  0)
  {
    xyzt[0] = xyzt[1] = 0.0;
    return M3DC1_SUCCESS;
  }

  if (m3dc1_model::instance()->modelType == 2)
  {
    apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
    if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
    {
      double norm_curv[3];
      m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
      xyzt[0] = norm_curv[0]; 
      xyzt[1] = norm_curv[1];
      return M3DC1_SUCCESS;
    }
    else 
    {
      std::cout << "Warning: In case of .dmg models, normal vectors are only provided at \n";
      std::cout << "predefined boundaries while creating mesh, make sure this function is \n";
      std::cout << "called at the correct model edges at correct boundary loop\n";

      // Dummy to return normals so code doesn't crash in case 
      xyzt[0] = 1.0;
      xyzt[1] = 0.0;
      return M3DC1_SUCCESS;
    }
  }

  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    xyzt[0] = norm_curv[0]; 
    xyzt[1] = norm_curv[1];
    return M3DC1_SUCCESS;
  }
  else
  { // if norm/curv is not attached, evaluate
    apf::Vector3 param(0,0,0);
    m3dc1_mesh::instance()->mesh->getParam(vt,param);
    // geo node avage on the connected edges
    if (gType == 0) // node is on the 
    {
      apf::Vector3 vcd_t;
      double vcd[3];
      m3dc1_mesh::instance()->mesh->getPoint(vt, 0, vcd_t);
      for (int i=0; i<3; i++) 
        vcd[i]=vcd_t[i];
      std::vector<gmi_ent*> gEdges;
      get_gv_edges(gent, gEdges);
      int numEdgePlane=0;
      double normalvec[3]={0.,0.,0.};
      xyzt[0]=xyzt[1]=xyzt[2]=0;
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: #adjEdge of gVertex="
		 <<gEdges.size()<<" (it should be minimum 2) \n";
      assert(gEdges.size()>=2);

      for (int i=0;i<gEdges.size();i++)
      {
        gmi_ent* pe = gEdges.at(i);
        double cd[3]={0,0,0};
        double paraRange[2];
        gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
        M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
        if (!pn) continue;
        numEdgePlane++;
        M3DC1::evalCoord(paraRange[0],cd, pn);
        if (checkSamePoint2D(vcd,cd))
          M3DC1::evalNormalVector(pn[0],pn[1], paraRange[0], normalvec);
        else
          evalNormalVector(pn[0],pn[1], paraRange[1], normalvec);
        xyzt[0]+=normalvec[0];
        xyzt[1]+=normalvec[1];
        xyzt[2]+=normalvec[2];
      }
      if (numEdgePlane!=2) 
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: numEdgePlane="
                 <<numEdgePlane<<" (it should be 2) \n";
      assert(numEdgePlane==2);

      double arclen=sqrt(xyzt[0]*xyzt[0]+xyzt[1]*xyzt[1]+xyzt[2]*xyzt[2]);
      assert(arclen>0);
      xyzt[0]=xyzt[0]/arclen;
      xyzt[1]=xyzt[1]/arclen;
      xyzt[2]=xyzt[2]/arclen;
    }
    else
    {
      apf::Vector3 param(0,0,0);
      m3dc1_mesh::instance()->mesh->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalNormalVector(pn[0],pn[1], param[0], xyzt);
    }
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);

  if (m3dc1_model::instance()->modelType == 2)
  {
    apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
    if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
    {
      double norm_curv[3];
      m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
      *curv = norm_curv[2];
      return M3DC1_SUCCESS;
    }
    else 
    {
      std::cout << "Warning: In case of .dmg models, curvatures are only provided at \n";
      std::cout << "predefined boundaries while creating mesh, make sure this function is \n";
      std::cout << "called at the correct model edges at correct boundary loop\n";

      // Dummy to return normals so code doesn't crash in case 
      *curv = 0.0;
      return M3DC1_SUCCESS;
    }
  }
  

  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    *curv = norm_curv[2]; 
    return M3DC1_SUCCESS;
  }

  *curv=0.0;
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if (gType==0)
  {
    apf::Vector3 vcd_t;
    double vcd[3];
    m3dc1_mesh::instance()->mesh->getPoint(vt, 0, vcd_t);
    for (int i=0; i<3; i++)
      vcd[i]=vcd_t[i];
    std::vector<gmi_ent*> gEdges;
    get_gv_edges(gent, gEdges);
    int numEdgesPlane=0;
    double curv_tmp;
    for (int i=0;i<gEdges.size();i++)
    {
      gmi_ent* pe = gEdges.at(i);
      double cd[3]={0,0,0};
      double paraRange[2];
      gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
      if (!pn) continue;
      numEdgesPlane++;
      evalCoord(paraRange[0], cd, pn);
      if (checkSamePoint2D(vcd,cd))
      {
        evalCurvature(pn[0],pn[1], paraRange[0], &curv_tmp); 
      }
      else
      {
        evalCurvature(pn[0],pn[1], paraRange[1], &curv_tmp);
      }
      *curv+=curv_tmp;
    }
    assert(numEdgesPlane==2);
    *curv/=numEdgesPlane;
  }
  else if (gType==1)
  {
      apf::Vector3 param(0,0,0);
      m3dc1_mesh::instance()->mesh->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalCurvature(pn[0],pn[1], param[0], curv);
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  *on_geom_bdry=(gType==0||gType==1);
  return M3DC1_SUCCESS;
}

//=========================================================================
void write_node(apf::Mesh2* m, const char* filename, int start_index)
{
  char node_filename[256];
  sprintf(node_filename,"%s-%d", filename, PCU_Comm_Self());
  FILE * np =fopen(node_filename, "w");

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    apf::Vector3 xyz;
    m->getPoint(e, 0, xyz);
    fprintf(np, "global ID: %d\t%lf\t%lf\t%lf\n", get_ent_globalid(m,e)+start_index, xyz[0],xyz[1],xyz[2]);
  } // while
  m->end(it);
  fclose(np);
}

//*******************************************************
int m3dc1_node_write (const char* filename, int* start_index)
//*******************************************************
{
  write_node(m3dc1_mesh::instance()->mesh, filename, *start_index);
  return M3DC1_SUCCESS;
}

// mesh region functions
//*******************************************************
int m3dc1_region_getoriginalface( int * /* in */ elm, int * /* out */ fac)
//*******************************************************
{
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->mesh, 3, *elm);
  apf::Downward downward;
  int num_adj_ent = m3dc1_mesh::instance()->mesh->getDownward(ent, 2, downward);
  assert(num_adj_ent==5);
  int triFace[2];
  int counter=0;

  for (int i=0; i<num_adj_ent; i++)
  {
    int num_adj_ent;
    apf::Downward downward2;
    int num_edge= m3dc1_mesh::instance()->mesh->getDownward(downward[i], 1, downward2);
    if (num_edge==3) triFace[counter++]= getMdsIndex(m3dc1_mesh::instance()->mesh,downward[i]);    
  }
  assert(counter==2);
  *fac = std::min(triFace[0],triFace[1]);
  return M3DC1_SUCCESS;
}

/** field manangement */
int fieldIdMax=0;
//*******************************************************
int m3dc1_field_getnewid ( FieldID* /*out*/field_id )
//*******************************************************
{
  *field_id = fieldIdMax+1;
  return M3DC1_SUCCESS;
}

// *scalar_type is either M3DC1_REAL or M3DC1_COMPLEX
int m3dc1_field_create (FieldID* /*in*/ field_id, const char* /* in */ field_name, int* /*in*/ num_values, 
int* /*in*/ scalar_type, int* /*in*/ num_dofs_per_value)
{
  if (!m3dc1_mesh::instance()->field_container)
    m3dc1_mesh::instance()->field_container=new std::map<FieldID, m3dc1_field*>;

  // shape evaluation will be performed outside the APF
  // only need to tell APF all dofs are attached to mesh vertex
  int components = (*num_values)*(*scalar_type+1)*(*num_dofs_per_value);
  apf::Field* f = createPackedField(m3dc1_mesh::instance()->mesh, field_name, components);
  m3dc1_mesh::instance()->field_container->insert(std::map<FieldID, m3dc1_field*>::value_type(*field_id, 
      new m3dc1_field(*field_id, f, *num_values, *scalar_type, *num_dofs_per_value)));
  apf::freeze(f); // switch dof data from tag to array

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": field "<<*field_id<<", #values "
             <<*num_values<<", #dofs "<<countComponents(f)<<", name "<<field_name<<"\n";

  if (*field_id>fieldIdMax) fieldIdMax=*field_id;
  double val[2]={0,0};
  m3dc1_field_assign(field_id, val, scalar_type);
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_field_mark4tx (FieldID* /*in*/ field_id)
//*******************************************************
{
  assert (m3dc1_mesh::instance()->field_container &&
          m3dc1_mesh::instance()->field_container->count(*field_id));
  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
#ifdef DEBUG
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": field "<<*field_id<<", name "<<getName(f)
             <<" marked for solution transfer\n";
#endif
  m3dc1_field* mf = (*m3dc1_mesh::instance()->field_container)[*field_id];
  mf->mark_for_solutiontransfer();
}

//*******************************************************
int m3dc1_field_delete (FieldID* /*in*/ field_id)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container)
    return M3DC1_FAILURE;
  if (!m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
#ifdef DEBUG
  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": field "<<*field_id<<", name "<<getName(f)<<"\n";
#endif

  m3dc1_mesh::instance()->mesh->removeField(f);
  destroyField(f);

  // remove f from field container
  delete (*m3dc1_mesh::instance()->field_container)[*field_id];
  m3dc1_mesh::instance()->field_container->erase(*field_id);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getinfo(FieldID* /*in*/ field_id, 
                        char* /* out*/ field_name, int* num_values, 
                        int* scalar_type, int* total_num_dof)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf ->get_field();
  strcpy(field_name, getName(f));
  *num_values = mf -> get_num_value();
  *scalar_type = mf ->get_value_type();
  *total_num_dof = countComponents(f);
  if (*scalar_type) *total_num_dof/=2;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_exist(FieldID* field_id, int * exist)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    *exist = 0;
  else
    *exist = 1;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_sync (FieldID* /* in */ field_id)
//*******************************************************
{
  if (PCU_Comm_Peers()==0) return 0;

#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  synchronize_field((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());

#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

// send non-owned copies' dof to owner copy and add them up
//*******************************************************
void accumulate_field(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;       

  int own_partid, n = countComponents(f);
  double* dof_data = new double[n];
  double* sender_data = new double[n];
  apf::MeshEntity* own_e;
  apf::MeshEntity* r;
  std::map<apf::MeshEntity*, std::vector<double> > save_map;

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid=get_ent_ownpartid(m, e);
    if (own_partid==PCU_Comm_Self() || pumi_ment_isGhost(e)) continue;
    assert(m->isShared(e));

    own_e = get_ent_owncopy(m, e);

    getComponents(f, e, 0, &(dof_data[0]));
      
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(dof_data[0]),n*sizeof(double));
  }
  m->end(it);

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while (! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      for (int i = 0; i < n; ++i)
        save_map[r].push_back(sender_data[i]);      
    }

  for (std::map<apf::MeshEntity*, std::vector<double> >::iterator mit=save_map.begin(); mit!=save_map.end(); ++mit)
  {
    e = mit->first;
    getComponents(f, e, 0, dof_data);
    int num_data = mit->second.size()/n;
    for (int i=0; i<num_data;++i)
    {
      for (int j=0; j<n; ++j)
        dof_data[j] += mit->second[i*n+j];
    }
    setComponents(f, e, 0, dof_data);
  } 
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_field_sum (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  accumulate_field((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());
  synchronize_field((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());

#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_sumsq (FieldID* /* in */ field_id, double* /* out */ sum)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();

#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  *sum=0.;
  int num_dof = countComponents(f);

  double* dof_data= new double[num_dof];
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    if (!is_ent_original(m3dc1_mesh::instance()->mesh,e)) continue;
    getComponents(f, e, 0, dof_data);
    for (int i=0; i<num_dof; ++i)
      *sum+=dof_data[i]*dof_data[i];
  }
  m3dc1_mesh::instance()->mesh->end(it);
  delete [] dof_data;
  return M3DC1_SUCCESS;
}

/** field dof functions */
//*******************************************************
int m3dc1_field_getlocaldofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  *start_dof_id=0;
  *end_dof_id_plus_one=num_dof*m3dc1_mesh::instance()->mesh->count(0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getowndofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0];
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  
  int start_id = num_own_ent;
  PCU_Exscan_Ints(&start_id,1);

  *start_dof_id=start_id*num_dof;
  *end_dof_id_plus_one=*start_dof_id+num_own_ent*num_dof;
  return M3DC1_SUCCESS;
}
 
//******************************************************* 
int m3dc1_field_getglobaldofid ( FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  assert(mf->get_num_value()*mf->get_dof_per_value()==num_dof);  

  *start_dof_id=0;
  *end_dof_id_plus_one=*start_dof_id+num_dof*m3dc1_mesh::instance()->num_global_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumlocaldof (FieldID* field_id, int* /* out */ num_local_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_local_dof = (m3dc1_mesh::instance()->mesh->count(0))*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumowndof (FieldID* field_id, int* /* out */ num_own_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_own_dof = (m3dc1_mesh::instance()->num_own_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumglobaldof (FieldID* field_id, int* /* out */ num_global_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_global_dof = (m3dc1_mesh::instance()->num_global_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumghostdof (FieldID* field_id, int* /* out */ num_ghost_dof)
//*******************************************************
{
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int num_ent = m3dc1_mesh::instance()->mesh->count(0)-m3dc1_mesh::instance()->num_local_ent[0];
  *num_ghost_dof = num_ent*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_field_getdataptr (FieldID* field_id, double** pts)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  if (!isFrozen(f)) freeze(f);
  *pts=getArrayData(f);
  return M3DC1_SUCCESS;
}

// add field2 to field1
//*******************************************************
int m3dc1_field_add(FieldID* /*inout*/ field_id1, FieldID* /*in*/ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value(); 
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;
  
  int dofPerEntDummy[2];
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for (int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)+=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  } 
#ifdef DEBUG
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_mult(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;

  double dofs[FIXSIZEBUFF], dofsNew[FIXSIZEBUFF];
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  int value_type = mf->get_value_type();
  assert(dofPerEnt<=sizeof(dofs)/sizeof(double)*(1+value_type));
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    if (*scalar_type==0)
    {
      for (int i=0; i<dofPerEnt*(1+value_type); i++)
dofsNew[i]=*fac*dofs[i];
    }
    else
    {
      for (int i=0; i<dofPerEnt; i++)
      {
dofsNew[2*i]=fac[0]*dofs[2*i]-fac[1]*dofs[2*i+1];
dofsNew[2*i+1]=fac[0]*dofs[2*i+1]+fac[1]*dofs[2*i];
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofsNew[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_assign(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;

  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()), fac[0]);
  if (*scalar_type)
    for (int i=0; i<dofPerEnt; i++)
      dofs.at(2*i+1)=fac[1];
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
  }
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_copy(FieldID* /* out */ field_id1, FieldID* /* in */ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  apf::Field* f1 =  mf1->get_field();
  apf::Field* f2 =  mf2->get_field();
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value();
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;

  int dofPerEntDummy[2];
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for (int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_retrieve (FieldID* /* in */ field_id, double * /*out*/ data, int * /* in */size)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  //m3dc1_field_print(field_id);
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(data, pts, *size*(1+value_type)*sizeof(double));
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_set (FieldID* /* in */ field_id, double * /*in*/ data, int * /* in */size)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(pts, data, *size*(1+value_type)*sizeof(double));
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_insert(FieldID* /* in */ field_id, int /* in */ * local_dof, 
         int * /* in */ size, double* /* in */ values, int * type, int * op)
//*******************************************************
{
  int num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
#ifdef DEBUG
  int num_local_dof;
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);
  assert(*local_dof<num_local_dof);
  if (!value_type) assert(!(*type)); // can not insert complex value to real vector
  // for (int i=0; i<*size*(1+(*type)); i++)
    // seol #1: crash on SCOREC linux clusters with 3d/3p
    // assert(values[i]==values[i]);
#endif
  std::vector<double> values_convert(*size*(1+value_type),0);
  if (!(*type)&&value_type) // real into complex
  {
    for (int i=0; i<*size; i++)
      values_convert.at(2*i)=values[i];
  }
  else
  {
    for (int i=0; i<*size*(1+value_type); i++)
      values_convert.at(i)=values[i];
  }
  double * dataptr;
  int ibegin=*local_dof*(1+value_type);
  m3dc1_field_getdataptr(field_id, &dataptr);
  if (*op==0) // set value
   for (int i=0; i<*size*(1+value_type); i++)
     dataptr[ibegin+i]=values_convert.at(i);
  else
    for (int i=0; i<*size*(1+value_type); i++)
      dataptr[ibegin+i]+=values_convert[i];
  return M3DC1_SUCCESS;
}

#define FIELDVALUELIMIT 1e100
bool value_is_nan(double val)
{
  return val!=val ||fabs(val) >FIELDVALUELIMIT;
}

//*******************************************************
int m3dc1_field_isnan(FieldID* /* in */ field_id, int * isnan)
//*******************************************************
{
  *isnan=0;
  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;

  int dofPerEnt;
  double dofs[FIXSIZEBUFF];
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    for (int i=0; i<dofPerEnt; i++)
      if (value_is_nan(dofs[i])) 
        *isnan=1;
  }
  return M3DC1_SUCCESS;
}

//=========================================================================
void write_vector(apf::Mesh2* m, m3dc1_field* mf, const char* filename, int start_index)
{
  std::string in(filename);
//  size_t p = in.rfind('.');
//  std::string base = in.substr(0,p);
//  std::string ext = in.substr(p+1);
  std::stringstream s;
  s <<in << "-"<<PCU_Comm_Self();
  std::string partFile = s.str();
  FILE * fp =fopen(partFile.c_str(), "w");

  apf::MeshEntity* e;
  apf::Field* f = mf->get_field();
  int num_dof=countComponents(f);
  double* dof_data = new double[num_dof];
  fprintf(fp, "name %s, #value %d, #dof/value %d, scalar type %d\n", getName(f), mf->get_num_value(), mf->get_dof_per_value(), mf->get_value_type());
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    getComponents(f, e, 0, dof_data);
    int ndof=0;
    for (int i=0; i<num_dof; ++i)
    {
      if (!m3dc1_double_isequal(dof_data[i], 0.0))
        ++ndof;
    }
    fprintf(fp, "\nglobal ID %d, local ID %d, #dof %d\n", get_ent_globalid(m,e)+start_index, getMdsIndex(m, e), ndof);
    for (int i=0; i<num_dof; ++i)
    {
      if (!m3dc1_double_isequal(dof_data[i], 0.0))
        fprintf(fp, "dof %d: %lf\n", i, dof_data[i]);
    }
  } // while
  m->end(it);
  fclose(fp);
  delete [] dof_data;
}

//*******************************************************
int m3dc1_field_write (FieldID* field_id, const char* filename, int* start_index)
//*******************************************************
{
  if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<"(field id "<<*field_id<<", file "<<filename<<")\n";
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  write_vector(m, (*(m3dc1_mesh::instance()->field_container))[*field_id], filename, *start_index);
  return M3DC1_SUCCESS;
}

void m3dc1_field_export(int* field_id, int* ts)
{ 
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  assert(mf);

  apf::Field* f;
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int num_dof;
  
  int file_id=*field_id;
    char fname[64];
    sprintf(fname, "fld-%d-ts%d-%d",file_id, *ts, PCU_Comm_Self());
    
    FILE * fp =fopen(fname, "w");
    
    f = mf->get_field();
    
    fprintf(fp, "%s\n", getName(f));
    fprintf(fp, "%d %d %d %d\n", // m3dc1_mesh::instance()->field_container->size(),
            mf->get_id(), mf->get_num_value(), mf->get_dof_per_value(), mf->get_value_type());
    
    if (!PCU_Comm_Self()) std::cout<<__func__<<" file "<<fname<<": field "<<mf->get_id()<<"\n";
    
    num_dof=countComponents(f);
    double* dof_data = new double[num_dof];
    
    it = m->begin(0);
    while ((e = m->iterate(it)))
    { 
      getComponents(f, e, 0, dof_data);
      int ndof=0;
      for (int i=0; i<num_dof; ++i)
      { 
        if (!m3dc1_double_isequal(dof_data[i], 0.0))
        ++ndof;
      }
      fprintf(fp, "%d ", ndof);
      for (int i=0; i<num_dof; ++i)
      { 
        if (!m3dc1_double_isequal(dof_data[i], 0.0))
          fprintf(fp, "%d %lf ", i, dof_data[i]);
      }
      fprintf(fp, "\n");
    } // while
    m->end(it);
    delete [] dof_data;
    fclose(fp);
}


void m3dc1_field_exportall()
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  m3dc1_field* mf;
  apf::Field* f;
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int num_dof;

  int file_id=0;
  std::map<FieldID, m3dc1_field*>::iterator fit=m3dc1_mesh::instance()->field_container->begin();
  for (; fit!=m3dc1_mesh::instance()->field_container->end(); ++fit)
  {
    char fname[64];
    sprintf(fname, "fld-%d-%d",file_id++, PCU_Comm_Self());

    FILE * fp =fopen(fname, "w");

    mf = fit->second;
    f = mf->get_field();

    fprintf(fp, "%s\n", getName(f));
    fprintf(fp, "%d %d %d %d\n", // m3dc1_mesh::instance()->field_container->size(),
            mf->get_id(), mf->get_num_value(), mf->get_dof_per_value(), mf->get_value_type());

    if (!PCU_Comm_Self()) std::cout<<__func__<<" file "<<fname<<": field "<<mf->get_id()<<"\n";

    num_dof=countComponents(f);
    double* dof_data = new double[num_dof];

    it = m->begin(0);
    while ((e = m->iterate(it)))
    {
      getComponents(f, e, 0, dof_data);
      int ndof=0;
      for (int i=0; i<num_dof; ++i)
      {
        if (!m3dc1_double_isequal(dof_data[i], 0.0))
        ++ndof;
      }
      fprintf(fp, "%d ", ndof);
      for (int i=0; i<num_dof; ++i)
      {
        if (!m3dc1_double_isequal(dof_data[i], 0.0))
          fprintf(fp, "%d %lf ", i, dof_data[i]);
      }
      fprintf(fp, "\n");
    } // while
    m->end(it);
    delete [] dof_data;
    fclose(fp);
  }
}

void m3dc1_field_import(int *field_id, int* ts)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  m3dc1_field* mf;
  apf::Field* f;

  int file_id=*field_id;
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int num_dof, non_zero_dof, k;
  int fid, num_val, num_dof_per_val, val_type;

  double zero_val[2]={0,0};

  char fname[64];
  sprintf(fname, "fld-%d-ts%d-%d", file_id, *ts, PCU_Comm_Self());
  FILE * fp =fopen(fname, "r");
  if (!fp)
  {
    if (!PCU_Comm_Self()) std::cout<<"field file "<<fname<<" doesn't exist\n";
    return;
  }

  char field_name[50];
  fscanf(fp, "%s\n", field_name);
  fscanf(fp, "%d %d %d %d\n", &fid, &num_val, &num_dof_per_val, &val_type);

  assert(fid==*field_id);

  f = m->findField(field_name);
  if (!f)
  {
    m3dc1_field_create (&fid, field_name, &num_val, &val_type, &num_dof_per_val);
    f = (*m3dc1_mesh::instance()->field_container)[fid]->get_field();
  }

    num_dof=countComponents(f);
    assert(num_dof==num_val*num_dof_per_val*(val_type+1));
    double* dof_data = new double[num_dof];

    m3dc1_field_assign(&fid, zero_val, &val_type);

    it = m->begin(0);
    while ((e = m->iterate(it)))
    {
      for (int i=0; i<num_dof; ++i)
        dof_data[i]=0.0;
      fscanf(fp, "%d\n", &non_zero_dof);
      for (int i=0; i<non_zero_dof; ++i)
      {
        fscanf(fp, "%d", &k);
        fscanf(fp, "%lf", &(dof_data[k]));
      }
      setComponents(f, e, 0, dof_data);
    } // while
    m->end(it);
    delete [] dof_data;
    fclose(fp);
}

void m3dc1_field_importall()
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  m3dc1_field* mf;
  apf::Field* f;

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int num_dof, non_zero_dof, k;
  int fid, num_val, num_dof_per_val, val_type;

  int file_id=0;
  double zero_val[2]={0,0};

  while (1)
  {
    char fname[64];
    sprintf(fname, "fld-%d-%d",file_id, PCU_Comm_Self());
    FILE * fp =fopen(fname, "r");
    if (!fp) 
      break;

    char field_name[50];
    fscanf(fp, "%s\n", field_name);
    fscanf(fp, "%d %d %d %d\n", &fid, &num_val, &num_dof_per_val, &val_type);

    if (!PCU_Comm_Self()) std::cout<<__func__<<": importing field "<<fname<<" (ID "<<fid<<")\n";

    f = m->findField(field_name);
    if (!f)
    {
      //if (!PCU_Comm_Self()) std::cout<<"field "<<field_name<<" doesn't exists\n";
      m3dc1_field_create (&fid, field_name, &num_val, &val_type, &num_dof_per_val);
      f = (*m3dc1_mesh::instance()->field_container)[fid]->get_field();
    }
 
    num_dof=countComponents(f);
    assert(num_dof==num_val*num_dof_per_val*(val_type+1));
    double* dof_data = new double[num_dof];
  
    m3dc1_field_assign(&fid, zero_val, &val_type);
     
    it = m->begin(0);
    while ((e = m->iterate(it)))
    {
      for (int i=0; i<num_dof; ++i) dof_data[i]=0.0;

      fscanf(fp, "%d\n", &non_zero_dof);
      for (int i=0; i<non_zero_dof; ++i)
      {
        fscanf(fp, "%d", &k);
        fscanf(fp, "%lf", &(dof_data[k]));
      }
      setComponents(f, e, 0, dof_data);
    } // while
    m->end(it);
    delete [] dof_data;
    fclose(fp);
    ++file_id;
  }
}


//*******************************************************
int m3dc1_field_print(FieldID* field_id)
//*******************************************************
{
  apf::Mesh2*  m = m3dc1_mesh::instance()->mesh;
  apf::Field*  f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  
  if (!f)
  {
    if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<" failed as field "<<*field_id<<" not found\n";
    return M3DC1_FAILURE;
  }

  double* field_data =getArrayData(f);

  if (!m->findField("node global id field"))
  {
    if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<" failed as global node id not found\n";
    return M3DC1_FAILURE;
  }

  apf::MeshEntity* e;
  int num_dof=countComponents(f);
  double* dof_data = new double[num_dof];
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    switch (num_dof)
    {
      case 1: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<"]\n";
        break;}
      case 2: {     
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<"]\n";
        break;}
      case 3: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<"]\n";
        break;}
    case 4: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<"]\n";
 
        break; }
      case 6: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<"]\n";
        break; }
      case 8: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<"]\n";
        break; }
      case 12: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<"]\n";
        break; }
      case 18: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<"]\n";
        break; }
      case 24: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<", "<<dof_data[18]
		     <<", "<<dof_data[19]
		     <<", "<<dof_data[20]
		     <<", "<<dof_data[21]
		     <<", "<<dof_data[22]
		     <<", "<<dof_data[23]
		     <<"]\n";
        break; }
      default: if (!PCU_Comm_Self()) std::cout<<__func__<<" failed for field "
               <<getName(f)<<": does support "<<num_dof<<" dofs\n";
               break;
    } // switch
  } // while
  m->end(it);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_compare(FieldID* field_id_1, FieldID* field_id_2)
//*******************************************************
{
  apf::Field* f_1 = (*(m3dc1_mesh::instance()->field_container))[*field_id_1]->get_field();
  double* field_data_1 =getArrayData(f_1);

  apf::Field* f_2 = (*(m3dc1_mesh::instance()->field_container))[*field_id_2]->get_field();
  double* field_data_2 =getArrayData(f_2);

  int num_dof_1=countComponents(f_1);
  int num_dof_2=countComponents(f_2);
  if (num_dof_1!=num_dof_2) 
  {
    if (!PCU_Comm_Self()) 
      cout<<"[M3D-C1 INFO] "<<__func__<<": #dof mismatch "<<getName(f_1)
          <<"- "<<num_dof_1<<", "<<getName(f_2)<<"- "<<num_dof_2<<"\n";
    return M3DC1_FAILURE;
  }
  int ierr = M3DC1_SUCCESS;
  for (int i=0; i<num_dof_1*m3dc1_mesh::instance()->mesh->count(0); ++i)
  {  if (!m3dc1_double_isequal(field_data_1[i], field_data_2[i])) 
    {
     cout<<"[M3D-C1 ERROR] "<<__func__<<": "<<getName(f_1)<<"["<<i<<"]="<<field_data_1[i]
          <<", "<<getName(f_2)<<"["<<i<<"]="<<field_data_2[i]<<"\n";
      ierr=M3DC1_FAILURE;
      break;
    }
  }
  int global_ierr;
  MPI_Allreduce(&ierr, &global_ierr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  if (global_ierr==M3DC1_FAILURE)
  {
    if (!PCU_Comm_Self())
      cout<<"[M3D-C1 INFO] "<<__func__<<": dof value mismatch of fields "<<getName(f_1)
          <<" and "<<getName(f_2)<<"\n";
    
    return M3DC1_FAILURE;
  }
  else
    return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                       int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);

  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;
  *start_dof_id = *ent_id*dof_per_node;
  *end_dof_id_plus_one = *start_dof_id +dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
         int* /* out */ start_global_dof_id, int* /* out */ end_global_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);

  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;

  int global_id = get_ent_globalid(m3dc1_mesh::instance()->mesh, e);
  *start_global_dof_id = global_id*dof_per_node;
  *end_global_dof_id_plus_one =*start_global_dof_id + dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getghostdofid (FieldID* field_id, 
    int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || 
      !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  *start_dof_id=num_dof*m3dc1_mesh::instance()->num_local_ent[0];
  *end_dof_id_plus_one=num_dof*m3dc1_mesh::instance()->mesh->count(0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getnumdof (int* /* in */ ent_dim, int* /* in */ ent_id, 
         FieldID* field_id, int* /* out */ num_dof)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_dof =  mf->get_num_value() * mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

#ifdef DEBUG
  assert(countComponents(f)==*num_dof*(1+mf->get_value_type()));
  for (int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
#endif
  setComponents(f, e, 0, dof_data);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(e);

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  getComponents(f, e, 0, dof_data);
  *num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
#ifdef DEBUG
  for (int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
  int start_dof_id,end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid(ent_dim, ent_id,field_id, &start_dof_id, &end_dof_id_plus_one);
  double* data;
  m3dc1_field_getdataptr(field_id, &data);
  int start=start_dof_id*(1+mf->get_value_type());
  for ( int i=0; i< *num_dof; i++)
    assert(data[start++]==dof_data[i]);
#endif
  return M3DC1_SUCCESS;
}

#ifdef M3DC1_PETSC
/** matrix and solver functions */
std::map<int, int> matHit;
int getMatHit(int id) { return matHit[id];};
void addMatHit(int id) { matHit[id]++; }

//*******************************************************
int m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID *field_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }

#ifdef DEBUG
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": matrix "<<*matrix_id<<", field "<<*field_id<<"\n";
#endif 

  if (*matrix_type==M3DC1_MULTIPLY) // matrix for multiplication
  {
    matrix_mult* new_mat = new matrix_mult(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }
  else 
  {
    matrix_solve* new_mat= new matrix_solve(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_assemble(int* matrix_id) 
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<"\n";
  if (!mat) 
  {   
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif
  mat->assemble();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_delete(int* matrix_id)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<"\n";
  if (!mat) 
  {   
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif
  typedef std::map<int, m3dc1_matrix*> matrix_container_map;
  m3dc1_solver::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  delete mat;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_matrix_reset(int* matrix_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<"\n";
  if (!mat)
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return;
  }
#endif
  mat->reset_values();
}


//*******************************************************
void m3dc1_matrix_update(int* matrix_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<"\n";
  if (!mat)
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return;
  }
#endif
  mat->update_values();
}


//*******************************************************
int m3dc1_matrix_insert(int* matrix_id, int* row, 
         int* col, int* scalar_type, double* val)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }

  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);

  int ent_id = *row/total_num_dof;
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, ent_id);
  assert(e);
  assert(!m3dc1_mesh::instance()->mesh->isGhost(e));
#endif

  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, INSERT_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, INSERT_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_add (int* matrix_id, int* row, int* col, 
                      int* scalar_type, double* val) //globalinsertval_
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat) 
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }

  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);

  int ent_id = *row/total_num_dof;
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, ent_id);
  assert(e);
  assert(!m3dc1_mesh::instance()->mesh->isGhost(e));
#endif

  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, ADD_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, ADD_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_setbc(int* matrix_id, int* row)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat) 
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif

  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);

#ifdef DEBUG
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, inode);
  assert(e);
  assert(!m3dc1_mesh::instance()->mesh->isGhost(e));

  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  int row_g = start_global_dof_id+*row%total_num_dof;
  (dynamic_cast<matrix_solve*>(mat))->set_bc(row_g);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_setlaplacebc(int * matrix_id, int *row,
         int * numVals, int *columns, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat) 
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }

  if (mat->get_type()!=M3DC1_SOLVE)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
#endif
  std::vector <int> columns_g(*numVals);
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);

#ifdef DEBUG
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, inode);
  assert(e);
  assert(!m3dc1_mesh::instance()->mesh->isGhost(e));
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif

  int row_g = start_global_dof_id+*row%total_num_dof;
  for (int i=0; i<*numVals; i++)
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;

  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
  return M3DC1_SUCCESS;
}

int m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol) //solveSysEqu_
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<*rhs_sol<<"\n";
  if (!mat) 
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }

  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
#endif

  (dynamic_cast<matrix_solve*>(mat))->solve(*rhs_sol);
  addMatHit(*matrix_id);
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_matrix_solve_with_guess(int* matrix_id, FieldID* rhs_sol, FieldID* xVec_guess)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<*rhs_sol<<"\n";
  if (!mat) 
  {  
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return;
  }

  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return;
  }
#endif

  (dynamic_cast<matrix_solve*>(mat))->solve_with_guess(*rhs_sol, *xVec_guess);
  addMatHit(*matrix_id);
}

//*******************************************************
int m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, 
         FieldID* outputvecid) 
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", in-field "<<*inputvecid<<", out-field "<<*outputvecid<<"\n";
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }

  if (mat->get_type()!=M3DC1_MULTIPLY)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for solving (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
#endif

  (dynamic_cast<matrix_mult*>(mat))->multiply(*inputvecid, *outputvecid);
  addMatHit(*matrix_id);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_getnumiter(int* matrix_id, int * iter_num)
//*******************************************************
{ 
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif
  *iter_num = (int)(dynamic_cast<matrix_solve*> (mat)->its);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_insertblock(int* matrix_id, int * ielm, 
          int* rowIdx, int * columnIdx, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif
  int field = mat->get_fieldOrdering();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = m3dc1_mesh::instance()->mesh->getDimension();
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, ielm_dim, *ielm);
  assert(e);
  if (m3dc1_mesh::instance()->mesh->isGhost(e)) return M3DC1_FAILURE;
  
  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;

  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowIdx<numVar && *columnIdx<numVar);
  PetscInt rows[1024], columns[1024];
  assert(sizeof(rows)/sizeof(int)>=dofPerVar*nodes_per_element);
  if (mat->get_type()==0)
  {
    matrix_mult* mmat = dynamic_cast<matrix_mult*> (mat);
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      if (mmat->is_mat_local()) m3dc1_ent_getlocaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      else m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for (int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
  }
  else
  {
    matrix_solve* smat = dynamic_cast<matrix_solve*> (mat);
    int nodeOwner[6];
    PetscInt columns_bloc[6], rows_bloc[6];
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnIdx;
      for (int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      if (nodeOwner[inode]!=PCU_Comm_Self()&&!m3dc1_solver::instance()->assembleOption)
        smat->add_blockvalues(1, rows_bloc+inode, nodes_per_element, columns_bloc, values+offset);
      else 
        smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      offset+=numValuesNode;
    }
  }
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_matrix_write(int* matrix_id, const char* filename, int* start_index)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!filename)
    return m3dc1_matrix_print(matrix_id);

#ifdef DEBUG
  if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<": matrix id "<<*matrix_id<<", file "<<filename<<"\n";
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif

  char matrix_filename[256];
  sprintf(matrix_filename,"%s-%d", filename, PCU_Comm_Self());
  FILE * fp =fopen(matrix_filename, "w");

  int row, csize, sum_csize=0, index=0;

  vector<PetscInt> rows;
  vector<int> n_cols;
  vector<PetscInt> cols;
  vector<double> vals;

  mat->get_values(rows, n_cols, cols, vals);
  for (int i=0; i<rows.size(); ++i)
    sum_csize += n_cols[i];
  assert(vals.size()==sum_csize);

  fprintf(fp, "%lu\t%lu\t%lu\n", rows.size(), n_cols.size(), vals.size());

  for (int i=0; i<rows.size(); ++i)
  {
    row = rows[i];
    csize = n_cols[i];
    for (int j=0; j<csize; ++j)
    {
      fprintf(fp, "%d\t%d\t%E\n", row+*start_index, cols[index]+*start_index,vals[index]);
      ++index;
    }
  }
  fclose(fp);  
  assert(index == vals.size());
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_print(int* matrix_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
#ifdef DEBUG
  if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<": matrix id "<<*matrix_id<<"\n";
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" does not exist\n";
    return M3DC1_FAILURE;
  }
#endif

  int row, csize, sum_csize=0, index=0;

  vector<PetscInt> rows;
  vector<int> n_cols;
  vector<PetscInt> cols;
  vector<double> vals;

  mat->get_values(rows, n_cols, cols, vals);
  for (int i=0; i<rows.size(); ++i)
    sum_csize += n_cols[i];
  assert(vals.size()==sum_csize);

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": printing matrix "<<*matrix_id<<"\n";

  for (int i=0; i<rows.size(); ++i)
  {
    row = rows[i];
    csize = n_cols[i];
    for (int j=0; j<csize; ++j)
    {
      std::cout<<"["<<PCU_Comm_Self()<<"]\t"<<row<<"\t"<<cols[index]<<"\t"<<vals[index]<<"\n";
      ++index;
    }
  }
  assert(index == vals.size());
  return M3DC1_SUCCESS;
}

// obsolete
//*******************************************************
int m3dc1_matrix_setassembleoption(int * op)
//*******************************************************
{
  return M3DC1_FAILURE;
}
#else
#include "m3dc1_ls.h"
#endif // #ifdef M3DC1_PETSC

//*******************************************************
int m3dc1_field_sum_plane (FieldID* /* in */ field_id)
//*******************************************************
{
  MPI_Comm icomm= m3dc1_model::instance()->getMPICommPlane();
  int num_vtx=m3dc1_mesh::instance()->mesh->count(0), num_dof=0;
  m3dc1_field_getnumlocaldof(field_id, &num_dof);
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);

  int data_size=total_num_dof*num_vtx*(1+value_type);
  assert(total_num_dof*num_vtx==num_dof);
  double * thevec = NULL;
  m3dc1_field_getdataptr(field_id, &thevec);
  double * sendbuf = new double [data_size];
  m3dc1_field_retrieve (field_id, sendbuf, &num_dof);
  MPI_Allreduce (sendbuf,thevec,data_size,MPI_DOUBLE,MPI_SUM,icomm) ;
  synchronize_field((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());
  delete [] sendbuf;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_spr_adapt (FieldID* field_id, int* index, int* ts,
    double* ar, double* max_size, int* refine_level, 
    int* coarsen_level, bool* update)
//*******************************************************
{
  char filename[256];
#ifdef DEBUG
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<" field id "<<*field_id<<" , name "<<get_field_name_from_id(*field_id)<<"\n";
#endif

  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  int np = m3dc1_model::instance()->num_plane;

  // in_filed will hold all the dofs of all the fields (num being the total number of fields)
  // at each vertex. e.g.
  // f1_1, f1_2, f1_3, f1_4, f1_5, f1_6, ! dofs of 1st field
  // f2_1, f2_2, f2_3, f2_4, f2_5, f2_6, ! dofs of 2nd field
  // ...
  //
  // findex_1, findex_2, findex_3, findex_4, findex_5, findex_6, ! dofs of index'th field
  // ...
  // fnum_1, fnum_2, fnum_3, fnum_5, fnum_5, fnum_6 ! dofs of num'th (last) field
  apf::Field* inField = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  PCU_ALWAYS_ASSERT_VERBOSE(inField, "pointer is empty!");
  if (!PCU_Comm_Self())
    std::cout << "received field with name " << apf::getName(inField) << "to run spr on" << std::endl;

  // the following call will extract the ones at index
  apf::Field* targetField = get_field_at_index(mesh, inField, *index, dofNode);

  // transfer targetField (for spr_computations) and all the fields
  // that need to be transfered for the next solve step  onto the master-plane
  // the vector pFields is only used in 3D
  // Preprocessing 3D mesh. This involves the following steps
  // 1- copy the fields that are needed for solution transfer onto the master plane
  // 2- convert the mesh to 2D
  

  // the vector pFields and zFields are only used for 3D
  // pFields holds the multi-plane fields that need to be transfered during adapt
  // zFields holds the other fields so they can be zero-ed out after adapt+3D mesh reconstruction
  std::vector<MultiField> pFields;
  pFields.clear();

  std::vector<apf::Field*> zFields;
  zFields.clear();

  if (np > 1) // 3D
  {
    MultiField targetMultiField;
    transfer_field_to_main(targetField, targetMultiField);
    pFields.push_back(targetMultiField);
    mesh->removeField(targetField);
    apf::destroyField(targetField);

    std::map<FieldID, m3dc1_field*>::iterator it = m3dc1_mesh::instance()->field_container->begin();
    while(it!=m3dc1_mesh::instance()->field_container->end())
    {
      apf::Field* field = it->second->get_field();
      int complexType = it->second->get_value_type();
      if (complexType) group_complex_dof(field, 1);
      if (isFrozen(field)) unfreeze(field);
      if (it->second->should_transfer())
      {
        MultiField mf;
        transfer_field_to_main(field, mf);
        pFields.push_back(mf);
      }
      else
        zFields.push_back(field);
      it++;
    }
    for (int i = 0; i < (int)targetMultiField.size() ; i++)
      zFields.push_back(targetMultiField[i]);

    m3dc1_mesh::instance()->remove3D();
    m3dc1_mesh::instance()->rebuildPointersOnNonMasterPlane(pFields, zFields);
    mesh = m3dc1_mesh::instance()->mesh;
  }

  int valueType = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_value_type();
  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();

  if (m3dc1_model::instance()->num_plane == 1) // 2D
  {
    while(it!=m3dc1_mesh::instance()->field_container->end())
    {
      apf::Field* field = it->second->get_field();
      int complexType = it->second->get_value_type();
      assert(valueType==complexType);
      if (complexType) group_complex_dof(field, 1);
      if (isFrozen(field)) unfreeze(field);
      if (it->second->should_transfer())
      {
        if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": field with name "<<getName(field)<< " is added to solution transfer\n";
        fields.push_back(field);
      }
      it++;
    }
  }
  else // 3D
  {
// ignoring the first one since that is only used for size calculations
// hence the reason for the following loop start from 1
    for (int i = 1; i < (int)pFields.size(); i++) {
      PCU_ALWAYS_ASSERT((int)pFields[i].size() == np);
      for (int j = 0; j < np; j++)
        fields.push_back(pFields[i][j]);
    }
  }

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": "<<fields.size()
             <<" fields have been added to solution transfer\n";
  for (int i = 0; i < (int)fields.size(); i++) 
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__<< " field with name " << apf::getName(fields[i]) 
               <<  " fields have been added to solution transfer\n";

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    apf::destroyNumbering(n);
  }  

  // compute the size field here and remove the ip and target fields afterwards
  apf::Field* size_field = 0;
  ma::Input* in;

  if (np == 1) // 2D
  {
    apf::Field* targetField0 = get_vectorComponent_of_field(mesh, targetField);
    apf::Field* ip = spr::getGradIPField(targetField0, "ip", 2);

    size_field = spr::getSPRSizeField(ip, *ar);
    process_size_field(mesh, size_field, *ts, *max_size, *refine_level, *coarsen_level);
    fields.push_back(size_field);

    mesh->removeField(ip);
    mesh->removeField(targetField);
    mesh->removeField(targetField0);
    destroyField(ip);
    destroyField(targetField);
    destroyField(targetField0);

    ReducedQuinticImplicit shape;
    ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
#ifdef OLDMA
    in = ma::configure(mesh,size_field,&slnTrans);
#else
    in = ma::makeAdvanced(ma::configure(mesh, size_field, &slnTrans));
#endif

    in->shouldSnap=false;
    in->shouldTransferParametric=false;

    if (m3dc1_model::instance()->snapping)
    {
      in->shouldSnap=true;
      in->shouldTransferParametric=true;
      if (!PCU_Comm_Self())
        std::cout<<"[M3D-C1 INFO] "<<__func__<<" snapping turned on\n";
    }
    in->shouldRunPostZoltan = true;
    in->goodQuality = 0.5;
    in->maximumIterations = (*refine_level) + 1;

    if (*coarsen_level < 0)
      in->shouldCoarsen=false;

    ma::adapt(in);
    mesh->removeField(size_field);
    apf::destroyField(size_field);
    apf::reorderMdsMesh(mesh);

    m3dc1_mesh::instance()->initialize();
    compute_globalid(m3dc1_mesh::instance()->mesh, 0);
    compute_globalid(m3dc1_mesh::instance()->mesh, m3dc1_mesh::instance()->mesh->getDimension());
  } 
  else // 3D
  {
    MPI_Comm groupComm;
    int groupSize = PCU_Comm_Peers()/np;
    int lpid = PCU_Comm_Self()/groupSize;
    int grnk = PCU_Comm_Self()%groupSize;
    MPI_Comm_split(m3dc1_model::instance()->oldComm, lpid, grnk, &groupComm);
    PCU_Switch_Comm(groupComm);
    if (m3dc1_model::instance()->local_planeid == 0)
    { 
      if (*update)
      { 
        int s = 1;
        int ns = -1;
        m3dc1_mesh_write("part", &s, &ns);
      }
      size_field = compute_multiplane_size_field(pFields[0], *ar);
      process_size_field(mesh, size_field, *ts, *max_size, *refine_level, *coarsen_level);
      fields.push_back(size_field);
      
      ReducedQuinticImplicit shape;
      ReducedQuinticTransfer slnTrans(mesh,fields, &shape);

#ifdef OLDMA
      ma::Input* in = ma::configure(mesh,size_field,&slnTrans);
#else
      ma::Input* in = ma::makeAdvanced(ma::configure(mesh, size_field, &slnTrans));
#endif
      
      in->shouldSnap=false;
      in->shouldTransferParametric=false;

      if (m3dc1_model::instance()->snapping)
      {
        in->shouldSnap=true;
        in->shouldTransferParametric=true;
        if (!PCU_Comm_Self())
          std::cout<<"[M3D-C1 INFO] "<<__func__<<" snapping turned on\n";
      }
      in->shouldFixShape = true;
      in->shouldRunPostZoltan = true;
      in->goodQuality = 0.5;
      in->maximumIterations = (*refine_level);
      
      if (*coarsen_level < 0)
        in->shouldCoarsen=false;
      
      ma::adapt(in);
      
      for (int i = 0; i < (int)zFields.size(); i++)
        zero_fields_on_master(zFields[i]);
      
      mesh->removeField(size_field);
      apf::destroyField(size_field);
      
      while (mesh->countNumberings())
      { 
        apf::Numbering* n = mesh->getNumbering(0);
        mesh->removeNumbering(n);
        apf::destroyNumbering(n);
      }
      
      for (int i = 1; i < (int)pFields.size(); i++)
        for (int j = 0; j < np; j++)
          synchronize_field(pFields[i][j]);
      apf::reorderMdsMesh(mesh);
    }

    PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
    MPI_Comm_free(&groupComm);

    m3dc1_mesh::instance()->restore3D();

    for (int i = 1; i < (int)pFields.size(); i++)
      zero_fields_on_non_master(pFields[i]);

    for (int i = 1; i < (int)pFields.size(); i++)
      transfer_field_from_main(pFields[i]);

    for (int i = 1; i < (int)pFields.size(); i++)
      for (int j = 0; j < (int)pFields[i].size(); j++) {
        mesh->removeField(pFields[i][j]);
        apf::destroyField(pFields[i][j]);
      }

    for (int i = 0; i < (int)zFields.size(); i++)
      apf::zeroField(zFields[i]);

    for (int i = 0; i < (int)pFields[0].size(); i++) {
      mesh->removeField(pFields[0][i]);
      apf::destroyField(pFields[0][i]);
    }
  }

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    if (!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first;
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    synchronize_field(field);

#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
}

int adapt_time=0;
int adapt_by_field (int * fieldId, double* psi0, double * psil)
{
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] running adaptation by post processed magnetic flux field\n";


  int shouldCoarsen = -1;
  set<int> field_keep;
  field_keep.insert(*fieldId);
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  std::ifstream ifs;
  ifs.open("sizefieldParam", std::ifstream::in);
  if (!ifs.is_open())
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] file \"sizefieldParam\" not found\n";
    return M3DC1_FAILURE;
  }

  // Note:
  // count will hold the total number of parameters read from "sizefieldParam" and it must be either 13 or 14. The expected
  // behaviour for each case is summarized below
  // a) if count == 13, coarsening will stay on
  // b) if count == 14, coarsening will be turned off (the actual value of the last (14th) parameter does not have an effect.
  double param[14];
  int count = 0;
  double p;
  while (ifs >> p && count<14)
  {
    param[count] = p;
    count++;
  }
  ifs.close();

  if (count==13)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] read 13 parameters from \"sizefieldParam\": coarsening will stay on\n";
    // turn on coarsening for this case
    shouldCoarsen = 1;
  }
  else if (count==14)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] read 14 parameters from \"sizefieldParam\"\n";
    shouldCoarsen = (int)param[13];
    if (!PCU_Comm_Self())
    {
      if (shouldCoarsen == 0)
	std::cout<<"[M3D-C1 INFO] the last parameter in \"sizefieldParam\" is 0 causing coarsening to be turned off\n";
      else if (shouldCoarsen == 1)
	std::cout<<"[M3D-C1 INFO] the last parameter in \"sizefieldParam\" is 1 causing coarsening to be turned on\n";
      else
      {
        std::cout<<"[M3D-C1 ERROR] the 14th value in \"sizefieldParam\" (it present) must be 0. or 1.\n";
        return M3DC1_FAILURE;
      }
    }
  }
  else
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] number of parameters in \"sizefieldParam\" must be 13 or 14\n";
    return M3DC1_FAILURE;
  }

  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] size field parameters: ";
  for (int i=0; i<count; ++i)
  {
    if (!PCU_Comm_Self()) std::cout<<std::setprecision(5)<<param[i]<<" ";
  }
  if (!PCU_Comm_Self()) std::cout<<"\n";

  apf::Field* psiField = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_field();

  // delete all the matrix
#ifdef M3DC1_TRILINOS
  while (m3dc1_ls::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_epetra*>::iterator mat_it = m3dc1_ls::instance()->matrix_container->begin();
    mat_it->second->destroy();
    delete mat_it->second;
    m3dc1_ls::instance()->matrix_container->erase(mat_it);
  }
#endif
#ifdef M3DC1_PETSC
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }
#endif

  int valueType = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_value_type();
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);
  double mmax[2], mmin[2];
  m3dc1_model_getmaxcoord(mmax,mmax+1);
  m3dc1_model_getmincoord(mmin,mmin+1);

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    assert(valueType==complexType);
    if (complexType) group_complex_dof(field, 1);
    if (isFrozen(field)) unfreeze(field);
    fields.push_back(field);
    it++;
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    apf::destroyNumbering(n);
  }
  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
#ifdef OLDMA
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
#else
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh,&sf,&slnTrans));
#endif
  in->maximumIterations = 9;

  in->shouldSnap=false;
  in->shouldTransferParametric=false;

  if (m3dc1_model::instance()->snapping)
  {
    in->shouldSnap=true;
    in->shouldTransferParametric=true;
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__<<" snapping turned on\n";
  }

#ifdef DISABLE_ZOLTAN
  in->shouldRunPostZoltan = false;
#else
  in->shouldRunPostZoltan = true;
#endif

  // set the coarsening options
  if (shouldCoarsen == 1)
    in->shouldCoarsen = true;
  else if (shouldCoarsen == 0)
    in->shouldCoarsen = false;
  else
    return M3DC1_FAILURE;


  ma::adapt(in);
  reorderMdsMesh(mesh);

  m3dc1_mesh::instance()->initialize();
  compute_globalid(m3dc1_mesh::instance()->mesh, 0);
  compute_globalid(m3dc1_mesh::instance()->mesh, m3dc1_mesh::instance()->mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    if (!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    synchronize_field(field);
#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }

  if (!PCU_Comm_Self())
    cout<<"#global_ent: V "<<m3dc1_mesh::instance()->num_global_ent[0]
        <<", E "<<m3dc1_mesh::instance()->num_global_ent[1]
        <<", F "<<m3dc1_mesh::instance()->num_global_ent[2]
        <<", R "<<m3dc1_mesh::instance()->num_global_ent[3]<<"\n";

  return M3DC1_SUCCESS;
}

double absSize[2]={0,1}, relSize[2]={0.3, 1.5};
int set_mesh_size_bound (double* abs_size, double * rel_size)
{
  for (int i=0; i<2; i++)
  {
    absSize[i] =  abs_size[i];
    relSize[i] = rel_size[i];
  }
  return M3DC1_SUCCESS;
}

void smooth_size_field (apf::Field* sizeField)
{
  int numVert=m3dc1_mesh::instance()->mesh->count(0);

  for (int i=0; i<numVert; i++)
  {
    vector<apf::MeshEntity*> nodes;
    apf::MeshEntity* e = apf:: getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    nodes.push_back(e);
    double sizeOrg=0;
    getComponents(sizeField, e, 0, &sizeOrg);
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,1,adjacent);
    for (int i=0; i<adjacent.getSize(); i++)
    {
      apf::Downward downward;
      m3dc1_mesh::instance()->mesh->getDownward(adjacent[i], 0, downward);
      nodes.push_back(downward[0]==e?downward[1]:downward[0]);
    }
    double size=0;
    for (int i=0; i<nodes.size(); i++)
    {
      double buff;
      getComponents(sizeField, nodes[i], 0, &buff);
      size+=1./buff;
    }
    size/=nodes.size();
    size=1./size;
    //if (size<0.85*sizeOrg)
    setComponents(sizeField, e, 0, &size);
  }
}

double p=4;
int set_adapt_p (double * pp) 
{
  p=*pp;
  return M3DC1_SUCCESS;
}

// Arguments of the function
// double* node_error: Node error data coming from the function node_error_3d_mesh()
// int num_planes: User defined number of planes
//    Must be greater than 1 to extrude the meshes between the planes
// double* errorAimed: Parameter "adapt_target_error" from the user input parameter file
//    Target discretization error on the adapted mesh
// double* max_adapt_node: Parameter "adapt_max_node" from the user input parameter file
//      Maximum node number in the adapted mesh. If the estimated mesh node number from adapt_target_error exceeds iadapt_max_node,
//    the target mesh size in the adapted mesh is scaled such that the mesh node number is below iadapt_max_node.
// int* option: Parameter "adapt_control" from the user input parameter file and is either 0 or 1
//    0: adapt_target_error is global (integral over the domain)
//    1: adapt_target_error is local (integral over the element)
int adapt_by_error_field (double * errorData, double * errorAimed, int * max_adapt_node, int * option)
{
  if (!PCU_Comm_Self()) 
  std::cout<<"[M3D-C1 INFO] running adaptation by error estimator\n";

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::Field* sizeField = createPackedField(m3dc1_mesh::instance()->mesh, "size_field", 1);
  SizeFieldError sf (m3dc1_mesh::instance()->mesh, sizeField, *errorAimed);

  int numVert=m3dc1_mesh::instance()->mesh->count(0);

  // first sum error ^ (2d/(2p+d))
  double d=2;
  double errorSum=0;
  if (*option)
  {
    for (int i=0; i<numVert; i++)
    {
      // Determine if the mesh is on the original poloidal plane or on ghost plane
      if (is_ent_original(m3dc1_mesh::instance()->mesh,getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i)))
        errorSum+=pow(errorData[i],d/(p+d/2.0));
    }
    double errorSumBuff=errorSum;
    MPI_Allreduce(&errorSumBuff, &errorSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errorSum = *errorAimed*(*errorAimed)/errorSum;
    errorSum = pow(errorSum,1./(2.*p));
  }
  else // default 
    errorSum=pow(*errorAimed,1./(p+d/2.));

  double size, targetSize, size_estimate=0;
  for (int i=0; i<numVert; ++i)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    // if not on original poloidal plane, ignore it
    if (!is_ent_original(mesh,e)) continue;
    size = sf.getSize(e);
    targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    size_estimate+=max(1.,1./targetSize/targetSize);
  }

  double size_estimate_buff=size_estimate;
  MPI_Allreduce(&size_estimate_buff, &size_estimate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int numNodeGlobl=0, dim=0;
  m3dc1_mesh_getnumglobalent(&dim, &numNodeGlobl);
  if (!pumi_rank())
    cout<<__func__<<": numVert "<<numNodeGlobl<<" size_estimate "<<size_estimate<<"\n";

  if (size_estimate>*max_adapt_node) errorSum*=sqrt(size_estimate/(*max_adapt_node));
  std::vector <double> target_size;
  int num_planes = m3dc1_model::instance()->num_plane;
  if (num_planes>1) // 3D
    target_size.resize(numVert);

  for (int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    size = sf.getSize(e);
    assert(errorData[i]==errorData[i]);
    targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    if (targetSize>relSize[1]) targetSize=relSize[1]; // not too much coarsening
    if (targetSize<relSize[0]) targetSize=relSize[0]; // not too much refining
    targetSize*=size;
    if (targetSize>absSize[1]) targetSize=absSize[1];
    if (targetSize<absSize[0]) targetSize=absSize[0];
    if (num_planes>1)
      target_size.push_back(targetSize);
    else // 2D 
      setComponents(sizeField, e, 0, &targetSize);
  }

  if (num_planes>1) // 3D
  {
    std::vector <double> target_sizefield;
    double final_target;
    // seol #3: if #vertex is N, #vertex in 3D is 2N. what is num_vert_on_plane for?
    int num_vert_on_plane = numVert/num_planes;
    for (int j=0; j<num_vert_on_plane; ++j)
    {
      for (int k=1; k<num_planes; ++k)
      {
        final_target = target_size[j];
        if (target_size[j+num_vert_on_plane*k] < final_target)
          final_target = target_size[j+num_vert_on_plane*k];
      }
      target_sizefield.push_back(final_target);     // This is the targetted sizefield for first poloidal plane
    }
  } // if (num_planes>1)
 
  // FIXME: why called twice?
  smooth_size_field(sizeField);
  smooth_size_field(sizeField);
  synchronize_field(sizeField);

  // delete all the matrix
#ifdef M3DC1_TRILINOS
  while (m3dc1_ls::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_epetra*>::iterator mat_it = m3dc1_ls::instance()->matrix_container->begin();
    mat_it->second->destroy();
    delete mat_it->second;
    m3dc1_ls::instance()->matrix_container->erase(mat_it);
  }
#endif
#ifdef M3DC1_PETSC
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }
#endif

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  fields.push_back(sizeField);
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 1);
    if (isFrozen(field)) unfreeze(field);
    //if (!PCU_Comm_Self()) std::cout<<"Solution transfer: add field "<<apf::getName(field)<<std::endl;
    fields.push_back(field);
    it++;
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    apf::destroyNumbering(n);
  }
  char filename[256];
  sprintf(filename,"before%d",adapt_time);
  //apf::writeVtkFiles(filename,mesh);

  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
#ifdef OLDMA
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
#else
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh,&sf,&slnTrans));
#endif
  in->maximumIterations = 5;
  in->shouldSnap=false;
  in->shouldTransferParametric=false;

  if (m3dc1_model::instance()->snapping)
  {
    in->shouldSnap=true;
    in->shouldTransferParametric=true;
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__<<" snapping turned on\n";
  }

#ifdef DISABLE_ZOLTAN
  in->shouldRunPostZoltan = false;
#else
  in->shouldRunPostZoltan = true;
#endif

  ma::adapt(in);
  reorderMdsMesh(mesh);

  m3dc1_mesh::instance()->initialize();
  compute_globalid(m3dc1_mesh::instance()->mesh, 0);
  compute_globalid(m3dc1_mesh::instance()->mesh, m3dc1_mesh::instance()->mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    if (!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    synchronize_field(field);

#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
  destroyField(sizeField);

  if (!PCU_Comm_Self())
    cout<<"#global_ent: V "<<m3dc1_mesh::instance()->num_global_ent[0]
        <<", E "<<m3dc1_mesh::instance()->num_global_ent[1]
        <<", F "<<m3dc1_mesh::instance()->num_global_ent[2]
        <<", R "<<m3dc1_mesh::instance()->num_global_ent[3]<<"\n";

  return M3DC1_SUCCESS;
}

int m3dc1_field_printcompnorm(FieldID* /* in */ field_id, char* info)
{
  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  int num_dof;
  m3dc1_field_getnumlocaldof(field_id, &num_dof);
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];

  apf::Field* f = mf ->get_field();
  int dof_per_node = countComponents(f);
  int num_comp=dof_per_node/C1TRIDOFNODE;
  vector<double> norms(dof_per_node/C1TRIDOFNODE);
  int j=0;
  for (int i=0; i<num_dof/C1TRIDOFNODE; i++)
  {
     for (int k=0; k<6; k++)
        norms.at(j)+=pts[i*C1TRIDOFNODE+k]*pts[i*C1TRIDOFNODE+k];
     j++;
     j%=num_comp;
  }
  vector<double> buff=norms;
  MPI_Allreduce(&buff[0],&norms[0], num_comp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  int psize;
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  if (PCU_Comm_Self() == psize-1)
  {
    std::cout<< "norm of vec "<<info;
    for (int i=0; i<num_comp; i++)
      std::cout<<" "<<std::sqrt(norms[i]);
    std::cout<<std::endl;
  }
  return M3DC1_SUCCESS;
}

int m3dc1_mesh_write(char* filename, int *option, int* timestep)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  char filename_buff[256];
  // vtk
  if (*option==0 ||*option==3)
  {
    if (*timestep >= 0)
      sprintf(filename_buff, "ts%04d-%s", *timestep, filename);
    else
      sprintf(filename_buff, "%s",filename);

    apf::MeshEntity* e;
    int dim=2, num_ent=mesh->count(2);

    vector<double> geoId (num_ent);
    apf:: MeshIterator* it = mesh->begin(dim);
    while ((e = mesh->iterate(it)))
    {
      int ent_id = getMdsIndex(mesh, e);
      int geom_class_dim,geom_class_id;
      m3dc1_ent_getgeomclass (&dim, &ent_id, &geom_class_dim, &geom_class_id);
      geoId.at(ent_id)=geom_class_id;
    }
    mesh->end(it);

    apf::writeVtkFiles(filename_buff, mesh);
    int one=1;
    if (*option==3) output_face_data (&one, &geoId[0], "geoId");

    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__<<": vtk folder \""<<filename_buff<<"\"\n";
  }
  else
  {
    if (*timestep >= 0)
      sprintf(filename_buff, "ts%d-%s.smb",*timestep,filename);
    else
      sprintf(filename_buff, "%s.smb",filename);

    int fieldID=12;
    double dofBuff[1024];
    m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[fieldID];
    apf::Field* f = mf ->get_field();
    int numDof = countComponents(f);
    apf::MeshTag* tag = mesh->createDoubleTag("field12", numDof);
    apf::MeshEntity* e;
    const int dim=0;
    apf:: MeshIterator* it = mesh->begin(dim);
    while ((e = mesh->iterate(it)))
    {
      apf::getComponents(f, e, 0, dofBuff);
      mesh->setDoubleTag(e,tag, dofBuff);
    }
    mesh->end(it);

    std::vector<apf::Field*> allFields;
    allFields.clear();
    while (mesh->countFields())
    {
      apf::Field* f = mesh->getField(0);
      apf::freeze(f);
      allFields.push_back(f);
      mesh->removeField(f);
    }
    mesh->writeNative(filename_buff);
    for (std::size_t i = 0; i < allFields.size(); i++) {
      apf::unfreeze(allFields[i]);
      mesh->addField(allFields[i]);
    }

    apf::removeTagFromDimension(mesh, tag, dim);
    mesh->destroyTag(tag);
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__<<": file \""<<filename_buff<<"\"\n";
  }
  return M3DC1_SUCCESS;
}

void print_mesh_info (apf::Mesh2* m)
{
  if (!PCU_Comm_Self()) std::cout<<"\n===== mesh size =====\n";

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  int* own_entity_count = new int[4*PCU_Comm_Peers()];

  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=own_entity_count[i]=0;

  apf::MeshEntity* e;
  int self = PCU_Comm_Self();

  for (int d=0; d<4;++d)
  {
    local_entity_count[4*self+d] = m->count(d);
    apf::MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it)))
    {
      if (m->getOwner(e)==self)
        ++own_entity_count[4*PCU_Comm_Self()+d];
    }
    m->end(it);
  }

  int* global_local_entity_count = new int[4*PCU_Comm_Peers()];
  int* global_own_entity_count = new int[4*PCU_Comm_Peers()];

  MPI_Allreduce(local_entity_count, global_local_entity_count, 4*PCU_Comm_Peers(),
                MPI_INT, MPI_SUM, PCU_Get_Comm());

  MPI_Allreduce(own_entity_count, global_own_entity_count, 4*PCU_Comm_Peers(),
                MPI_INT, MPI_SUM, PCU_Get_Comm());

  if (!PCU_Comm_Self())
  {
    int* global_entity_count = new int[4];
    global_entity_count[0]=global_entity_count[1]=global_entity_count[2]=global_entity_count[3]=0;
    for (int d=0; d<4;++d)
    {
      for (int p=0; p<PCU_Comm_Peers();++p)
        global_entity_count[d] += global_own_entity_count[p*4+d];
    }

    for (int p=0; p<PCU_Comm_Peers(); ++p)
      std::cout<<"(p"<<p<<") # local ent: v "<<global_local_entity_count[p*4]
        <<", e "<<global_local_entity_count[p*4+1]
        <<", f "<<global_local_entity_count[p*4+2]
        <<", r "<<global_local_entity_count[p*4+3]<<"\n";
    std::cout<<"\n";
    for (int p=0; p<PCU_Comm_Peers(); ++p)
      if (global_own_entity_count[p*4])
        std::cout<<"(p"<<p<<") # own ent: v "<<global_own_entity_count[p*4]
          <<", e "<<global_own_entity_count[p*4+1]
          <<", f "<<global_own_entity_count[p*4+2]
          <<", r "<<global_own_entity_count[p*4+3]<<"\n";
    std::cout<<"\n";

    std::cout<<"# global ent: v "<<global_entity_count[0]<<", e "<<global_entity_count[1]
              <<", f "<<global_entity_count[2]<<", r "<<global_entity_count[3]<<"\n\n";

    delete [] global_entity_count;

    }

  delete [] local_entity_count;
  delete [] global_local_entity_count;
  delete [] own_entity_count;
  delete [] global_own_entity_count;
}

void m3dc1_mesh_verify()
{
  apf::verify(m3dc1_mesh::instance()->mesh, false);
  print_mesh_info(m3dc1_mesh::instance()->mesh);
}

int sum_edge_data (double * data, int* size)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int num_edge=m3dc1_mesh::instance()->mesh->count(1), edg_dim=1;

  PCU_Comm_Begin();
  for (int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(data[i*(*size)]),(*size)*sizeof(double));
  }

  PCU_Comm_Send();
  double* receive_buff = new double [*size];
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]+=receive_buff[i];
    }

  PCU_Comm_Begin();
  for (int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;

    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]=receive_buff[i];
    }
  delete [] receive_buff;
  return M3DC1_SUCCESS;
}

int get_node_error_from_elm (double * elm_data, int * size, double* nod_data)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int num_node=m3dc1_mesh::instance()->mesh->count(0);
  int num_elm=m3dc1_mesh::instance()->mesh->count(2);
  int nod_dim=0;

  PCU_Comm_Begin();
  double* buff = new double[*size];
  double* area = new double[num_node];
  for (int i=0; i<num_node; i++)
    area[i]=0.;
  for (int i=0; i<(*size)*num_elm; i++ )
    for (int j=0; j<*size; j++)
    {
      assert(elm_data[(*size)*i+j]==elm_data[(*size)*i+j]);
      assert(elm_data[(*size)*i+j]>=0);
    }
  for (int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
    apf::Adjacent adjacent;
    m->getAdjacent(e,2,adjacent);
    for (int j=0; j<adjacent.getSize(); j++)
    {
       apf::MeshElement* me = createMeshElement(m, adjacent[j]);
       double s = apf::measure(me);
       int ielm = getMdsIndex(m, adjacent[j]);
       assert(ielm>=0 &&ielm<num_elm);
       for (int k=0; k<*size; k++)
          nod_data[i*(*size)+k]+=s*elm_data[(*size)*ielm+k];
       area[i]+=s;
       destroyMeshElement(me);
    }

    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_COMM_PACK(own_partid, area[i]);
    PCU_Comm_Pack(own_partid, &(nod_data[(*size)*i]), sizeof(double)*(*size));
  }

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      int inode = getMdsIndex(m, node);
      double s;
      PCU_COMM_UNPACK(s);
      area[inode]+=s;
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]+=buff[i];
    }

  for (int i=0; i<num_node; i++)
  {
    for (int j=0; j<*size; j++)
      nod_data[i*(*size)+j]/=area[i];
  }
  PCU_Comm_Begin();
  for (int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(nod_data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      int inode = getMdsIndex(m, node);
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]=buff[i];
    }
  delete [] buff;
  delete [] area;  
  return M3DC1_SUCCESS;
}

int m3dc1_field_max (FieldID* field_id, double * max_val, double * min_val)
{
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (mf->get_value_type()) dofPerEnt *= 2;
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=m3dc1_mesh::instance()->mesh->count(0);
  int vertex_type=0;

  std::vector<double> maxVal(dofPerEnt, -1e30), minVal(dofPerEnt,1e30), dofs(dofPerEnt);
  for (int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for (int i=0; i<dofPerEnt; i++)
    {
      if (maxVal[i]<dofs[i]) maxVal[i]=dofs[i];
      if (minVal[i]>dofs[i]) minVal[i]=dofs[i];
    }
  }
  MPI_Allreduce(&(maxVal[0]), max_val, dofPerEnt, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(minVal[0]), min_val, dofPerEnt, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return M3DC1_SUCCESS;
}

void setSizeFieldOnVertex(ma::Entity* mV, double& xSize, double& ySize,SizeFieldPsi sf, double* dirVector)
{ 
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int indx= getMdsIndex(m, mV);
  ma::Matrix M;
  ma::Vector H;
  sf.getValue(mV,M,H);
  xSize = H[0];
  ySize = H[1];
  dirVector[indx*3+0] = M[0][0];
  dirVector[indx*3+1] = M[0][1];
  dirVector[indx*3+2] = 0.0;
}

void sizeFieldTransition( int refIndx, int indx, double& xSize, double& ySize, double* dirVector,
     int* field_id1, int* field_id2, SizeFieldPsi sf, int mode)
{ 
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  double data1, data2;
  int size_data=1;
  m3dc1_node_getfield(&refIndx, field_id1, &data1, &size_data);
  m3dc1_node_getfield(&refIndx, field_id2, &data2, &size_data);
        
  if (mode == 0)
  {     
        xSize = data1*2.0;
        ySize = data2*2.0;
  }
  else if (mode == 1)
  {     
        xSize = data1*0.5;
        ySize = data2*0.5;
  }
  
  ma::Entity* mV = getMdsEntity(m, 0, indx);
  ma::Matrix M;
  ma::Vector H;
  sf.getValue(mV,M,H);
  dirVector[indx*3+0] = M[0][0];
  dirVector[indx*3+1] = M[0][1];
  dirVector[indx*3+2] = 0.0;
}

int setSizeFieldOnVertices(int* field_id1, int* field_id2, SizeFieldPsi sf, double* dir,int adaptFaceNumber)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  gmi_model* model = m3dc1_model::instance()->model;
  double data_1, data_2;
  int size_data=1;

  apf::MeshTag* sizeFieldSet = m->createIntTag("sfDefined", 1);
  std::vector <int> adaptModelFaces;
  adaptModelFaces.push_back(adaptFaceNumber);
  for (int nFace = 0; nFace < adaptModelFaces.size(); ++nFace)
  {
    int modelFaceIndx = adaptModelFaces[nFace];
    for (int i=0; i<m->count(0); ++i)
    {
      ma::Entity* mV = getMdsEntity(m, 0, i);
      gmi_ent* gent= (gmi_ent*)(m->toModel(mV));
      int gTag = gmi_tag(model, gent);
      int gType = gmi_dim(model,gent);
      int vTag = 0;
      if (gType == 2 && gTag == modelFaceIndx)
      {
        setSizeFieldOnVertex(mV,data_1, data_2, sf, dir);
        m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
        m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);
        vTag = 1;
        m->setIntTag(mV, sizeFieldSet, &vTag);
      }
      else
        m->setIntTag(mV, sizeFieldSet, &vTag);
    }
    std::vector <int> adaptModelEdges;
    get_gent_adj(2, modelFaceIndx, 1, adaptModelEdges);
    adaptModelEdges.push_back(1);
    for (int nEdge = 0; nEdge < adaptModelEdges.size(); ++nEdge)
    {
      int modelEdgeIndx =  adaptModelEdges[nEdge];
      gmi_ent* ge = gmi_find(model, 1, modelEdgeIndx);
      gmi_set* gvOnEdge;
      gvOnEdge = gmi_adjacent(model, ge, 0);
      int modelNodeIndx = gmi_tag(model,gvOnEdge->e[0]);
      std::vector <apf::MeshEntity*> meshV, meshVTemp;
      for (int i=0; i<m->count(0); ++i)
      {
        ma::Entity* mV = getMdsEntity(m, 0, i);
        gmi_ent* gent= (gmi_ent*)(m->toModel(mV));
        int gTag = gmi_tag(model, gent);
        int gType = gmi_dim(model,gent);
        if ((gType == 1 && gTag == modelEdgeIndx) || (gType == 0 && gTag == modelNodeIndx))
        {
          setSizeFieldOnVertex(mV,data_1, data_2, sf, dir);
          m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
          m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);
          m->removeTag(mV, sizeFieldSet);
          int vTag = 1;
           m->setIntTag(mV, sizeFieldSet, &vTag);
           meshV.push_back(mV);
        }
      }
      int iterCount = 0;
      while (meshV.size() > 0)
      {
        apf::MeshEntity* mV = meshV[0];
        double d_1, d_2;
        int vId= getMdsIndex(m, mV);
        m3dc1_node_getfield(&vId, field_id1, &d_1, &size_data);
        m3dc1_node_getfield(&vId, field_id2, &d_2, &size_data);
        double desiredLength = sqrt((d_1*d_1)+(d_2*d_2));
        apf::Adjacent edges;
        m->getAdjacent(mV,1, edges);
        for (int i = 0; i < edges.getSize(); ++i)
        {
          apf::MeshEntity* edge = edges[i];
          gmi_ent* gent= (gmi_ent*)(m->toModel(edge));
          int gTag = gmi_tag(model, gent);
          int gType = gmi_dim(model,gent);
          if ((gType == 2 && gTag == modelFaceIndx) || (gType == 1 && gTag == modelEdgeIndx))
            continue;
          else
          {
            ma::Entity* otherV = getEdgeVertOppositeVert(m, edge, mV);
            int vIndx = getMdsIndex(m, otherV);
            int vTagSet = -1;
            m->getIntTag(otherV, sizeFieldSet, &vTagSet);
            if (vTagSet == 1)
              continue;

            if (apf::measure(m,edge) > 1.5*desiredLength)
              sizeFieldTransition(vId, vIndx, data_1, data_2, dir, field_id1, field_id2,sf,0);
            else
              setSizeFieldOnVertex(otherV,data_1, data_2, sf, dir);
            m3dc1_node_setfield(&vIndx, field_id1, &data_1, &size_data);
            m3dc1_node_setfield(&vIndx, field_id2, &data_2, &size_data);
            m->removeTag(otherV, sizeFieldSet);
            int vTag = 1;
            m->setIntTag(otherV, sizeFieldSet, &vTag);
            meshVTemp.push_back(otherV);
          }
        }
        meshV.erase(meshV.begin());
        if (meshV.size() == 0 && iterCount < 2)
        {
                for (int i = 0; i < meshVTemp.size(); ++i)
                   meshV.push_back(meshVTemp[i]);
                meshVTemp.clear();
                iterCount++;
        }
      }
    }
  }
  for (int i=0; i<m->count(0); ++i)
  {
    ma::Entity* mV = getMdsEntity(m, 0, i);
    int vTagSet = -1;
    m->getIntTag(mV, sizeFieldSet, &vTagSet);

    if (vTagSet == 1)
      continue;

    setSizeFieldOnVertex(mV,data_1, data_2, sf, dir);
    m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
    m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);
  }

  m3dc1_field_sync(field_id1);
  m3dc1_field_sync(field_id2);

  return M3DC1_SUCCESS;
}

int adapt_model_face(int * fieldId, double* psi0, double * psil, int* iadaptFaceNumber)
{
  if (!PCU_Comm_Self())
  std::cout<<"[M3D-C1 INFO] running adaptation by post processed magnetic flux field\n";

  FILE *fp = fopen("sizefieldParam", "r");
  if (!fp)
  {
    std::cout<<"[M3D-C1 ERROR] file \"sizefieldParam\" not found\n";
    return M3DC1_FAILURE;
  }
  double param[13];
  set<int> field_keep;
  field_keep.insert(*fieldId);
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  gmi_model* model = m3dc1_model::instance()->model;

  gmi_ent* gf;
  gmi_iter* gfIter = gmi_begin(model, 2);
  bool faceFound = false;
  std::set <int> faceIds;
  while((gf = gmi_next(m3dc1_model::instance()->model, gfIter)) )
  {
    if (gmi_tag(model, gf) == *iadaptFaceNumber)
      faceFound = true;
    faceIds.insert(gmi_tag(model, gf));
  }
  if (!faceFound)
  {
    if (!PCU_Comm_Self())
    {
        std::cout << "*************************************** WARNING **********************************************\n";
        std::cout << "The model face to adapt (ID = " << *iadaptFaceNumber << ") given in C1Input doesn't exist in model.\n";
        std::cout << "Choose a face ID from the following range: "
                  << *next(faceIds.begin(),0) << " - " << *next(faceIds.begin(),faceIds.size()-1) << "\n";
        std::cout << "**********************************************************************************************\n";
    }
    exit(0);
  }
  else
  {
    if (!PCU_Comm_Self())
      std::cout << "The model face to adapt with ID " << *iadaptFaceNumber << " found in the model" << "\n";
  }
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] size field parameters: ";

  for (int i=0; i<13; ++i)
  {
    fscanf(fp, "%lf ", &param[i]);
    if (!PCU_Comm_Self()) std::cout<<std::setprecision(5)<<param[i]<<" ";
  }
  fclose(fp);
  if (!PCU_Comm_Self()) std::cout<<"\n";

  apf::Field* psiField = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_field();

  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  int valueType = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_value_type();
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);

  int fid_size1, fid_size2;
  m3dc1_field_getnewid (&fid_size1);

  fid_size2 = fid_size1+1;
  int scalar_type=0;
  int num_value=1;

  m3dc1_field_create (&fid_size1, "size 1", &num_value, &scalar_type, &num_value);
  m3dc1_field_create (&fid_size2, "size 2", &num_value, &scalar_type, &num_value);

  double* dir = new double[mesh->count(0)*3];
  int adaptFaceNumber = *iadaptFaceNumber;

  setSizeFieldOnVertices(&fid_size1, &fid_size2, sf, dir,adaptFaceNumber);

  apf::MeshTag* doNotAdaptFace = mesh->createIntTag("doNotAdapt", 1);
  for (int i=0; i<mesh->count(2); ++i)
  {
     apf::MeshEntity* ele = getMdsEntity(mesh, 2, i);
     gmi_ent* gent= (gmi_ent*)(mesh->toModel(ele));
     int gTag = gmi_tag(model, gent);
     int gType = gmi_dim(model,gent);
     int elemTag = -1;
     apf::MeshEntity*  vertices[3];
     mesh->getDownward(ele, 0, vertices);
     int count = 0;
     apf::MeshTag* vTag = mesh->findTag("sfDefined");
     for (int j =0; j < 3; ++j)
     {
       int tagNode = -1;
       mesh->getIntTag(vertices[j], vTag, &tagNode);
       if (tagNode == 1)
         count++;
     }
     if (count == 3)
       elemTag = 0;
     else
       elemTag = 1;

     mesh->setIntTag(ele, doNotAdaptFace, &elemTag);
  }

  m3dc1_mesh_adapt(&fid_size1, &fid_size2, dir);
  return M3DC1_SUCCESS;
}

// Snapping Operation during mesh adapt is not supported with current analytical model
// The new vertices after adapt are not reparameterized on the model edge since snapping 
// is disabled. Without parametric coordinates, the current methods return wrong normal 
// vectors. This function evaluates the normal vector on the new vertices at the boundary 
// after the mesh adapt.
// Will be obselete when snapping will be supported (in plans for near future)
void m3dc1_node_getNormVecOnNewVert(apf::MeshEntity* v, double* normalVec)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int entType = m->getType(v);
  assert(entType == 0);	
  
  bool v0_found = false;
  bool v1_found = false;
  apf::Adjacent edges;
  m->getAdjacent(v, 1, edges);
  std::vector <apf::MeshEntity*> adjVert;
  for (int i = 0; i < edges.getSize(); ++i)
  { 
    apf::MeshEntity* edge = edges[i];
    apf::MeshEntity* otherV = getEdgeVertOppositeVert(m, edge, v);
    gmi_ent* gent= (gmi_ent*)(m->toModel(otherV));
    int gType = gmi_dim(m3dc1_model::instance()->model,gent);
    if (gType != 1)
      continue;
    adjVert.push_back(otherV);
  }
  apf::Vector3 pt0(0.0,0.0,0.0);
  apf::Vector3 pt1(0.0,0.0,0.0);
  apf::MeshEntity* v1 = adjVert[0]; 
  m->getPoint(v, 0, pt0);
  m->getPoint(v1,0,pt1);
  double dx = pt0[0] - pt1[0];
  double dy = pt0[1] - pt1[1];
  double len = sqrt(dx*dx + dy*dy);

  // Calculate the center of domain to evalaute direction of normal vector
  double xMin = m3dc1_model::instance()->boundingBox[0];
  double yMin = m3dc1_model::instance()->boundingBox[1];
  double xMax = m3dc1_model::instance()->boundingBox[2];
  double yMax = m3dc1_model::instance()->boundingBox[3];
  double center[2] = {(xMin+xMax)/2, (yMin+yMax)/2};
  
  
  // Find the vector from the vertex v to center point
  double vec[2] = {center[0] - pt0[0], center[1] - pt0[1]}; 
  
  // Find the normal vector (temporaray since direction needs to be adjusted)
  double normTemp[2];
  normTemp[0] = dy/len;
  normTemp[1] = -1.0*dx/len;

  // Find the direction of normal vector (inward or outward)
  double dir = normTemp[0]*vec[0] + normTemp[1]*vec[1];

  // Find the final normal vector
  if (dir > 0) 
  {
    normalVec[0] = -normTemp[0];
    normalVec[1] = -normTemp[1];
  }
  else
  {
    normalVec[0] = normTemp[0];
    normalVec[1] = normTemp[1];	
  }
}



#ifdef M3DC1_TRILINOS
#include <Epetra_MultiVector.h>
#include <AztecOO.h>
#include <Epetra_Version.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include "m3dc1_ls.h"

void verifyFieldEpetraVector(apf::Field* f, Epetra_MultiVector* x)
{
  double* field_data =getArrayData(f);
  assert(countComponents(f)*m3dc1_mesh::instance()->mesh->count(0)==x->MyLength());

  for (int i=0; i<x->MyLength(); ++i)
  { 
    assert(!value_is_nan((*x)[0][i]) && !value_is_nan(field_data[i]));

    if (!(m3dc1_double_isequal((*x)[0][i], field_data[i])))
      std::cout<<"[p"<<PCU_Comm_Self()<<"] x["<<i<<"]="<<(*x)[0][i]
                <<", field_data["<<i<<"]="<<field_data[i]<<"\n";
      assert(m3dc1_double_isequal((*x)[0][i], field_data[i]));
  }
}
#endif

int m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }
  m3dc1_ls::instance()->add_matrix(*matrix_id, new m3dc1_epetra(*matrix_id, *matrix_type, *scalar_type, *field_id));
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_delete(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  typedef std::map<int, m3dc1_epetra*> matrix_container_map;
  m3dc1_ls::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  mat->destroy();
  delete mat;
  return M3DC1_SUCCESS;
#endif
}

void m3dc1_epetra_reset(int* matrix_id)
{
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}


int m3dc1_epetra_insert(int* matrix_id, int* row, int* col, int* scalar_type, double* val)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  assert(*scalar_type==M3DC1_REAL);

  int err = mat->epetra_mat->ReplaceGlobalValues(*row, 1, val, col);
  if (err) {
    err =mat->epetra_mat->InsertGlobalValues(*row, 1, val, col);
    assert(err == 0);
  }
  return M3DC1_SUCCESS;
#endif
}

void print_elem (int elem_id)
{
  int ielm_dim=m3dc1_mesh::instance()->mesh->getDimension();
  apf::Mesh2* m=m3dc1_mesh::instance()->mesh;
  
  int num_node_per_element;
  apf::Downward downward;

  apf::MeshEntity* e = getMdsEntity(m, ielm_dim, elem_id);
  num_node_per_element = m->getDownward(e, 0, downward);
 
  int *id = new int[num_node_per_element];
  for (int i=0; i<num_node_per_element; ++i)
    id[i] = get_ent_globalid(m, downward[i]);

  switch (num_node_per_element)
  { 
    case 3: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<"\n";
            break;
    case 4: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<"\n";
            break;
    case 5: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<" "<<id[4]<<"\n";
            break;
    case 6: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<" "<<id[4]<<" "<<id[5]<<"\n";
            break;
    default: break;
  }
  delete [] id;
}

#ifdef M3DC1_TRILINOS
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  double val[1];
  int col[1];
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    for (int j=0; j<csize; j++)
    {
      col[0] = columns[j];
      val[0] = values[i*csize+j];
      int ierr = mat->SumIntoGlobalValues(rows[i], 1, val, col);
      if (ierr) 
        ierr =mat->InsertGlobalValues(rows[i], 1, val, col);
      assert(!ierr);
    }
  } // for i
}

// seol -- this does weird thing so shouldn't be used
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values_wrong(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    int ierr = mat->SumIntoGlobalValues(rows[i], csize, &values[i*csize], columns);
    if (ierr) 
      ierr =mat->InsertGlobalValues(rows[i], csize, &values[i*csize], columns);
    assert(!ierr);
  } // for i
}

#endif


int m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (!mat)
    return M3DC1_FAILURE;

  int field = mat->get_field_id();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = m3dc1_mesh::instance()->mesh->getDimension();
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;
  int start_global_dof,end_global_dof_id;
  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowVarIdx<numVar && *columnVarIdx<numVar);
  int* rows = new int[dofPerVar*nodes_per_element];
  int* columns = new int[dofPerVar*nodes_per_element];

  if (mat->matrix_type==M3DC1_MULTIPLY)
  {
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for (int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    //FIXME: mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
    epetra_add_values(mat->epetra_mat, dofPerVar*nodes_per_element, 
                      rows,dofPerVar*nodes_per_element, columns, values);     
  }
  else //M3DC1_SOLVE
  {
    int nodeOwner[6];
    int columns_bloc[6], rows_bloc[6];
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowVarIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnVarIdx;
      for (int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for (int inode=0; inode<nodes_per_element; inode++)
    {
      // FIXME: smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      epetra_add_values(mat->epetra_mat, dofPerVar, rows+dofPerVar*inode, 
                       dofPerVar*nodes_per_element, columns, values+offset);
      offset += numValuesNode;
    }
  }
  delete [] rows;
  delete [] columns;

  return M3DC1_SUCCESS;
#endif
}


int m3dc1_epetra_setbc(int* matrix_id, int* row)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  global_ordinal_type col[1]; col[0] = row_g;
  double val[1]; val[0]=1.0; 
 
  // MatSetValue(*A, row, row, 1.0, ADD_VALUES);

  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, 1, val, col);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, 1, val, col);
  assert(err == 0);

  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  std::vector <global_ordinal_type> columns_g(*numVals);
  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  for (int i=0; i<*numVals; i++)
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;
//  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, *numVals, values, &columns_g[0]);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, *numVals, values, &columns_g[0]);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_print(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),/*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();
  A.Print(cout);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_write(int* matrix_id, const char* filename, int* skip_zero, int* start_index)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (!filename)
    return m3dc1_epetra_print(matrix_id);

  char matrix_filename[256];
  sprintf(matrix_filename,"%s-%d",filename, PCU_Comm_Self());
  if (*skip_zero==0)
    write_matrix(mat->epetra_mat, matrix_filename,false,*start_index);
  else
    write_matrix(mat->epetra_mat, matrix_filename,true,*start_index);

  // assemble matrix
  Epetra_Export exporter(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

  sprintf(matrix_filename,"assembled-%s-%d",filename, PCU_Comm_Self());

  if (*skip_zero==0)
    write_matrix(&A, matrix_filename,false,*start_index);
  else
    write_matrix(&A, matrix_filename,true,*start_index);

  return M3DC1_SUCCESS;
#endif
}

#ifdef M3DC1_TRILINOS
void copyEpetraVec2Field(Epetra_MultiVector* x, apf::Field* f)
{
  int start_global_dofid, num_dof = countComponents(f);
  std::vector<double> dof_data(num_dof);
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;

  int index=0;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (get_ent_ownpartid(m, e)!=PCU_Comm_Self()) continue;
    for (int i=0; i<num_dof; ++i)
      dof_data.at(i)=(*x)[0][index++];
    setComponents(f, e, 0, &(dof_data[0]));
  }
  m->end(it);
  assert(index == num_dof*m3dc1_mesh::instance()->num_own_ent[0]);
  synchronize_field(f);
}
#endif

int m3dc1_solver_aztec_old(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, int* num_iter, 
  double* tolerance,const char* krylov_solver, const char* preconditioner)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  else
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<* x_fieldid<<" (tol "<<*tolerance<<")\n";

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),
			 /*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();

//  const char* mfilename = mfilename_str.c_str();
//  EpetraExt::RowMatrixToMatlabFile(mfilename, A);

  // copy field to vec  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  synchronize_field(b_field);
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1);
  b.Export(b_field_vec,exporter,Insert);

  // vector for solution
  Epetra_MultiVector x(*(mat->_owned_map), 1);

  Epetra_LinearProblem problem(&A,&x,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_output,1);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  int overlap = 1;
  solver.SetAztecOption(AZ_overlap, overlap);

  solver.Iterate(*num_iter,*tolerance);
  mat->num_solver_iter = solver.NumIters();
  
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  copyEpetraVec2Field(&x, x_field);
  return M3DC1_SUCCESS;
#endif
}


int m3dc1_solver_aztec_unused(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, int* num_iter, double* tolerance)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  else
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<* x_fieldid<<" (tol "<<*tolerance<<")\n";

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),
			 /*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

//  const char* mfilename = mfilename_str.c_str();
//  EpetraExt::RowMatrixToMatlabFile(mfilename, A);

  // copy field to vec  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  synchronize_field(b_field);
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1);
  b.Export(b_field_vec,exporter,Insert);

  // vector for solution
  Epetra_MultiVector x(*(mat->_owned_map), 1);

  Epetra_LinearProblem problem(&A,&x,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_output,1);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  int overlap = 1;
  solver.SetAztecOption(AZ_overlap, overlap);

  solver.Iterate(*num_iter,*tolerance);
  mat->num_solver_iter = solver.NumIters();
  
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  copyEpetraVec2Field(&x, x_field);
  return M3DC1_SUCCESS;
#endif
}

bool isNotAlnum(char c) {
    return isalnum(c) == 0;
}


// For filtering spaces and control characters in Trilinos
// options
bool invalidChar (char c) 
{  
  return !((c > 65 && c < 90) ||
	   (c > 97 && c < 122) ||
	   (c == '_'));
} 


int m3dc1_solver_aztec(int* matrix_id, FieldID* x_fieldid, FieldID*
		       b_fieldid, int* num_iter, double* tolerance,
		       const char* krylov_solver, const char*
		       preconditioner, const char* sub_dom_solver,
		       int* overlap, int* graph_fill, double*
		       ilu_drop_tol, double* ilu_fill, double*
		       ilu_omega, int* poly_ord)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  else
    if (!PCU_Comm_Self())
	std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<*
	matrix_id<<", field "<<* x_fieldid<<" (tol "<<*tolerance<<")\n";

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),
			 /*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

  // copy field to vec  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  synchronize_field(b_field);
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1);
  b.Export(b_field_vec,exporter,Insert);

  // vector for solution
  Epetra_MultiVector x(*(mat->_owned_map), 1);

  Epetra_LinearProblem problem(&A,&x,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_output,1);

  // Setup solver from input/default

  // Convert const char* to string for comparison
  std::string krylov_solver_s = krylov_solver;
  std::string preconditioner_s = preconditioner;
  std::string sub_dom_solver_s = sub_dom_solver;

  krylov_solver_s.erase(std::remove_if (krylov_solver_s.begin(),
				       krylov_solver_s.end(),
				       invalidChar),
			krylov_solver_s.end());
  preconditioner_s.erase(std::remove_if (preconditioner_s.begin(),
					preconditioner_s.end(),
					invalidChar),
			 preconditioner_s.end());
  sub_dom_solver_s.erase(std::remove_if (sub_dom_solver_s.begin(),
					sub_dom_solver_s.end(),
					invalidChar),
			 sub_dom_solver_s.end());
  
  if (krylov_solver_s == "cg")
    solver.SetAztecOption(AZ_solver, AZ_cg);

  if (krylov_solver_s == "cg_condnum")
    solver.SetAztecOption(AZ_solver, AZ_cg_condnum);

  if (krylov_solver_s == "gmres")
    solver.SetAztecOption(AZ_solver, AZ_gmres);

  if (krylov_solver_s == "gmres_condnum")
    solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);

  if (krylov_solver_s == "cgs")
    solver.SetAztecOption(AZ_solver, AZ_cgs);

  if (krylov_solver_s == "tfqmr")
    solver.SetAztecOption(AZ_solver, AZ_tfqmr);

  // Setup preconditioner from input/default
  if (preconditioner_s == "none")
    solver.SetAztecOption(AZ_precond, AZ_none);

  if (preconditioner_s == "Jacobi")
    solver.SetAztecOption(AZ_precond, AZ_Jacobi);

  if (preconditioner_s == "Neumann")
    solver.SetAztecOption(AZ_precond, AZ_Neumann);

  if (preconditioner_s == "ls")
    solver.SetAztecOption(AZ_precond, AZ_ls);

  if (preconditioner_s == "sym_GS")
    solver.SetAztecOption(AZ_precond, AZ_sym_GS);

  if (preconditioner_s == "dom_decomp")
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  
  // Setup subdomain solver from input/default
  if (preconditioner_s == "dom_decomp")
    {
      if (sub_dom_solver_s == "ilu")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);

      if (sub_dom_solver_s == "lu")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_lu);

      if (sub_dom_solver_s == "ilut")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);

      if (sub_dom_solver_s == "rilu")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_rilu);

      if (sub_dom_solver_s == "bilu")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_bilu);

      if (sub_dom_solver_s == "icc")
	solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);
      
      // Set Aztec options from input for dom_decomp
      solver.SetAztecOption(AZ_overlap, *overlap);
      solver.SetAztecOption(AZ_graph_fill, *graph_fill);

      // Setup Aztec parameters from input/default
      solver.SetAztecParam(AZ_tol, *tolerance);
      solver.SetAztecParam(AZ_drop, *ilu_drop_tol);
      solver.SetAztecParam(AZ_ilut_fill, *ilu_fill);
      if (sub_dom_solver_s == "rilu")
	solver.SetAztecParam(AZ_omega, *ilu_omega);
    }
  
  // Setup alternate preconditioner options from input/default
  if (preconditioner_s == "Jacobi" ||
      preconditioner_s == "Neumann" ||
      preconditioner_s == "ls" ||
      preconditioner_s == "sym_GS")
    solver.SetAztecOption(AZ_poly_ord, *poly_ord);

  // Now perform the solve
  solver.Iterate(*num_iter,*tolerance);
  mat->num_solver_iter = solver.NumIters();
  
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  copyEpetraVec2Field(&x, x_field);
  return M3DC1_SUCCESS;
#endif
}

// solve Ax=b
#ifdef M3DC1_TRILINOS
//#include "Amesos2.hpp"
//#include "Amesos2_Version.hpp"
#endif

int m3dc1_solver_amesos(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, const char* solver_name)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: check if amesos2 is available\n";
  return M3DC1_FAILURE;

  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  Epetra_Export exporter(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  // copy field to vec
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1 );

  A.Export(*(mat->epetra_mat),exporter,Add);
  b.Export(b_field_vec,exporter,Insert);

  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

// using SuperLUDIST

  // Before we do anything, check that the solver is enabled
/*
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if ( !Amesos2::query(solver_name) )
  {
    if (!PCU_Comm_Self()) 
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": "<< solver_name << " not enabled.  Exiting..." << std::endl;
    return M3DC1_FAILURE;	// Otherwise CTest will pick it up as failure, which it isn't really
  }
  else
    if (!PCU_Comm_Self())  *fos <<__func__<<": "<<Amesos2::version() <<" with "<<solver_name<< std::endl;
  // Constructor from Factory
  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;

  Teuchos::RCP<MAT> rcp_A = Teuchos::rcp(&A);
//  A->Print(*fos);
//    *fos <<"["<<PCU_Comm_Self()<<"] GlobalNumNonZero="<<rcp_A->NumGlobalNonzeros()<<"\n";
//  int nrows = A->NumGlobalEntries();
//  if (!PCU_Comm_Self()) 
//      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": nrows="<<nrows<<std::endl;

    int numVecs=1;
    Epetra_MultiVector x(*(mat->_owned_map), 1 );
    Teuchos::RCP<MV> X = Teuchos::rcp(&x);
    Teuchos::RCP<MV> B = Teuchos::rcp(&b);

    // copy field to vec
    for (int i=0; i<B->MyLength(); ++i)
      (*B)[0][i] = b_field_data[i];

    // Solve A*Xhat = B for Xhat using the Superlu solver
    Teuchos::RCP<Amesos2::Solver<MAT,MV> > solver;
    try 
    {
      solver = Amesos2::create<MAT,MV>(solver_name, rcp_A, X, B );
    }
    catch (std::invalid_argument e)
    {
      if (!PCU_Comm_Self()) *fos <<"[M3D-C1 ERROR] "<<__func__<<": "<< e.what() << std::endl;
      return M3DC1_FAILURE;
    }

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    //solver->printTiming(*fos);
    //X.Describe(*fos, Teuchos::VERB_EXTREME);

  //  X->Print(*fos); //, Teuchos::VERB_EXTREME);
  // print residual
//  double norm_data[1];
//  x.Norm2(norm_data);
//  *fos << "["<<PCU_Comm_Self()<<"] Norm2 of Ax - b = " << norm_data[0] << std::endl;

  // get solution
  Epetra_Import importer(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_MultiVector sol_x(*(mat->_overlap_map),1);
  sol_x.Import(x, importer, Add);
 
  double** s;
  sol_x.ExtractView(&s);
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  double* x_field_data = getArrayData(x_field);
  for (int i=0; i<sol_x.MyLength(); ++i)
    x_field_data[i] =s[0][i];
*/
  return M3DC1_SUCCESS;

#endif
}

int m3dc1_solver_getnumiter(int* matrix_id, int * num_iter)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  *num_iter = mat->num_solver_iter;
  return M3DC1_SUCCESS;
#endif

}

// local matrix multiplication
// do accumulate(out_field) for global result
int m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_MULTIPLY)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  if (!mat->epetra_mat->Filled())
    mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());

  apf::Field* in_field = (*(m3dc1_mesh::instance()->field_container))[*in_fieldid]->get_field();
  double* in_field_data =getArrayData(in_field);

  Epetra_MultiVector x(mat->epetra_mat->RowMap(), 1);
  // copy field to vec
  for (int i=0; i<x.MyLength(); ++i)
    x[0][i] = in_field_data[i];
  Epetra_MultiVector b(mat->epetra_mat->RowMap(), 1);
  EPETRA_CHK_ERR(mat->epetra_mat->Multiply(false, x, b));
  apf::Field* out_field = (*(m3dc1_mesh::instance()->field_container))[*out_fieldid]->get_field();
  double* out_field_data =getArrayData(out_field);
  b.ExtractCopy(out_field_data, b.MyLength());
  accumulate(out_field);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_assemble(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());
  return M3DC1_SUCCESS;
#endif
}

// for solver debugging

//*******************************************************
int m3dc1_matrix_getstatus (int* matrix_id, int* status) // checkMatrixStatus_
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  *status = mat->get_status();
  return M3DC1_SUCCESS;
}

// getMatrixLocalDofNum_
//*******************************************************
void m3dc1_matrix_getlocalnumdof(int* matrix_id, int *num_own_dof) 
//*******************************************************
{
/*
  int NumberingId;
  Matrix *mat = getMatrix(*matrix_id);
  NumberingId = mat->getNumberingid();
  NumberingContainer * ncptr = getNumberingContainer(NumberingId);
  // this is equivalent to mat_dim
  *ldb = ncptr->proclastdofplusone - ncptr->procfirstdof;  
*/
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);

  int num_own_ent=m3dc1_mesh::instance()->num_own_ent[0];
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[mat->fieldOrdering];
  *num_own_dof = (m3dc1_mesh::instance()->num_own_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
}

// getMatrixGlobalDofs_
//*******************************************************
void m3dc1_matrix_getglobalnumdof(int *matrix_id, int *num_global_dof)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  int NumberingId = mat->getNumberingid();
  NumberingContainer * ncptr = getNumberingContainer(NumberingId);
  *numglobaldofs =  ncptr->numglobaldofs;
*/
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);

  int num_own_ent=m3dc1_mesh::instance()->num_global_ent[0];
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[mat->fieldOrdering];
  *num_global_dof = (m3dc1_mesh::instance()->num_global_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
}

// getMatrixPetscDnnzOnnz_
//*******************************************************
void m3dc1_matrix_getpetscdnnzonnz(int *matrix_id, int *valType, int *d_nnz, int *o_nnz)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  solveMatrix *smat = (solveMatrix*) mat;
  // a member function of the solve matrix 
  smat->assembleDnnzOnnz(valType, d_nnz, o_nnz); 
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// getMatrixNNZRowSize_
//*******************************************************
void gm3dc1_matrix_getnnzrowsize(int *matrix_id, int *valType, int *rowSize)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  solveMatrix *smat = (solveMatrix*) mat;
  
  *rowSize = smat->getRowSize(*valType);
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// getMatrixNNZRowId_
//*******************************************************
void m3dc1_matrix_getnnzrowid(int *matrix_id, int *valType, int* ith, int *rowId)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  solveMatrix *smat = (solveMatrix*) mat;
  *rowId = smat->getRowId(*valType, *ith)+1;
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// getMatrixNNZColSize_
//*******************************************************
void  m3dc1_matrix_getnnzcolsize(int *matrix_id, int *valType, int *rowId, int *colSize)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  solveMatrix *smat = (solveMatrix*) mat;
  *colSize = smat->getColSize(*valType, *rowId-1);
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// getMatrixNNZValues_
//*******************************************************
void m3dc1_matrix_getnnzvalue(int *matrix_id, int *valType, int *rowId, int *colId, double *dvalues)
//*******************************************************
{
/*
  Matrix *mat =  getMatrix(*matrix_id);
  solveMatrix *smat = (solveMatrix*) mat;
  smat->getNNZValues(*valType, *rowId-1, colId, dvalues);
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// getMatrixFirstDof_
//*******************************************************
void m3dc1_matrix_getmatrixfirstdof(int *matrix_id, int *firstdof)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrix_id);
  int NumberingId = mat->getNumberingid();

  NumberingContainer * ncptr = getNumberingContainer(NumberingId);
  *firstdof = ncptr->procfirstdof+1;
*/
  // FIXME: to be provided
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";
}

// cleanMatrixValues_ is identical to m3dc1_matrix_reset

// setMatrixSoln_
//*******************************************************
void m3dc1_matrix_setsoln(int *matrix_id, int *valType, double *soln)
//*******************************************************
{
/*
  Matrix *mat = getMatrix(*matrixid);
  int NumberingId = mat->getNumberingid();
  NumberingContainer * ncptr = getNumberingContainer(NumberingId);
  ncptr->shareInfo(soln, 1+(*valType));
  mat->setStatus(4);
*/
  if (!PCU_Comm_Self())
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported\n";

}

#define m3dc1_matrix_setsoln setMatrixSoln_
