/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_MODEL_H
#define M3DC1_MODEL_H
#include "apf.h"
#include "apfMesh.h"
#include "gmi.h"
#include "gmi_mesh.h"
#include <mpi.h>
#include <map>
#include <vector>

void get_gent_adj(int gent_dim, int gent_id, int adj_dim, std::vector<int>& adj_ids);

int get_prev_plane_partid(int plane_id);
int get_next_plane_partid(int plane_id);
const double tol=1e-4;
void faceFunction(double const p[2], double x[3], void * data);
void vertexFunction(double const p[2], double x[3], void * data);

inline int checkSamePoint (double a[3],double b[3])
{
  double res=fabs(a[0]-b[0])+fabs(a[1]-b[1])+fabs(a[2]-b[2]);
  return res<tol;
}
inline int checkSamePoint2D (double* a,double* b)
{
  double res=fabs(a[0]-b[0])+fabs(a[1]-b[1]);
  return res<tol;
}
inline double getDist2D(double * cd1, double* cd2)
{
  double dist=0;
  dist=(cd1[0]-cd2[0])*(cd1[0]-cd2[0]);
  dist+=(cd1[1]-cd2[1])*(cd1[1]-cd2[1]);
  return sqrt(dist);
}

// for tokamak geometry with general components
enum curveType{LIN=1,POLY,BSPLINE, OTHER};
void create_loop ( int* loopId, int * numEdges, int* EdgeList);
void create_edge ( int* edge, int* startingvtx, int* endingvtx);
void create_vtx ( int* vertex, double* position);

// used by test/convert_polar
gmi_ent* create_model_vertex( int id, double* position);
gmi_ent* create_b_spline_curve( int id, int order, int numPts, double* ctrlPts, double * knots, double* weight);

void attach_linear_curve ( int* edge);
void attach_polynomial_curve ( int* edge, int* order, double* coefficients );
void attach_cubic_hermite_curve ( int* edge, double * positions, double * tangents);
void attach_natural_cubic_curve ( int* edge, int * numPts, double * points);
void attach_periodic_cubic_curve ( int* edge, int * numPts, double * points);
void attach_piecewise_linear_curve ( int* edge, int * numPts, double * points);
void attach_b_spline_curve ( int* edge, int * order, int* numPts, double* ctrlPts, double * knots, double* weight);
void eval_position ( int* edge, double* para, double* position );
void eval_normal ( int* edge, double* para, double* normal );
void eval_curvature ( int* edge, double* para, double* curvature );

void load_model(const char* filename);
void offset_point (double * src, double * normal, double * thickness, double * des);

class m3dc1_model
{
public:
// functions
  m3dc1_model();
  ~m3dc1_model();
  static m3dc1_model* instance();

  void load_analytic_model(const char *name);

  void create3D(); // create 3D region

  void set_num_plane(int);
  void set_phi(double min_val, double max_val);
  double get_phi(int part_id);
  void set_phi(int plane_id, double phi);
  void get_phi(int plane_id, double* phi);
  gmi_ent* geomEntNextPlane(gmi_ent* gent); 
  gmi_ent* geomEntBtwPlane(gmi_ent* gent); 

// data
  gmi_model* model;
 
  int xperiodic, yperiodic; // if periodic in x/y direction

  MPI_Comm oldComm;

  int group_size; // # processes per group
  int num_plane; // # total planes
  int local_planeid; // the local plane id
  int prev_plane_partid; // id of corresponding part in the prev. plane
  int next_plane_partid; // id of corresponding part in the next plane
  bool snapping; // support for snapping

  int modelType = 1;  // = 1 for analytical, = 2 for PUMI (.dmg) model
  double* phi;
  int numEntOrig[3];
  double boundingBox[4];
  void caculateBoundingBox();

// for m3dc1_field_sum_plane
  void setupCommGroupsPlane();
  MPI_Comm & getMPICommPlane();
  std::map<int,MPI_Comm> PlaneGroups;
// storage of 3D model tags 
  int** ge_tag;
  void save_gtag();
  void restore_gtag();
private:
  static m3dc1_model* _instance;
  std::map<gmi_ent*, std::pair<gmi_ent*,gmi_ent*> > newModelEnts;
};
#endif
