/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_model.h"
#include "m3dc1_scorec.h"
#include "PCU.h"
#include "m3dc1_mesh.h"
#include "apfMDS.h"
#include "gmi_analytic.h"
#include "gmi_base.h"
#include "Expression.h"
#include "CMODGeoExpression.h"
#include "PolyNomial.h"
#include "BSpline.h"
#include <utility>
#include <string.h>
#include <assert.h>
#include <iostream>

int separatrixLoop=-1, innerWallLoop=-1, outerWallLoop=-1, vacuumLoop=-1;
std::map< int, std::vector<int> > loopContainer;
std::map<int, std::pair<int, int> > edgeContainer;
std::map<int, int> edgeType;
std::map<int, std::vector<double> > vtxContainer;
std::vector<void*> data2Clean;

void get_gent_adj(int gent_dim, int gent_id, int adj_dim, vector<int>& adj_ids)
{ 
  switch (gent_dim)
  { 
    case 1: // model edge
          adj_ids.push_back(edgeContainer[gent_id].first);
          adj_ids.push_back(edgeContainer[gent_id].second);
          break;
    case 2: // model face
           {
             std::map< int, std::vector<int> >::iterator it=loopContainer.begin();
             for (int i=1; i<gent_id; ++i)
               ++it;
             if (adj_dim==1)
             { 
               for (unsigned int j=0; j<it->second.size(); ++j)
                 adj_ids.push_back((it->second)[j]);
             }
             else // adj_dim==0
             { 
               int edge;
               std::set<int> vids;
               for (unsigned int j=0; j<it->second.size(); ++j)
               { 
                 edge=(it->second)[j];
                 vids.insert(edgeContainer[edge].first);
                 vids.insert(edgeContainer[edge].second);
               }
               for (std::set<int>::iterator sit=vids.begin(); sit!=vids.end(); ++sit)
                 adj_ids.push_back(*sit);
             }
             break;
           } 
    default: // model vtx or region
            std::cout<<__func__<<" FATAL ERROR: model adjacency not supported\n";
            break;
  }
}

void interpolateCubicBSpline( vector<double>& points,vector<double>& knots, vector<double> &ctrlPoints, int bc);
void faceFunction(double const p[2], double x[3], void * data) {}
void vertexFunction(double const p[2], double x[3], void * data) {}

// **********************************************
int get_prev_plane_partid(int partid)
// **********************************************
{
  int prev_plane_partid = (partid-m3dc1_model::instance()->group_size)%PCU_Comm_Peers();
  if (prev_plane_partid<0)
    return prev_plane_partid+PCU_Comm_Peers();
  else 
    return prev_plane_partid;
}

// **********************************************
int get_next_plane_partid(int partid)
// **********************************************
{
  return (partid+m3dc1_model::instance()->group_size)%PCU_Comm_Peers();
}


// m3dc1_model
m3dc1_model::m3dc1_model()
{
  model=NULL;
  phi=NULL;
  xperiodic=yperiodic=0;
  local_planeid=0;
  num_plane=1;
  snapping=false;
  group_size = PCU_Comm_Peers();
  oldComm = PCU_Get_Comm();
  ge_tag=NULL;
}

m3dc1_model::~m3dc1_model()
{
  delete [] phi;
  for(int i=0; i<data2Clean.size(); i++)
  {
    M3DC1::Expression** ptr=(M3DC1::Expression**)data2Clean.at(i);
    delete ptr[0];
    delete ptr[1];
    delete []ptr;
  }
  
  if (ge_tag)
  {
    for(int i = 0; i < 3; ++i)
      delete [] ge_tag[i];
    delete [] ge_tag;
  }
//  PUMI_Geom_Del(model);
}

m3dc1_model* m3dc1_model::_instance=NULL;
m3dc1_model* m3dc1_model::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_model();
  return _instance;
}

void edgeFunction(double const p[2], double *xyz, void*  data)
{
  M3DC1::evalCoord(p[0], xyz, data);
  xyz[2]=0.;
}

// **********************************************
void reparam_zero(double const from[2], double to[2], void*)
// **********************************************
{
  to[0] = 0;
  to[1] = 0;
}

// **********************************************
void reparam_one(double const from[2], double to[2], void*)
// **********************************************
{
  to[0] = 1;
  to[1] = 0;
}

// **********************************************
agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
// **********************************************
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

// **********************************************
agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
// **********************************************
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

// **********************************************
void make_edge_topo(gmi_model* m, gmi_ent* e, int v0tag, int v1tag)
// **********************************************
{
  agm_bdry b = add_bdry(m, e);
  agm_use u0 = add_adj(m, b, v0tag);
  gmi_add_analytic_reparam(m, u0, reparam_zero, 0);
  agm_use u1 = add_adj(m, b, v1tag);
  gmi_add_analytic_reparam(m, u1, reparam_one, 0);
}

// **********************************************
void make_face_topo (gmi_ent* face, std::vector <int> tagsOfEdgesOnFace)
// **********************************************
{
  agm_bdry b = add_bdry(m3dc1_model::instance()->model, face);
  for (int i = 0; i < tagsOfEdgesOnFace.size(); ++i)
  {
    agm_use u1 = add_adj(m3dc1_model::instance()->model, b, tagsOfEdgesOnFace[i]);
    gmi_add_analytic_reparam(m3dc1_model::instance()->model, u1, faceFunction, 0);
  }
}

// **********************************************
void m3dc1_model::load_analytic_model(const char *name)
// **********************************************
{
  std::string filename(name);
  std::string analyticmodel("AnalyticModel");
  if (!filename.compare(analyticmodel))
  {
    FILE *fp = fopen(name, "r");
    assert(fp);
    double para[5];
    for(int i=0; i<5; i++)
      fscanf(fp, "%lf", para+i);  
    fclose(fp);
  
    M3DC1::Expression** data=new M3DC1::Expression*[2];
    data[0] = new M3DC1::CMODExpressionR(para[0],para[1],para[2]);
    data[1] = new M3DC1::CMODExpressionZ(para[3],para[4]);
    data2Clean.push_back(data);
    int edgePeriodic = 1;
    double edgeRange[2] = {0, 2*M3DC1_PI};
    gmi_add_analytic(model, 1, model->n[1]+1, edgeFunction, &edgePeriodic, &edgeRange, data);

    int facePeriodic[2] = {0, 0};
    double faceRanges[2][2] = {{0,0},{0,0}};
    gmi_add_analytic(model, 2, model->n[2]+1, faceFunction, facePeriodic, faceRanges, NULL);
  }
  else 
    load_model(name);

  return;
}

// **********************************************
void load_model(const char* filename)
// **********************************************
{
  FILE* fp= fopen(filename, "r");
  int numL,separatrixLoop, innerWallLoop, outerWallLoop, vacuumLoop;
  fscanf(fp,"%d %d %d %d %d\n", &numL, &separatrixLoop, &innerWallLoop, &outerWallLoop, &vacuumLoop);

  // Identify outermost loop. We will always have inner loop even if numL==1.
  int outerMostLoop = -1;
  if (outerWallLoop > 0 || vacuumLoop > 0)
  {
    if (vacuumLoop > 0)
      outerMostLoop = vacuumLoop;
    else
      outerMostLoop = outerWallLoop;
  }

  std::vector<int> loop_ids;
  loop_ids.resize(numL);

  for (int i=0; i< numL; i++)
  {
    int numE;
    int loop; // loop ID
    fscanf(fp,"%d %d\n", &loop, &numE);
    loop_ids[i] = loop;
    // first read all vtx on the loop
    for( int j=0; j<numE; j++)
    {
      double xyz[3];
      int vertex;
      fscanf(fp,"%d %lf %lf %lf\n", &vertex, xyz, xyz+1, xyz+2);
      create_vtx(&vertex,xyz);
    }
    int * edges= new int[numE];
    for(int i=0; i<numE; i++)
    {
      int edge, beginvtx, endvtx,edgeType;
      fscanf(fp,"%d %d %d %d\n", &edge, &beginvtx, &endvtx, &edgeType);
      if (loop == innerWallLoop)
        m3dc1_model::instance()->innerLoop.push_back(edge);
      if (loop == outerMostLoop)
        m3dc1_model::instance()->outerLoop.push_back(edge);

      if (edgeContainer.find(edge)!=edgeContainer.end()) // edge already created
      {
        edges[i]=edge;
        switch (edgeType)
        {
          case BSPLINE:
               {
                 int order, numPts;
                 double dummy1, dummy2;
                 fscanf(fp,"%d %d ", &order, &numPts);
                 for (int k=0; k < order+numPts; k++)
                   fscanf(fp,"%lf ", &dummy1);
                 for (int k=0; k<numPts; k++)
                   fscanf(fp,"%lf %lf ", &dummy1, &dummy2);
                 fscanf(fp, "\n");
               }
          default: break;
        }
        continue;
      }

      create_edge(&edge,&beginvtx, &endvtx);
      edges[i]=edge;
      switch( edgeType )
      {
        case LIN:
          attach_linear_curve(&edge);
          break;
        case POLY:
        case BSPLINE:
        {
          int order;
          int numPts;
          fscanf(fp,"%d %d ", &order, &numPts);
          double * knots= new double[order+numPts];
          double * ctrlPts= new double[2*numPts];
          for( int k=0; k<order+numPts; k++)
          {
            fscanf(fp,"%lf ",knots+k);
          }
          for( int k=0; k<numPts; k++)
          {
            fscanf(fp,"%lf %lf ",ctrlPts+2*k,ctrlPts+2*k+1);
          }
          fscanf(fp, "\n");
          attach_b_spline_curve(&edge, &order, &numPts,ctrlPts, knots, NULL);
          delete []knots;
          delete []ctrlPts;
          break;
        }
        default: std::cout<<"[M3DC1 ERROR] "<<__func__<<": unsupported curve type "<<std::endl; 
                 throw 1; 
      }
    }
    // loop doesn't exist
    assert (loopContainer.find(loop)==loopContainer.end());
    
    for( int i=0; i<numE; i++)
      loopContainer[loop].push_back(edges[i]);
    delete [] edges;
  }

  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};

  // Delete this block of code in first if statement when snapping is set to default
  // for every case
  int numFaces = 0;
  fscanf(fp,"%d\n", &numFaces);
  
  //if (!PCU_Comm_Self()) cout<<"numFaces in model file "<<numFaces<<"\n";

  if (!numFaces) // Faces actually exist but not provided in output model file (.txt file)
  {
    for (int i=1; i<=numL; ++i)
      gmi_ent* gf=gmi_add_analytic(m3dc1_model::instance()->model, 2, i,
                     faceFunction, facePeriodic, faceRanges, NULL);
      if (!PCU_Comm_Self()) cout<<"[M3DC1_INFO] Regenerate model/mesh files using the latest mesh generation program\n";  
  }

  if (numFaces>0)  
  {
    std::vector <int> edgeIdsOnFace;
    for (int i = 1; i <= numFaces; ++i)
    {
      int faceNumber, numLoops;
      fscanf(fp,"%d %d", &faceNumber, &numLoops);
      gmi_ent* gf=gmi_add_analytic(m3dc1_model::instance()->model, 2, faceNumber,
                      faceFunction, facePeriodic, faceRanges, NULL);
      std::vector <int> edgeIdsOnLoop;
      for (int j = 0;j < numLoops; ++j)
      {
        int loopNumber;
        fscanf(fp,"%d", &loopNumber);
        for (int k = 0; k < loopContainer[loopNumber].size(); ++k)
          edgeIdsOnFace.push_back(loopContainer[loopNumber][k]);
      }
      fscanf(fp, "\n");

      make_face_topo(gf, edgeIdsOnFace);
      edgeIdsOnFace.clear();
    }
  }
  fclose(fp);

}

// **********************************************
void create_edge( int* edge, int* startingvtx, int* endingvtx)
// **********************************************
{
  if(edgeContainer.find(*edge)!=edgeContainer.end())
  {
    std::cout<<"[M3DC1 ERROR] "<<__func__<<": edge "<<*edge<<" has been created previously "<<std::endl;
    throw 1;
  }
  edgeContainer[*edge]=std::make_pair(*startingvtx, *endingvtx); 
}

// **********************************************
void create_vtx( int* vertex, double* position)
// **********************************************
{
  vector<double> vtx_coord(3);
  vtx_coord[0]=position[0];
  vtx_coord[1]=position[1];
  vtx_coord[2]=0.;
  vtxContainer[*vertex]=vtx_coord;
  gmi_add_analytic(m3dc1_model::instance()->model, 0, *vertex, vertexFunction, NULL, NULL, NULL);
}

// **********************************************
void attach_linear_curve (int* edge)
// **********************************************
{
  edgeType[*edge]=LIN;
  std::pair< int, int> vtx=edgeContainer[*edge];
  std::vector<double>& cd1=  vtxContainer[vtx.first];
  std::vector<double>& cd2 = vtxContainer[vtx.second];
  int order=2;
  double coefficients[4];
  coefficients[0]=cd2[0]-cd1[0];
  coefficients[1]=cd2[1]-cd1[1];
  coefficients[2]=cd1[0];
  coefficients[3]=cd1[1]; 
  attach_polynomial_curve ( edge, &order, coefficients );
}

// **********************************************
void attach_polynomial_curve ( int* edge, int* order, double* coefficients )
// **********************************************
{
  if( *order>2) edgeType[*edge]=POLY;
  std::pair< int, int> vtx=edgeContainer[*edge];

  std::vector<double> X_p, Y_p;
  for (int i=0; i<*order;i++)
  {
    X_p.push_back(coefficients[2*i]);
    Y_p.push_back(coefficients[2*i+1]);
  }

  M3DC1::PolyNomial** data=new M3DC1::PolyNomial*[2];
  data[0] = new M3DC1::PolyNomial(*order,X_p);
  data[1] = new M3DC1::PolyNomial(*order,Y_p);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  double edgeRange[2]= {0.0, 1.0};
//  gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  gmi_ent* gedge = gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<*edge<<" - vtx("<<vtx.first<<", "<<vtx.second<<")\n";

}
// use clamped b-spline as underlying representation
// in fact bezier curve
// **********************************************
void attach_cubic_hermite_curve ( int* edge, double * positions, double * tangents)
// **********************************************
{
  edgeType[*edge]=BSPLINE;
  double knots[]={0,0,0,0,1,1,1,1};
  std::pair< int, int> vtx=edgeContainer[*edge];
  int order_p=4;
  vector<double> ctrlPointsX(4),ctrlPointsY(4),weight;
  vector<double> knots_vec(knots,knots+8);
  ctrlPointsX.at(0)=positions[0];
  ctrlPointsY.at(0)=positions[1];
  ctrlPointsX.at(1)=positions[0]+1.0/3.0*tangents[0];
  ctrlPointsY.at(1)=positions[1]+1.0/3.0*tangents[1];
  ctrlPointsX.at(2)=positions[2]-1.0/3.0*tangents[2];
  ctrlPointsY.at(2)=positions[3]-1.0/3.0*tangents[3];
  ctrlPointsX.at(3)=positions[2];
  ctrlPointsY.at(3)=positions[3];
  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  data[0] = new M3DC1::BSpline(order_p,ctrlPointsX,knots_vec, weight);
  data[1] = new M3DC1::BSpline(order_p,ctrlPointsY,knots_vec, weight);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  double edgeRange[2] = {0.0, 1.0};
//  gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  gmi_ent* gedge = gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<*edge<<" - vtx("<<vtx.first<<", "<<vtx.second<<")\n";
}

// use clamped b-spline as underlying representation
// **********************************************
void attach_natural_cubic_curve ( int* edge, int * numPts, double * points)
// **********************************************
{
  edgeType[*edge]=BSPLINE;
  std::pair< int, int> vtx=edgeContainer[*edge];
 
  int order_p=4; 
  int knotsize=2*order_p+*numPts-2;
  vector<double> knots(knotsize,0.);
  vector<double> ctrlPointsX(*numPts+2),ctrlPointsY(*numPts+2),weight;
  for( int i=0; i<order_p; i++)
  {
    knots.at(knotsize-i-1)=1.0;
  }
  double increment=1.0/(*numPts-1);
  for (int i=0; i<*numPts-2; i++)
  {
    //double increment=inter_len.at(i)/len;
    knots.at(order_p+i)=knots.at(order_p+i-1)+increment;
  }
  vector<double> pointsX(*numPts),pointsY(*numPts);
  for( int i=0; i<*numPts; i++)
  {
    pointsX.at(i)=points[2*i];
    pointsY.at(i)=points[2*i+1];
  }
  interpolateCubicBSpline(pointsX,knots,ctrlPointsX,0);
  interpolateCubicBSpline(pointsY,knots,ctrlPointsY,0);
  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  data[0] = new M3DC1::BSpline(order_p,ctrlPointsX,knots, weight);
  data[1] = new M3DC1::BSpline(order_p,ctrlPointsY,knots, weight);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  double edgeRange[2] = {0.0, 1.0};
//  gmi_ent* ae=gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  gmi_ent* gedge=gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<*edge<<" vtx("<<vtx.first<<", "<<vtx.second<<")\n";
}

// **********************************************
void attach_piecewise_linear_curve ( int* edge, int * numPts, double * points)
// **********************************************
{
  edgeType[*edge]=BSPLINE;
  std::pair< int, int> vtx=edgeContainer[*edge];

  int order_p=2;
  int knotsize=2*order_p+*numPts-2;
  vector<double> knots(knotsize,0.);
  vector<double> ctrlPointsX(*numPts),ctrlPointsY(*numPts),weight;
  for( int i=0; i<order_p; i++)
  {
    knots.at(knotsize-i-1)=1.0;
  }
  double increment=1.0/(*numPts-1);
  for (int i=0; i<*numPts-2; i++)
  {
    //double increment=inter_len.at(i)/len;
    knots.at(order_p+i)=knots.at(order_p+i-1)+increment;
  }

  for( int i=0; i<*numPts; i++)
  {
    ctrlPointsX.at(i)=points[2*i];
    ctrlPointsY.at(i)=points[2*i+1];
  }
  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  data[0] = new M3DC1::BSpline(order_p,ctrlPointsX,knots, weight);
  data[1] = new M3DC1::BSpline(order_p,ctrlPointsY,knots, weight);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  double edgeRange[2] = {0.0, 1.0};
  //gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  gmi_ent* gedge = gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<*edge<<" vtx("<<vtx.first<<", "<<vtx.second<<")\n";
}

// **********************************************
void attach_b_spline_curve( int* edge, int * order, int* numPts, double* ctrlPts, double * knots, double* weight)
// **********************************************
{
  edgeType[*edge]=BSPLINE;
  std::pair< int, int> vtx=edgeContainer[*edge];
  vector<double> X_p, Y_p, knots_p;
  for (int i=0; i<*numPts;i++)
  {
    X_p.push_back(ctrlPts[2*i]);
    Y_p.push_back(ctrlPts[2*i+1]);
  }
  for( int i=0; i<*numPts+*order; i++)
  {
    knots_p.push_back(knots[i]);
  }

  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  vector<double> wt;
  data[0] = new M3DC1::BSpline(*order,X_p,knots_p, wt);
  data[1] = new M3DC1::BSpline(*order,Y_p,knots_p, wt);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  if(vtx.first==vtx.second) edgePeriodic=1;
  double edgeRange[2] = {0.0, 1.0};
  gmi_ent* gedge = gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  //if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<*edge<<" vtx("<<vtx.first<<", "<<vtx.second<<")\n";
}

// **********************************************
void set_gent_tag(gmi_model* model, gmi_ent* ge, int new_tag)
// **********************************************
{
  gmi_base_set_tag (model, ge, new_tag);
}

// *********************************************************
void m3dc1_model::create3D() // construct 3D model out of 2D
// *********************************************************
{
  if (newModelEnts.size()) 
  {
/*
    for (std::map<gmi_ent*, std::pair<gmi_ent*,gmi_ent*> >::iterator it=newModelEnts.begin(); it!=newModelEnts.end(); it++)
    {
      std::cout<<"[M3D-C1 INFO]::"<<__func__<<" 3D model Info"<<std::endl;
      std::cout<<" rank "<<PCU_Comm_Self()
      <<" ent_org (dim  "<<gmi_dim(model, it->first)<<", id "<<gmi_tag(model, it->first)<<") "
      <<" ent_btw (dim  "<<gmi_dim(model, it->second.first)<<", id "<<gmi_tag(model, it->second.first)<<") "
      <<" ent_next (dim  "<<gmi_dim(model, it->second.second)<<", id "<<gmi_tag(model, it->second.second)<<")\n";
    }
*/
    return;
  } 

  if (model->n[0]==0)
  {
    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO]::"<<__func__<<": model has 0 vertex, 1 edge, 1 face\n";
    assert(model->n[1]==1);
    assert(model->n[2]==1); 
    gmi_ent* ae_org = gmi_find(model, 1, 1);
    void* edgeData = gmi_analytic_data(model, ae_org);
    set_gent_tag(model, ae_org, 1+local_planeid);
    gmi_ent* af_org = gmi_find(model, 2, 1);
    set_gent_tag(model, af_org, 1+2*local_planeid);

    int edgePeriodic = 1;
    double edgeRange[2] = {0, 2*M3DC1_PI};
    int facePeriodic[2] = {0, 0};
    double faceRanges[2][2] = {{0,0},{0,0}};

    for (int i=0; i<num_plane; i++)
    {
      int planeId= (local_planeid+i)%num_plane;
      int ae_next_tag=1+(planeId+1)%num_plane;
      int af_btw_tag=2+2*planeId;
      int af_next_tag=1+2*((planeId+1)%num_plane);
      int ar_tag=1+planeId;

      gmi_ent* ae_next=NULL;
      if (i<num_plane-1)
        gmi_add_analytic(model, 1, ae_next_tag, edgeFunction, &edgePeriodic, &edgeRange, edgeData);
      ae_next=gmi_find(model, 1, ae_next_tag);
      assert(ae_next);

      gmi_ent* af_next=NULL;
      if (i<num_plane-1)
        gmi_add_analytic(model, 2, af_next_tag, faceFunction, facePeriodic, faceRanges, NULL);
      af_next=gmi_find(model, 2, af_next_tag);
      assert(af_next);

      gmi_add_analytic(model, 2, af_btw_tag, faceFunction, facePeriodic, faceRanges, NULL);
      gmi_ent* af_btw = gmi_find(model, 2, af_btw_tag); 
      assert(af_btw);

      gmi_add_analytic_region (model, ar_tag);
      gmi_ent* ar = gmi_find(model, 3, ar_tag); 
      assert(ar);

     // if (i==0)
      {
        newModelEnts[ae_org]=std::make_pair(af_btw, ae_next);
        newModelEnts[af_org]=std::make_pair(ar, af_next);
      }
      
      af_org=af_next;
      ae_org=ae_next;
    }
  }
  else
  {
    int edgePeriodic = 0;
    double edgeRange[2]= {0.0, 1.0};
    int facePeriodic[2] = {0, 0};
    double faceRanges[2][2] = {{0,0},{0,0}};
    int numL=loopContainer.size();
    int numV=model->n[0];
    int numE=model->n[1];
    int numF=model->n[2];
    assert(numL==numF);
    assert(numV==numE);
    vector<gmi_ent*> af_org_vec, ae_org_vec, av_org_vec;

    // first update the tag of the original plane
    for (int iloop=1; iloop<=numL; iloop++)
    {
      int af_org_tag= iloop+local_planeid*(numF+numE);
      gmi_ent* af_org=gmi_find(model, 2,iloop);
      af_org_vec.push_back(af_org);
      set_gent_tag(model, af_org, af_org_tag);

      for (int iedge=0; iedge<loopContainer[iloop].size(); iedge++)
      {
        int edgeTag=loopContainer[iloop].at(iedge);
        gmi_ent* ae_org = gmi_find (model, 1,edgeTag); 
        ae_org_vec.push_back(ae_org);
        assert(gmi_analytic_data(model, ae_org));
        int ae_org_tag= edgeTag + local_planeid*(numE+numV);   
        set_gent_tag(model, ae_org, ae_org_tag);
        assert(gmi_analytic_data(model, ae_org));

        gmi_ent* av_org = gmi_find (model, 0,edgeTag);
        int av_org_tag= edgeTag + local_planeid*numV;
        av_org_vec.push_back(av_org);
        set_gent_tag(model, av_org, av_org_tag);
      }
    }

    for (int i=0; i<num_plane; i++)
    {
      int planeId= (local_planeid+i)%num_plane;
      int next_plane_id= (planeId+1)%num_plane;
      for(int iloop=1; iloop<=numL; iloop++)
      {
        int af_next_tag=iloop+next_plane_id*(numF+numE);
        int ar_tag=iloop+planeId*numF;
        // first create vertex
        int numVLoop=loopContainer[iloop].size();
        for (int ivtx=0; ivtx<numVLoop; ivtx++)
        {
          int edgeTag=loopContainer[iloop].at(ivtx);
          int av_next_tag=edgeTag + next_plane_id*numV;
          if (i<num_plane-1)
          { 
            assert(!gmi_find(model,0, av_next_tag));
            gmi_add_analytic(model, 0, av_next_tag, vertexFunction, NULL, NULL, NULL);
          }
          gmi_ent* av_next=gmi_find(model,0, av_next_tag);
          assert(av_next);

          // if (i==0)
          {
            gmi_ent* gv_org = av_org_vec.at(edgeTag-1);
            newModelEnts[gv_org].second=av_next;
          }
        }
        // create edges 
        for(int iedge=0; iedge<numVLoop; iedge++)
        {
          int edgeTag=loopContainer[iloop].at(iedge);
          int ae_next_tag = edgeTag+ next_plane_id*(numV+numE);
          gmi_ent* ae_next = NULL;
          if (i<num_plane-1) 
          {
            gmi_ent* ae_org = ae_org_vec.at(edgeTag-1);
            void* data = gmi_analytic_data(model, ae_org);
            assert(data);
            assert(!gmi_find(model,1, ae_next_tag));
            gmi_ent* gedge = gmi_add_analytic(model, 1, ae_next_tag, edgeFunction, &edgePeriodic, &edgeRange, data);
            // added by seol
            std::pair<int, int> vtx=edgeContainer[edgeTag];
            make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first+ next_plane_id*numV, vtx.second+ next_plane_id*numV);
            if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] new next edge "<<ae_next_tag<<" vtx("<<vtx.first+ next_plane_id*numV<<", "<<vtx.second+ next_plane_id*numV<<")\n";

          }
          ae_next = gmi_find(model,1, ae_next_tag);
          assert(ae_next);
          
          int ae_btw_tag = edgeTag + numE + planeId*(numE+numV);
          assert(!gmi_find(model,1, ae_btw_tag));
          gmi_ent* gedge = gmi_add_analytic(model, 1, ae_btw_tag, edgeFunction, &edgePeriodic, &edgeRange, NULL);
          // added by seol
          std::pair<int, int> vtx=edgeContainer[edgeTag];
          make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first+ planeId*numV, vtx.first+ next_plane_id*numV);
          if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] new btw edge "<<ae_next_tag<<" vtx("<<vtx.first+ planeId*numV<<", "<<vtx.first+ next_plane_id*numV<<")\n";

          // if (i==0)
          {
            gmi_ent* gv_org = av_org_vec.at(edgeTag-1);
            newModelEnts[gv_org].first=gmi_find(model,1,ae_btw_tag);
          }
        } 
        // create faces between
        for(int iedge=0; iedge<numVLoop; iedge++)
        {
          int edgeTag=loopContainer[iloop].at(iedge);
          int ae_next_tag = edgeTag+ next_plane_id*(numV+numE);
          gmi_ent* ae_next = gmi_find(model,1, ae_next_tag);
          assert(ae_next);
          int af_btw_tag=edgeTag+ numF + planeId*(numE+numF);
          assert(!gmi_find(model,2, af_btw_tag));
          gmi_add_analytic(model, 2, af_btw_tag, faceFunction, facePeriodic, faceRanges, NULL);
          gmi_ent* af_btw =gmi_find(model, 2, af_btw_tag);
          assert(af_btw);
          // if (i==0)
          {
            gmi_ent* ae_org = ae_org_vec.at(edgeTag-1);
            newModelEnts[ae_org]=std::make_pair(af_btw,ae_next);
          }
        }
        // create face on next plane
        gmi_ent *af_next = NULL;
        if (i<num_plane-1)
        {
          assert(!gmi_find(model,2, af_next_tag));
          gmi_add_analytic(model, 2, af_next_tag, faceFunction, facePeriodic, faceRanges, NULL);
        }
        af_next= gmi_find(model,2, af_next_tag);
        assert(af_next);
        // create region
        assert(!gmi_find(model,3, ar_tag));
        gmi_add_analytic_region (model, ar_tag);
        gmi_ent* ar = gmi_find(model,3, ar_tag);
        // if (i==0)
        {
          gmi_ent* af_org=af_org_vec.at(iloop-1);
          newModelEnts[af_org]=std::make_pair(ar, af_next);
        }
      }
    }
  }
}

void gmi_print(gmi_model* model)
{
  gmi_ent* g;
  gmi_iter* gi;
  for (int dim=0; dim<3; ++dim)
  {
    gi = gmi_begin(model, dim);
    while( (g = gmi_next(model, gi)) )
    {
      switch (dim)
      { 
        case 0: std::cout<< PCU_Comm_Self()<<"] geom vertex "<< gmi_tag(model, g)<<"\n";
                break;
        case 1: std::cout<< PCU_Comm_Self()<<"] geom edge "<< gmi_tag(model, g)<<"\n";
                break;
        case 2: std::cout<< PCU_Comm_Self()<<"] geom face "<< gmi_tag(model, g)<<"\n";
                break;
        default: break;
      }
    } 
  }
  gmi_end(model, gi); // end the iterator
}

// *********************************************************
void m3dc1_model::set_num_plane(int factor)
// *********************************************************
{
  num_plane = factor;
  int self = PCU_Comm_Self();
  group_size = PCU_Comm_Peers()/factor;
  local_planeid = self/group_size; // divide

  prev_plane_partid = (PCU_Comm_Self()-group_size)%PCU_Comm_Peers();
  if (prev_plane_partid<0)
    prev_plane_partid = prev_plane_partid+PCU_Comm_Peers();
  next_plane_partid = (self+group_size)%PCU_Comm_Peers();

  if (factor==1) return;

  int groupRank = self%group_size; // modulo

  MPI_Comm groupComm;
  MPI_Comm_split(oldComm, local_planeid, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}


// *********************************************************
void m3dc1_model::set_phi(double min_val, double max_val)
// *********************************************************
{
  if (!phi) phi = new double[num_plane];
  for(int i=0;i<num_plane;++i)
    phi[i] = i*(max_val-min_val)/(num_plane-1)+min_val;
}

// *********************************************************
double m3dc1_model::get_phi(int part_id)  // internal use
// *********************************************************
{
  if (!phi) return 0.0;
  return phi[part_id/group_size];
}

// *********************************************************
void m3dc1_model::set_phi(int plane_id, double p)   // api use
// *********************************************************
{
  if (!phi)
  {
    phi = new double[num_plane];
    for(int i=0;i<num_plane;++i)
      phi[i] = 0.0;
  }
  if (plane_id<0 || plane_id>num_plane) // plane_id not specified
    phi[PCU_Comm_Self()/group_size] = p;
  else
    phi[plane_id] = p;
}

// *********************************************************
void m3dc1_model::get_phi(int plane_id, double* p)   // api use
// *********************************************************
{
  if (!phi)
  {
    *p = 0.0;
  }
  else
  {
    if (plane_id<0 || plane_id>num_plane) // plane_id not specified
      *p = phi[PCU_Comm_Self()/group_size];
    else
      *p = phi[plane_id];
  }
}

// **********************************************
gmi_ent* m3dc1_model :: geomEntNextPlane(gmi_ent* gent)
// **********************************************
{
  return newModelEnts[gent].second;
} 

// **********************************************
gmi_ent* m3dc1_model:: geomEntBtwPlane(gmi_ent* gent)
// **********************************************
{
  return newModelEnts[gent].first;
}

// **********************************************
void m3dc1_model:: caculateBoundingBox()
// **********************************************
{
  const double largevalue=1.0e30;
  double max[3], min[3];
  max[0]=max[1]=max[2]=-1.0*largevalue;
  min[0]=min[1]=min[0]=largevalue;
  gmi_iter* giter = gmi_begin(model, 1);
  while(gmi_ent* gedge = gmi_next(model, giter))
  {
    double paraRange[2];
    gmi_range(model, gedge, 0, paraRange);
    double delta_param = (paraRange[1]-paraRange[0])/101.;
    int i;
    double xyz[3]={0,0,0};
    for(i=0; i<102; i++)
    {
      double parm = paraRange[0] + delta_param*(double)i;
      gmi_eval(model, gedge, &parm,xyz);
      if(xyz[0]>max[0]) max[0]=xyz[0];
      if(xyz[1]>max[1]) max[1]=xyz[1];
      if(xyz[2]>max[2]) max[2]=xyz[2];
      if(xyz[0]<min[0]) min[0]=xyz[0];
      if(xyz[1]<min[1]) min[1]=xyz[1];
      if(xyz[2]<min[2]) min[2]=xyz[2];
    }
  }
  gmi_end(model, giter);
  m3dc1_model::instance()->boundingBox[0]=min[0];
  m3dc1_model::instance()->boundingBox[1]=min[1];
  m3dc1_model::instance()->boundingBox[2]=max[0];
  m3dc1_model::instance()->boundingBox[3]=max[1];
}
/* DGESV prototype */
extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                 double* b, int* ldb, int* info );

// **********************************************
void interpolateCubicBSpline( vector<double>& points,vector<double>& knots, vector<double> &ctrlPoints, int bc)
// **********************************************
{
  int numPts=points.size();
  int order_p=4;
  assert(numPts>1);
  ctrlPoints.resize(numPts+2);
  vector< double> coeffs((numPts+2)*(numPts+2),0.); // numPts+2 * numPts+2 linear system
  // first find the constraint to coninside with points
  // 2 natural cubic spline at the boundary
  for( int i=0; i<numPts+2; i++)
  {
    vector<double> points_tmp(numPts+2,0.);
    points_tmp.at(i)=1.0;
    vector<double>  weight_p;
    M3DC1::BSpline  basis(order_p, points_tmp, knots, weight_p);
    for( int j=0; j<numPts; j++)
    {
      double para=knots.at(order_p+j-1);
      double res= basis.eval(para);
      double secondDeriv0=basis.evalSecondDeriv(0);
      double secondDeriv1=basis.evalSecondDeriv(1);
      coeffs.at(i*(numPts+2)+j)=res;
    } 
    double secondDeriv0=basis.evalSecondDeriv(0);
    double secondDeriv1=basis.evalSecondDeriv(1);
    if(bc==0) // natural
    {
      coeffs.at(i*(numPts+2)+numPts)=secondDeriv0;
      coeffs.at(i*(numPts+2)+numPts+1)=secondDeriv1;
    }
    else // periodic
    {
      double firstDeriv0=basis.evalFirstDeriv(0);
      double firstDeriv1=basis.evalFirstDeriv(1);
      coeffs.at(i*(numPts+2)+numPts)=firstDeriv0-firstDeriv1;
      coeffs.at(i*(numPts+2)+numPts+1)=secondDeriv0-secondDeriv1; 
    }
  }

  // set up the linear system and solve
  vector<double> rhs(numPts+2,0.0);
  for( int i=0; i<numPts; i++)
    rhs.at(i)=points.at(i);

  int info,one=1, dim=numPts+2;
  vector<int> ipiv(dim,0);
  dgesv_( &dim, &one,& (coeffs.at(0)), &dim, &(ipiv.at(0)), &(rhs.at(0)), &dim, &info );
  assert( info==0);
  for ( int i=0; i<numPts+2; i++)
    ctrlPoints.at(i)=rhs.at(i);  
}

// **********************************************
void m3dc1_model::setupCommGroupsPlane()
// **********************************************
{
  /**get planeId where the current processor is */
  int planeId=local_planeid;
  int rank = PCU_Comm_Self();
  /**get the localrank of the current processor is */
  int localrank=rank-planeId*group_size;
  /** split MPI_COMM_WORLD, put the processors of the same planeId into one CommWorld*/
  MPI_Comm_split(MPI_COMM_WORLD,localrank,planeId, &(PlaneGroups[localrank]));
  //MPI_Barrier(MPI_COMM_WORLD);
}

// **********************************************
MPI_Comm & m3dc1_model:: getMPICommPlane()
// **********************************************
{
  int rank = PCU_Comm_Self();
  int planeId=local_planeid;
  /**get the localrank of the current processor is */
  int localrank=rank-planeId*group_size;
  return PlaneGroups[localrank];
}

// use clamped b-spline as underlying representation
// **********************************************
void attach_periodic_cubic_curve ( int* edge, int * numPts, double * points)
// **********************************************
{
  edgeType[*edge]=BSPLINE;
  std::pair< int, int> vtx=edgeContainer[*edge];
  assert(vtx.first==vtx.second);
 
  int order_p=4; 
  int knotsize=2*order_p+*numPts-2;
  vector<double> knots(knotsize,0.);
  vector<double> ctrlPointsX(*numPts+2),ctrlPointsY(*numPts+2),weight;
  for( int i=0; i<order_p; i++)
  {
    knots.at(knotsize-i-1)=1.0;
  }
  double increment=1.0/(*numPts-1);
  for (int i=0; i<*numPts-2; i++)
  {
    //double increment=inter_len.at(i)/len;
    knots.at(order_p+i)=knots.at(order_p+i-1)+increment;
  }
  vector<double> pointsX(*numPts),pointsY(*numPts);
  for( int i=0; i<*numPts; i++)
  {
    pointsX.at(i)=points[2*i];
    pointsY.at(i)=points[2*i+1];
  }
  interpolateCubicBSpline(pointsX,knots,ctrlPointsX, 1);
  interpolateCubicBSpline(pointsY,knots,ctrlPointsY, 1);
  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  data[0] = new M3DC1::BSpline(order_p,ctrlPointsX,knots, weight);
  data[1] = new M3DC1::BSpline(order_p,ctrlPointsY,knots, weight);
  data2Clean.push_back(data);
  int edgePeriodic = 1;
  double edgeRange[2] = {0.0, 1.0};
  gmi_add_analytic(m3dc1_model::instance()->model, 1, *edge, edgeFunction, &edgePeriodic, &edgeRange, data);
}

// **********************************************
void offset_point (double * src, double * normal, double * thickness, double * des)
// **********************************************
{
  for(int i=0; i<2; i++)
  {
    des[i]=src[i]+*thickness*normal[i];
  }
}

// **********************************************
void eval_position( int* edge, double* para, double* position )
// **********************************************
{
  gmi_model* model = m3dc1_model::instance()->model;
  gmi_ent* ae = gmi_find(model, 1, *edge);
  M3DC1::evalCoord( *para, position, gmi_analytic_data(model,ae));
}

// **********************************************
void eval_normal( int* edge, double* para, double* normal )
// **********************************************
{
  gmi_model* model = m3dc1_model::instance()->model;
  gmi_ent* ae = gmi_find(model, 1, *edge);
  M3DC1::Expression** expr = (M3DC1::Expression**) gmi_analytic_data(model, ae);
  M3DC1::evalNormalVector(expr[0],expr[1],*para, normal);
}

// **********************************************
void eval_curvature( int* edge, double* para, double* curvature )
// **********************************************
{
  gmi_model* model = m3dc1_model::instance()->model;
  gmi_ent* ae = gmi_find(model, 1, *edge);
  M3DC1::Expression** expr = (M3DC1::Expression**) gmi_analytic_data(model, ae);
  M3DC1::evalCurvature(expr[0],expr[1],*para, curvature);
}

// FUNCTIONS USED BY TEST/CONVERT_POLAR
// **********************************************
gmi_ent* create_model_vertex( int id, double* xyz)
// **********************************************
{
  vector<double> vtx_coord(3);
  vtx_coord[0]=xyz[0];
  vtx_coord[1]=xyz[1];
  vtx_coord[2]=0.;
  vtxContainer[id]=vtx_coord;
  return gmi_add_analytic(m3dc1_model::instance()->model, 0, id, vertexFunction, NULL, NULL, NULL);
}

// **********************************************
gmi_ent* create_b_spline_curve( int id, int order, int numPts, double* ctrlPts, double * knots, double* weight)
// **********************************************
{
  edgeType[id]=BSPLINE;
  std::pair< int, int> vtx=edgeContainer[id];
  vector<double> X_p, Y_p, knots_p;
  for (int i=0; i<numPts;i++)
  {
    X_p.push_back(ctrlPts[2*i]);
    Y_p.push_back(ctrlPts[2*i+1]);
  }
  for( int i=0; i<numPts+order; i++)
  {
    knots_p.push_back(knots[i]);
  }

  M3DC1::BSpline** data=new M3DC1::BSpline*[2];
  vector<double> wt;
  data[0] = new M3DC1::BSpline(order,X_p,knots_p, wt);
  data[1] = new M3DC1::BSpline(order,Y_p,knots_p, wt);
  data2Clean.push_back(data);
  int edgePeriodic = 0;
  if(vtx.first==vtx.second) edgePeriodic=1;
  double edgeRange[2] = {0.0, 1.0};
  gmi_ent* gedge = gmi_add_analytic(m3dc1_model::instance()->model, 1, id, edgeFunction, &edgePeriodic, &edgeRange, data);
  make_edge_topo(m3dc1_model::instance()->model,gedge, vtx.first, vtx.second);
  //if (!PCU_Comm_Self()) std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": new edge "<<id<<" vtx("<<vtx.first<<", "<<vtx.second<<")\n";
  return gedge;
}

void m3dc1_model::save_gtag()
{
  ge_tag = new int*[3];
  gmi_iter* g_it;
  gmi_ent* ge;
  int cnt;

  for (int i = 0; i < 3; ++i)
    ge_tag[i] = new int[model->n[i]];

  for (int dim=0; dim<=2; ++dim)
  {
    g_it = gmi_begin(model, dim);
    cnt=0;
    while ((ge = gmi_next(model, g_it)))
    {
      ge_tag[dim][cnt++]= gmi_tag(model,ge);
      gmi_base_set_tag (model, ge, cnt);
    }
    gmi_end(model, g_it);
  }
}

void m3dc1_model::restore_gtag()
{
    gmi_iter* g_it;
    gmi_ent* ge;

    for (int dim=0; dim<=2; ++dim)
    {
      g_it = gmi_begin(model, dim);
      int cnt=0;
      while ((ge = gmi_next(model, g_it)))
        gmi_base_set_tag (model, ge, ge_tag[dim][cnt++]);
      gmi_end(model, g_it);
    }
} 
