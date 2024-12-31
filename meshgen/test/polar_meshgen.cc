/****************************************************************************** 

  (c) 2005-2023 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <iostream>
#include "apf.h"
#include "PCU.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include <gmi_analytic.h>
#include <apfMDS.h>
#include <assert.h>
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "SimUtil.h"
#include "SimModel.h"
#include "SimAdvModel.h"
#include "MeshSim.h"
#include "BSpline.h"
#include <PCU.h>
#include <SimPartitionedMesh.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_sim.h>
#include "gmi_mesh.h"
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <cstdlib>
#include <iostream>
#include <string.h>

void save_sim_model ();

#ifdef LICENSE
#include <SimLicense.h>
#ifdef STELLAR
char simLic[128]="/home/PPPL/simmetrix/license/simmetrix.lic";
#endif
#ifdef PPPL
char simLic[128]="/usr/pppl/Simmetrix/simmodsuite.lic";
#endif
#ifdef SDUMONT
char simLic[128]="/scratch/ntm/software/Simmetrix/license/simmodsuite.lic";
#endif
#else
#ifdef MIT
  char simLic[128]="/orcd/nese/psfc/001/software/simmetrix/RLMServer-14/server.lic";
#else // scorec
  char simLic[128]="/net/common/meshSim/license/license.txt";
#endif
#endif

char meshSizeFnName[]="meshSizeFn";
double meshSizeFn(const double gpt[3], void* data);

void messageHandler(int type, const char *msg);
void M_writeFMDB(pMesh theMesh, char *name, int snap);
using namespace std;
double thickness=0.02, width=2.5, height=4.;
double offsetX=0, offsetY=0;
char modelName[128]="model";
char meshName[128]="mesh";
char in_polarFile[128];
char outFile[128];
int modelType=0;
int reorder=0;
int useVacuumParams=0;
int useOriginal = 1;
int numInterPts =20;
double meshSizes[] = {0.05, 0.05, 0.05};
int useAnaltyicSize=0;
double meshGradationRate=0.3;

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    printf("Info: %s\n",msg);
    break;
  case Sim_DebugMsg:
    printf("Debug: %s\n",msg);
    break;
  case Sim_WarningMsg:
    printf("Warning: %s\n",msg);
    break;
  case Sim_ErrorMsg:
    printf("Error: %s\n",msg);
    break;
  }
  return;
}

using namespace std;
void create_edge(apf::Mesh2* m, vector<int>& g_edge_ids, apf::ModelEntity* g_face, apf::MeshEntity** ev)
{
// apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find  
  apf::MeshEntity* found = findUpward(m,apf::Mesh::EDGE,ev);
  if (!found) 
  {
    int g_dim1 = gmi_dim(m3dc1_model::instance()->model,(gmi_ent*)(m->toModel(ev[0])));
    int g_dim2 = gmi_dim(m3dc1_model::instance()->model,(gmi_ent*)(m->toModel(ev[1])));
    int vtx1, vtx2;
    if (g_dim1==0 && g_dim2==0)  
    {
      // create g geom edge 
      int g_tag1=gmi_tag(m3dc1_model::instance()->model,(gmi_ent*)(m->toModel(ev[0])));
      int g_tag2=gmi_tag(m3dc1_model::instance()->model,(gmi_ent*)(m->toModel(ev[1])));      
      int id = g_edge_ids.size();
      create_edge(&id, &g_tag1, &g_tag2);
      int order_p=2;
      int numCtrPts=2;
      double knots_p[]={0,0,1,1};
      apf::Vector3 xyz_1(0.0, 0.0, 0.0);
      apf::Vector3 xyz_2(0.0, 0.0, 0.0);
      m->getPoint(ev[0], 0, xyz_1);
      m->getPoint(ev[1], 0, xyz_2);
      double ctrlPts_p[]={xyz_1[0],xyz_1[1],xyz_2[0], xyz_2[1]};
      gmi_ent* g_edge = create_b_spline_curve(id,2,2,&ctrlPts_p[0],&knots_p[0],NULL);
      g_edge_ids.push_back(id);
      // create a mesh edge classified on the model edge
      m->createEntity(apf::Mesh::EDGE, (apf::ModelEntity*)g_edge, ev);
    }
    else // create a mesh edge classified on the model face
      m->createEntity(apf::Mesh::EDGE, g_face, ev);
  }
}

int main( int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();

  if (argc!=2)
  {
    cout<<"[m3dc1_meshgen ERROR] usage: ./polar_meshgen input_filename\n";
    MPI_Finalize();
    return 1;
  }

  FILE* infp= fopen(argv[1],"r");
  if (!infp)
  {
    cout<<"File \""<<argv[1]<<"\" not found"<<endl;
    MPI_Finalize();
    return 1;
  }
  char namebuff[512];
  while(EOF!=fscanf(infp,"%s",namebuff))
  {
    if (strcmp(namebuff,"reorder")==0) fscanf(infp,"%d",&reorder);
    if (strcmp(namebuff,"inFile")==0) fscanf(infp, "%s", in_polarFile);
    if (strcmp(namebuff,"outFile")==0)
    {
      fscanf(infp, "%s", modelName);
      strcpy(meshName, modelName);
    }

    if (strcmp(namebuff,"meshSize")==0)
    {
      fscanf(infp, "%lf",meshSizes);
      meshSizes[1] = meshSizes[2] = meshSizes[0];
    }
  }
  fclose(infp);


  FILE *fp = fopen(in_polarFile, "r");
  if (!fp)
  {
    std::cout<<"> convert_polar ERROR: failed to open the file \"POLAR\"\n";
    PCU_Comm_Free();
    MPI_Finalize();    
    return 1;
  }

  int n_theta, n_radial; 
  int i, j;

  gmi_ent* gent;
  apf::MeshEntity* ment;
  
  apf::Mesh2* m = apf::makeEmptyMdsMesh(m3dc1_model::instance()->model, 2, false);
 
  // create model face
  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};
  apf::ModelEntity* g_face = (apf::ModelEntity*)
               (gmi_add_analytic(m3dc1_model::instance()->model, 2, 1 /*id*/, 
               faceFunction, facePeriodic, faceRanges, NULL));

  apf::Vector3 param(0,0,0);
  apf::Vector3 xyz_coords(0.0, 0.0, 0.0);

  fscanf(fp, "%d %d\n", &n_radial, &n_theta);
  cout<<"> convert_polar: reading \""
      <<in_polarFile<<"\" - #radial = "<<n_radial<<", #theta = "<<n_theta<<endl;

  apf::MeshEntity* v[2][n_theta];
  apf::MeshEntity* fv[3];  
  apf::MeshEntity* ev[2];

  // read the first two layers
  for(i=0; i<n_theta; i++)
    fscanf(fp, "%lf %lf", &xyz_coords[0], &xyz_coords[1]);

  fv[0] = m->createVertex(g_face, xyz_coords, param);

  for(i=0; i<n_theta; i++)
  {
    fscanf(fp, "%lf %lf", &xyz_coords[0], &xyz_coords[1]);
    v[1][i] = m->createVertex(g_face, xyz_coords, param); 
  }

  for(i=0; i<n_theta; i++) 
  {
    fv[1] = v[1][i];
    fv[2] = v[1][(i+1)%(n_theta)];
    // create edge 0
    ev[0] = fv[0]; ev[1] = fv[1];
    apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find

    // create edge 1
    ev[0] = fv[1]; ev[1] = fv[2];
    apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
    // create edge 2
    ev[0] = fv[0]; ev[1] = fv[2];
    apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
    // create face
    apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
  }
  
  // for the rest layers
  for(j=2; j<n_radial-1; j++) 
  {    
    // read the nodes
    for(i=0; i<n_theta; i++) 
    {
      fscanf(fp, "%lf %lf", &xyz_coords[0], &xyz_coords[1]);
      v[j%2][i] = m->createVertex(g_face, xyz_coords, param); 	
    }
    
    if((j%2)) 
    {
      // create the elements
      for(i=0; i<n_theta; i++) 
      {
	fv[0] = v[1][(i)%(n_theta)];
        fv[1] = v[1][(i+1)%(n_theta)]; 
        fv[2] = v[0][(i+1)%(n_theta)];
        // create edge 0
        ev[0] = fv[0]; ev[1] = fv[1];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find
        // create edge 1
        ev[0] = fv[1]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create edge 2
        ev[0] = fv[0]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create face
	apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
	
	fv[0] = v[0][(i+1)%(n_theta)]; 
        fv[1] = v[0][(i)%(n_theta)]; 
        fv[2] = v[1][(i)%(n_theta)];
        // create edge 0
        ev[0] = fv[0]; ev[1] = fv[1];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find  
        // create edge 1
        ev[0] = fv[1]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create edge 2
        ev[0] = fv[0]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create face
	apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
      }
    }
    else
    {
      // create the elements
      for(i=0; i<n_theta; i++) 
      {
	fv[0] = v[0][(i)%(n_theta)]; 
        fv[1] = v[0][(i+1)%(n_theta)]; 
        fv[2] = v[1][(i+1)%(n_theta)];
        // create edge 0
        ev[0] = fv[0]; ev[1] = fv[1];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find  
        // create edge 1
        ev[0] = fv[1]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create edge 2
        ev[0] = fv[0]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create face

        apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);

	fv[0] = v[1][(i+1)%(n_theta)]; 
        fv[1] = v[1][(i)%(n_theta)]; 
        fv[2] = v[0][(i)%(n_theta)];
        // create edge 0
        ev[0] = fv[0]; ev[1] = fv[1];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev); // make or find  
        // create edge 1
        ev[0] = fv[1]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create edge 2
        ev[0] = fv[0]; ev[1] = fv[2];
        apf::buildElement(m, g_face, apf::Mesh::EDGE, ev);
        // create face
	apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
      }
    }
  }

  apf::MeshTag* norm_curv_tag = m->createDoubleTag("norm_curv", 3);

  // create the outer layer 
  double norm_curv[3];
//  vector<double> interpolate_points;
//  double mid[]={0.,0.}; 
//  interpolate_points.resize(2*n_theta);
//  double bdbox[4]={1e99,1e99, -1e99,-1e99};
  apf::ModelEntity* g_vertex;
  for(i=0; i<n_theta; i++) 
  {
    fscanf(fp, "%lf %lf", &xyz_coords[0], &xyz_coords[1]);
//    interpolate_points.at(2*i)=xyz_coords[0];
//    interpolate_points.at(2*i+1)=xyz_coords[1];
//    if(xyz_coords[0]<bdbox[0]) bdbox[0]=xyz_coords[0];
//    if(xyz_coords[0]>bdbox[2]) bdbox[2]=xyz_coords[0];
//    if(xyz_coords[1]<bdbox[1]) bdbox[1]=xyz_coords[1];
//    if(xyz_coords[1]>bdbox[3]) bdbox[3]=xyz_coords[1];
    // create a model vertex
    g_vertex = (apf::ModelEntity*)(create_model_vertex(i,&xyz_coords[0]));
    // create mesh vertex
    v[j%2][i] = m->createVertex(g_vertex, xyz_coords, param); 
    fscanf(fp, "%lf %lf %lf", norm_curv, norm_curv+1, norm_curv+2);
    m->setDoubleTag(v[j%2][i], norm_curv_tag, &norm_curv[0]);
  }
//  mid[0]=(bdbox[0]+bdbox[2])/2.;
//  mid[1]=(bdbox[1]+bdbox[3])/2.;
//  double left[]={mid[0]-width/2.,mid[1]}, right[]={mid[0]+width/2.,mid[1]}; 
  vector<int> g_edge_ids;

  if ((j%2)) 
  {
    // create the elements
    for(i=0; i<n_theta; i++) 
    {
      fv[0] = v[1][(i)%(n_theta)]; 
      fv[1] = v[1][(i+1)%(n_theta)]; 
      fv[2] = v[0][(i+1)%(n_theta)];
      // create edge 0
      ev[0] = fv[0]; ev[1] = fv[1];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 1
      ev[0] = fv[1]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 2
      ev[0] = fv[0]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create face
      apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
      
      fv[0] = v[0][(i+1)%(n_theta)]; 
      fv[1] = v[0][(i)%(n_theta)]; 
      fv[2] = v[1][(i)%(n_theta)];
      // create edge 0
      ev[0] = fv[0]; ev[1] = fv[1];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 1
      ev[0] = fv[1]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 2
      ev[0] = fv[0]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create face
      apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
    }
  }
  else
  {
    // create the elements
    for(i=0; i<n_theta; i++)
    {
      fv[0] = v[0][(i)%(n_theta)]; 
      fv[1] = v[0][(i+1)%(n_theta)]; 
      fv[2] = v[1][(i+1)%(n_theta)];
      // create edge 0
      ev[0] = fv[0]; ev[1] = fv[1];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 1
      ev[0] = fv[1]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 2
      ev[0] = fv[0]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create face
      apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
      
      fv[0] = v[1][(i+1)%(n_theta)]; 
      fv[1] = v[1][(i)%(n_theta)]; 
      fv[2] = v[0][(i)%(n_theta)];
      // create edge 0
      ev[0] = fv[0]; ev[1] = fv[1];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 1
      ev[0] = fv[1]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create edge 2
      ev[0] = fv[0]; ev[1] = fv[2];
      create_edge (m, g_edge_ids, g_face, ev);
      // create face
      apf::buildElement(m, g_face, apf::Mesh::TRIANGLE, fv);
    }
  }
  fclose(fp);

  int one=1;
  assert(n_theta == g_edge_ids.size());
  create_loop(&one,&n_theta,&g_edge_ids[0]);
  set_inner_wall_boundary (&one);

  m->acceptChanges();
  m->verify();

  save_sim_model();

/*
  cout<<"> convert_polar: writing normal/curvature in \"norm_curv\"\n";
  fp = fopen("norm_curv", "w"); 
  gmi_iter* ge_it = gmi_begin(m3dc1_model::instance()->model, 1);
  while ((gent = gmi_next(m3dc1_model::instance()->model, ge_it))) 
  {
    // get mesh vertices classified on the model vertex
    apf::MeshIterator* mv_it;
    mv_it = m->begin(0);
    while ((ment=m->iterate(mv_it))!=0) 
    {
      if ((gmi_ent*)(m->toModel(ment)) != gent) continue;
      assert(m->hasTag(ment, norm_curv_tag));
      m->getDoubleTag(ment, norm_curv_tag, &norm_curv[0]);
      m->getPoint(ment, 0, xyz_coords);
      fprintf(fp, "%d %lf %lf %lf %lf %lf \n",  getMdsIndex(m, ment), 
              xyz_coords[0], xyz_coords[1], norm_curv[0], norm_curv[1], norm_curv[2]);
    }
  }
  gmi_end(m3dc1_model::instance()->model, ge_it);
  fclose(fp);

  cout<<"> convert_polar: model #v "<<m3dc1_model::instance()->model->n[0]
      <<" #e "<<m3dc1_model::instance()->model->n[1] 
      <<" #f "<<m3dc1_model::instance()->model->n[2]<<"\n";

  // destroy tag and mesh m
  apf::removeTagFromDimension(m, norm_curv_tag, 0);
  m->destroyTag(norm_curv_tag);
*/

  m->destroyNative();
  apf::destroyMesh(m);
  
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

int make_sim_model (pGModel& sim_model)
{
  std::map< int, std::vector<int> > loopContainer;
  std::map<int, std::pair<int, int> > edgeContainer;
  std::map<int, int> edgeType;
  std::map<int, std::vector<double> > vtxContainer;
  export_model_data(vtxContainer, edgeType, edgeContainer, loopContainer);

  std::map<int, pGVertex> vertices;
  std::map<int, pGEdge> edges;
  gmi_model* model = m3dc1_model::instance()->model;

#ifdef SIM12
  pGIPart part = GM_part(sim_model);
#else // SIMMODSUITE_MAJOR_VERSION >= 15
  pGIPart part = GM_rootPart(sim_model);
#endif

  pGRegion outerRegion = GIP_outerRegion(part);
  // Now we'll add loops
  vector<pGEdge> innerLoop;
  vector<pGEdge> currentLoop;
  for( std::map<int, vector<int> >:: iterator it=loopContainer.begin(); it!=loopContainer.end(); it++)
  {
    int numE=it->second.size();
    // First we'll add the vertices on the loop
    for( int i=0; i<numE; i++)
    {
      int edge=it->second[i];
      std::pair<int, int> vtx=edgeContainer[edge];
      std::vector<double>& xyz= vtxContainer.at(vtx.first);
#ifdef SIM12
      vertices[vtx.first] = GIP_insertVertexInRegion(part, &xyz[0], outerRegion);
#else 
      vertices[vtx.first] = GR_createVertex(GIP_outerRegion(part), &xyz[0]);
#endif
    }

    for( int j=0; j<numE; j++)
    {
      pGVertex startVert, endVert;
      int edge = it->second[j];
      gmi_ent* ae = gmi_find(model, 1, edge);

      int curvType=edgeType.at(edge);
      std::pair<int, int> vtx=edgeContainer[edge];
      pCurve curve=NULL;
      startVert=vertices[vtx.first];
      endVert=vertices[vtx.second];
      if ( startVert== NULL || endVert==NULL)
      {
        std::cout<<" invalid edge vertex not created edge tag "<<edge<<std::endl;
        throw 1;
      }
          {
            double xyz1[3], xyz2[3];
            GV_point(startVert,xyz1);
            GV_point(endVert,xyz2);
            curve = SCurve_createLine(xyz1, xyz2);
          }
#ifdef SIM12
      pGEdge pe = GIP_insertEdgeInRegion(part, startVert, endVert, curve, 1, outerRegion);
#else
      pGEdge pe = GR_createEdge(GIP_outerRegion(part), startVert, endVert, curve, 1);
#endif
      edges[edge]=pe;
      currentLoop.push_back( pe);
    }

    // Now add the faces
    double corner[3]={0,0,0}, xPt[3]={1,0,0}, yPt[3]={0,1,0};  // the points defining the surface of the face
    vector<pGEdge> faceEdges;
    vector<int> faceDirs;
    int loopDef[2] = {0,0};
    int numloops=1;
    if (GM_numFaces(sim_model)!=0)
    {
      numloops=2;
      int innerEdgeNum=innerLoop.size();
      loopDef[1]=innerEdgeNum;
      for(int k=innerEdgeNum-1; k>=0; k--)
      {
        faceEdges.push_back(innerLoop.at(k));
        faceDirs.push_back(0);
      }
    }
    int currentEdgeNum=currentLoop.size();
    for(int k=0; k<currentEdgeNum; k++)
    {
      faceEdges.push_back(currentLoop.at(k));
      faceDirs.push_back(1);
    }
    innerLoop=currentLoop;
    currentLoop.clear();
    pSurface planarSurface;

    planarSurface = SSurface_createPlane(corner,xPt,yPt);
#ifdef SIM12
    GIP_insertFaceInRegion(part,faceEdges.size(),&(faceEdges[0]),&(faceDirs[0]),
      numloops,loopDef,planarSurface,1,outerRegion);
#else
    GR_createFace(GIP_outerRegion(part), faceEdges.size(),&(faceEdges[0]),&(faceDirs[0]),numloops,loopDef,planarSurface,1);
#endif
  }
  printf("Number of vertices in Simmetrix model: %d\n",GM_numVertices(sim_model));
  printf("Number of edges in Simmetrix model: %d\n",GM_numEdges(sim_model));
  printf("Number of faces in Simmetrix model: %d\n",GM_numFaces(sim_model));
  printf("Number of regions in Simmetrix model: %d\n",GM_numRegions(sim_model));
}


void save_sim_model ()
{
  gmi_ent* gent;
  apf::MeshEntity* ment;

  char filename[128];
  char model_filename[128];
  char mesh_filename[128];

#ifdef LICENSE
  SimLicense_start("geomsim_core,geomsim_adv,meshsim_surface,meshsim_adv",simLic);
#else
  // for SCOREC
  Sim_readLicenseFile(simLic);
#endif

  gmi_register_mesh();
  Sim_logOn("m3dc1_meshgen.log");
  SimModel_start();

  // Tessellation of GeomSim geometry requires Meshing to have started
  MS_init();

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  //pGModel sim_model = gmi_export_sim(model);
  //sprintf(model_filename,"%s.smd",filename);
  //GM_write(sim_model,model_filename,0,0);
#ifdef SIM12
  pGModel sim_model = GM_new();
#else
  pGModel sim_model = GM_new(1);
#endif

  make_sim_model (sim_model);

  //pParMesh sim_mesh = PM_new(0,sim_model,PMU_size());
  pMesh sim_mesh = M_new(0,sim_model);
  pACase meshCase = MS_newMeshCase(sim_model);

  int iface=0;
  GFIter faces = GM_faceIter(sim_model);
  int nGFace = GM_numFaces(sim_model);
  while( pGFace face= GFIter_next(faces))
  {
    MS_setMeshSize(meshCase,face,1, meshSizes[iface++],NULL);
    cout<<"set mesh size of face "<<iface-1<<" to "<<meshSizes[iface-1]<<"\n";
  }
  assert(iface<=3);
  GFIter_delete(faces);

  pSurfaceMesher surfMesh = SurfaceMesher_new(meshCase,sim_mesh);
  MS_setGlobalSizeGradationRate(meshCase, meshGradationRate);

  SurfaceMesher_execute(surfMesh,progress);
  SurfaceMesher_delete(surfMesh);

  std::cout<<"\n[INFO] # model entities: V "<<GM_numVertices(sim_model)
                          <<", E "<<GM_numEdges(sim_model)
                          <<", F "<<GM_numFaces(sim_model)<<"\n";
  std::cout<<"[INFO] # mesh  entities: V "<<M_numVertices(sim_mesh)
                          <<", E "<<M_numEdges(sim_mesh)
                          <<", F "<<M_numFaces(sim_mesh) <<"\n\n";

  // convert simmetrix to PUMI
  gmi_sim_start();
  gmi_register_sim();
 
  sprintf(model_filename,"%s.dmg", modelName);
  cout <<"polar_meshgen: model is saved in PUMI-readable \""<<model_filename<<"\"\n";
  gmi_write_dmg(m3dc1_model::instance()->model, model_filename);

  sprintf(model_filename,"%s.txt", modelName);
  cout <<"polar_meshgen: model is saved in M3DC1-readable \""<<model_filename<<"\"\n";
  save_model(model_filename);

  sprintf(model_filename,"%s.smd", modelName);
  cout <<"polar_meshgen: model is saved in Simmetrix-readable \""<<model_filename<<"\"\n";
  GM_write(sim_model,model_filename,0,0);

  int num_mfaces=M_numFaces(sim_mesh);

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.sms", meshName, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.sms", meshName, num_mfaces);
  M_write(sim_mesh,mesh_filename,0,0);

  cout<<"\npolar_meshgen: mesh is saved in Simmetrix-readable \""<<mesh_filename<<"\"\n";

  MS_deleteMeshCase(meshCase);
  M_release(sim_mesh);
  Progress_delete(progress);

  progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pParMesh sim_pmesh = PM_load(mesh_filename, sim_model, progress);

  apf::Mesh* simApfMesh = apf::createMesh(sim_pmesh);
  apf::Mesh2* mesh = apf::createMdsMesh(m3dc1_model::instance()->model, simApfMesh); //, reorder);

  // FIXME: reorder doesn't work yet
  /*
  if (reorder)
  {
    std::cout<<"\n[INFO] adjacency-based mesh reordering in PUMI is ON\n\n";
    std::cout <<"m3dc1_meshgen: initial mesh before reordering is saved in \"initial-mesh.vtk\" for Paraview\n";

    apf::Numbering* nn=apf::numberOwnedDimension(simApfMesh, "elem", 2);
    writeVtkFiles("initial-mesh.vtk", simApfMesh);
    destroyNumbering(nn);
  }
  else
    std::cout<<"\n[INFO] adjacency-based mesh reordering in PUMI is OFF\n\n";
  */

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.smb", meshName, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.smb", meshName, num_mfaces);

  std::cout <<"polar_meshgen: mesh is saved in PUMI-readable \""<<mesh_filename<<"\"\n";
  mesh->writeNative(mesh_filename);

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.vtk", meshName, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.vtk", meshName, num_mfaces);
  std::cout <<"polar_meshgen: mesh is saved in \""<<mesh_filename<<"\" for Paraview\n";

  apf::Numbering* nr=apf::numberOwnedDimension(mesh, "elem", 2);
  writeVtkFiles(mesh_filename, mesh);
  destroyNumbering(nr);
  std::cout<<"\n<< Check the attribute \"elem\" in Paraview for the element order! >>\n";

  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);

  cout<<"\n<< Continue with \"simmodeler\" for more meshing control! >>\n\n";

/*

  cout<<"> convert_polar: writing normal/curvature in \"norm_curv\"\n";
  FILE *fp = fopen("norm_curv", "w");
  gmi_iter* ge_it = gmi_begin(m3dc1_model::instance()->model, 1);
  while ((gent = gmi_next(m3dc1_model::instance()->model, ge_it)))
  {
    // get mesh vertices classified on the model vertex
    apf::MeshIterator* mv_it;
    mv_it = mesh->begin(0);
    while ((ment=mesh->iterate(mv_it))!=0)
    {
      if ((gmi_ent*)(mesh->toModel(ment)) != gent) continue;
      assert(mesh->hasTag(ment, norm_curv_tag));
      mesh->getDoubleTag(ment, norm_curv_tag, &norm_curv[0]);
      mesh->getPoint(ment, 0, xyz_coords);
      fprintf(fp, "%d %lf %lf %lf %lf %lf \n",  getMdsIndex(mesh, ment),
              xyz_coords[0], xyz_coords[1], norm_curv[0], norm_curv[1], norm_curv[2]);
    }
  }
  gmi_end(m3dc1_model::instance()->model, ge_it);
  fclose(fp);
  
  // destroy tag and mesh m
  apf::removeTagFromDimension(m, norm_curv_tag, 0);
  m->destroyTag(norm_curv_tag);
*/

  gmi_sim_stop();

 // clean-up
  GM_release(sim_model);
  M_release(sim_pmesh);
  Progress_delete(progress);

  MS_exit();
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_logOff();
  Sim_unregisterAllKeys();

#ifdef LICENSE
  SimLicense_stop();
#endif
}
