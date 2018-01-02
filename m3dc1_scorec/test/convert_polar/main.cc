/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
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
#ifdef M3DC1_MESHGEN
#include "SimUtil.h"
#include "SimModel.h"
#include "SimError.h"
#include "SimAdvModel.h"
#include "MeshSim.h"
#include "BSpline.h"
void save_sim_model ();
char simLic[50] = "/net/common/meshSim/license/license.txt";
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
#endif

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

  FILE *fp = fopen("POLAR", "r");
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
  apf::ModelEntity* g_face = (apf::ModelEntity*)(gmi_add_analytic(m3dc1_model::instance()->model, 2, 1/*id*/, faceFunction, facePeriodic, faceRanges, NULL));

  apf::Vector3 param(0,0,0);
  apf::Vector3 xyz_coords(0.0, 0.0, 0.0);

  fscanf(fp, "%d %d\n", &n_radial, &n_theta);
  cout<<"> convert_polar: reading \"POLAR\" - N_RADIAL = "<<n_radial<<", N_THETA = "<<n_theta<<endl;

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

  // white model and mesh file
  cout<<"> convert_polar: writing model in PUMI-readable \"model.txt\"\n";
  save_model("model.txt");
#ifdef M3DC1_MESHGEN
  cout<<"> convert_polar: writing model in Simmetrix-readable \"model.smd\"\n";
  save_sim_model();
#endif
  cout<<"> convert_polar: writing mesh in PUMI-readable \"mesh.smb\"\n";
  m->writeNative("mesh.smb");
  cout<<"> convert_polar: writing mesh in Paraview data files\n";
  apf::writeVtkFiles("mesh",m);

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

  // destroy tag and mesh m
  apf::removeTagFromDimension(m, norm_curv_tag, 0);
  m->destroyTag(norm_curv_tag);
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

#ifdef M3DC1_MESHGEN
void save_sim_model ()
{
// contact seols@rpi.edu to get the code
}
#endif
