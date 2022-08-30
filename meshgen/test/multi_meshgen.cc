/****************************************************************************** 

  (c) 2005-2022 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "BSpline.h"
#include "PolyNomial.h"
#include "PCU.h"
#include "gmi_analytic.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <iomanip>

#include "SimUtil.h"
#include "SimModel.h"
#include "SimError.h"
#include "SimAdvModel.h"
#include "SimDisplay.h"
#include "MeshSim.h"

#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <gmi.h>
#include <gmi_sim.h>
#include "gmi_mesh.h"
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include "apfShape.h"
#include <ma.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

// for PPPL
#ifdef PPPL
#include <SimLicense.h>
char simLic[128]="/usr/pppl/Simmetrix/simmodsuite.lic";
#else
// for SCOREC
char simLic[128]="/net/common/meshSim/license/license.txt";
#endif

extern pGVertex GE_insertVertex(pGEdge, double);

int make_sim_model (pGModel& sim_model, vector< vector<int> >& rgn_bdry);
int make_sim_model_old(pGModel& sim_model, vector< vector<int> >& rgn_bdry);

void messageHandler(int type, const char *msg);
using namespace std;
double thickness=0.02, width=2.5, height=4.;
double offsetX=0, offsetY=0;
char modelName[256];
char in_bdryFile[256];
int numBdry=0, numRgn=1;
bool addVacuum = false;
char bdryFile[128];
vector<string> bdry_files;
char out_bdryFile[256];
int modelType=0;
int reorder=0;
int useVacuum = 0;
int useVacuumParams=1;
int useOld=0;
int useOriginal = 1;
int numInterPts =20;
double* meshSizes;
double meshGradationRate=0.3;
double a_param, b_param, c_param, d_param, e_param;
double x, y;
double bdbox[4]={1e99,1e99, -1e99,-1e99};
void aexp(double u, double *xyz)
{
  xyz[0] = a_param + b_param*(cos(u + c_param*sin(u)));
  xyz[1] = d_param + e_param*sin(u);
  xyz[2] = 0.;
}

int gv1_id=1, gv2_id=2, ge1_id=1, ge2_id=2, loop_id=1, num_ge=2;
double knots[]={0,0,0,0,0,0,1,1,1,1,1,1};
int order_p=6;
double ctr_pts[12];
int numCtrPts=6;
int num_pts;
double mid[]={0.,0.};
double pt_pre[2], mergeTol;
double dir1[2], dir2[2], normal1[2],normal2[2];
double paras[]={0.,1.};
double ctr_pts_o[12];

vector<double> interpolate_points;

void get_vacuum_rgn();
void get_multi_rgn();

void inner_wall_pts(int num_pts);

void inner_outer_wall_pts();

bool curveOrientation(std::vector <double> &curvePts);
int outerLoop(std::vector <int> loops, std::map< int, std::vector<pGEdge> > edgesOnLoop);
void setParallelVertices(int loopInner, int loopOuter, std::map< int, std::vector<pGEdge> > edgesOnLoop);
int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();

  if (argc<2) 
  {    
    cout<<"[m3dc1_meshgen ERROR] usage: ./m3dc1_meshgen input_filename\n"; 
    MPI_Finalize();
    return 1;
  }

  if (argc==3)
    useOld = atoi(argv[2]);

  FILE* infp= fopen(argv[1],"r");
  if (!infp)
  {
    cout<<"File \""<<argv[1]<<"\" not found"<<endl;
    MPI_Finalize();
    return 1;
  }
  char namebuff[512];
  vector< vector<int> > rgn_bdry;
  int index=0, numLoop, loopID;  

  while(EOF!=fscanf(infp,"%s",namebuff))
  {
    if (strcmp(namebuff,"numBdry")==0) 
      fscanf(infp,"%d",&numBdry);
    if (strcmp(namebuff,"numRgn")==0) 
    {
      fscanf(infp,"%d",&numRgn); 
      rgn_bdry.resize(numRgn);
      meshSizes = new double[numRgn];
      for (int i=0; i<numRgn; ++i)
        meshSizes[i] = 0.05;
    }
    if (strcmp(namebuff,"bdryFile")==0) 
    {
      fscanf(infp, "%s", bdryFile);
      bdry_files.push_back(bdryFile);
    }

   if (strcmp(namebuff,"rgnBdry")==0)
   {
     fscanf(infp,"%d",&numLoop);  
     for (int i=0; i<numLoop; ++i)
     {
       fscanf(infp,"%d",&loopID);
       rgn_bdry[index].push_back(loopID);
     }
     ++index;
   }
   
    if (strcmp(namebuff,"modelType")==0) fscanf(infp,"%d",&modelType);
    if (strcmp(namebuff,"reorder")==0) fscanf(infp,"%d",&reorder);
    if (strcmp(namebuff,"numVacuumPts")==0) fscanf(infp,"%d",&numInterPts);
    if (strcmp(namebuff,"inFile")==0)
    {
      fscanf(infp, "%s", bdryFile);
      bdry_files.push_back(bdryFile);
    }
    if (strcmp(namebuff,"outFile")==0) fscanf(infp, "%s", modelName);
    if (strcmp(namebuff,"resistive-width")==0) fscanf(infp, "%lf", &thickness); // ignored if modelType==4
    if (strcmp(namebuff,"vacuum-width")==0) fscanf(infp,"%lf",&width);
    if (strcmp(namebuff,"vacuum-height")==0) fscanf(infp,"%lf",&height);
    if (strcmp(namebuff,"plasma-offsetX")==0) fscanf(infp,"%lf",&offsetX);
    if (strcmp(namebuff,"plasma-offsetY")==0) fscanf(infp,"%lf",&offsetY);
    if (strcmp(namebuff,"meshSize")==0) 
    {
      for(int i=0; i<numRgn; i++)
        fscanf(infp, "%lf",meshSizes+i);
    }
    if (strcmp(namebuff,"useVacuumParams")==0) fscanf(infp, "%d",&useVacuumParams);
    if (strcmp(namebuff,"vacuumParams")==0) 
      fscanf(infp, "%lf %lf %lf %lf %lf\n", &a_param, &b_param, &c_param, &d_param, &e_param);
    if (strcmp(namebuff,"meshGradationRate")==0) fscanf(infp, "%lf",&meshGradationRate);
    if (strcmp(namebuff,"useVacuum")==0) fscanf(infp, "%d",&useVacuum);
  }
  fclose(infp);
  assert(numBdry==bdry_files.size());

  mergeTol=0.3*thickness;

  // *************************
  //     interpolate points
  // *************************
  if (numBdry==0) // no boundary files available
  {
    //  five doubles x0, x1, x2, z0, z1 for analytic expression should be provided
    // x = x0 + x1*cos(theta+x2*sin(theta))
    // z = z0 + z1*sin(theta)
    assert(useVacuumParams==1); // five doubles should be provided
    get_vacuum_rgn();
  }
  else
  {
    get_multi_rgn();
 
    // increase ID's for geometric vertices and edges
    gv1_id=gv1_id+1;
    gv2_id=gv1_id+1;
    ge1_id=ge1_id+1;
    ge2_id=ge1_id+1;

    if (!useVacuumParams && useVacuum) 
    {
      double left[]={mid[0]-width/2.,mid[1]}, right[]={mid[0]+width/2.,mid[1]}; 
      double offset_t=0.3*height;
      double ctrPtsBottom[]={left[0],left[1],left[0],left[1]-0.3*height-offset_t,
               left[0],left[1]-0.6*height-offset_t,right[0],right[1]-0.6*height,
               right[0],right[1]-0.3*height,right[0],right[1]};
        double ctrPtsTop[]={right[0],right[1],right[0],right[1]+0.3*height,
               right[0],right[1]+0.6*height,left[0],left[1]+0.6*height+offset_t,
               left[0],left[1]+0.3*height+offset_t,left[0],left[1]};
        create_vtx(&gv1_id,left);
        cout<<"create_vtx "<<gv1_id<<"\n";
        create_vtx(&gv2_id,right);
        cout<<"create_vtx "<<gv2_id<<"\n";
        create_edge(&ge1_id, &gv1_id, &gv2_id);
        cout<<"create_edge "<<ge1_id<<" ("<<gv1_id<<", "<<gv2_id<<")\n";
        create_edge(&ge2_id, &gv2_id, &gv1_id);
        cout<<"create_edge "<<ge2_id<<" ("<<gv2_id<<", "<<gv1_id<<")\n";
        attach_b_spline_curve(&ge1_id,&order_p,&numCtrPts,ctrPtsBottom,knots,NULL);
        attach_b_spline_curve(&ge2_id,&order_p,&numCtrPts,ctrPtsTop,knots,NULL);
        int vacuumEdges[]={ge1_id,ge2_id};
        num_ge = 2;
        create_loop(&loop_id,&num_ge,vacuumEdges);
        cout<<"create_loop "<<loop_id<<" ("<<ge1_id<<", "<<ge2_id<<")\n";
        set_vacuum_boundary (&loop_id);
         ++loop_id;
    }
    else if (useVacuumParams && useVacuum)
    {
      vector<double> interpolate_points_vacuum;
      bdbox[0]=a_param-b_param;
      bdbox[2]=a_param+b_param;
      bdbox[1]=d_param-e_param;
      bdbox[3]=d_param+e_param;

      int num_pts=numInterPts;
      interpolate_points_vacuum.resize(2*num_pts);
      for(int i=0; i<num_pts; i++)
      {
        double para =2*M3DC1_PI/(num_pts-1)*i;
        double xyz[3];
        aexp(para, xyz);
        interpolate_points_vacuum.at(2*i)=xyz[0];
        interpolate_points_vacuum.at(2*i+1)=xyz[1];
      }

      create_vtx(&gv1_id,&(interpolate_points_vacuum.at(0)));
      cout<<"create_vtx "<<gv1_id<<"\n";
      create_edge(&ge1_id, &gv1_id, &gv1_id);
      cout<<"create_edge "<<ge1_id<<" ("<<gv1_id<<", "<<gv1_id<<")\n";
      attach_periodic_cubic_curve(&ge1_id,&num_pts,&(interpolate_points_vacuum.at(0)));
      int vacuumEdges[]={ge1_id};
      num_ge=1;
      create_loop(&loop_id,&num_ge,vacuumEdges);
      cout<<"create_loop "<<loop_id<<" ("<<ge1_id<<")\n";
      set_vacuum_boundary (&loop_id);
      ++loop_id;
    }
  } // numBdry>0

/*
  int numEdge = m3dc1_model::instance()->model->n[1];
  cout<<"numEdge = "<<numEdge<<"\n";
  int numplot=5000;
  double diff=0;
  FILE* fp_t=fopen("plot_geo","w");
  for( int iedge=1; iedge<=numEdge; iedge++)
  {
    for(int j=0; j<numplot+1; j++)
    {
      double xyz[2];
      double para=double(j)/double(numplot);
      eval_position(&iedge,&para,xyz);
      double normal[2];
      eval_normal(&iedge,&para,normal);
      double curv;
      eval_curvature(&iedge,&para,&curv);
      fprintf(fp_t,"%f %f %f %f %f\n", xyz[0],xyz[1],normal[0],normal[1],curv);
      if (modelType==0)
      {
        double xyz_phy[3];
        aexp(para*2*M3DC1_PI, xyz_phy);
        diff=max(diff,getDist2D(xyz,xyz_phy));
      }
    }
  }
  fclose(fp_t);

  if (useVacuumParams) cout<<"[m3dc1_meshgen INFO] interpolation error "<<diff<<endl;
*/

  char filename[256];
  char model_filename[256];
  char mesh_filename[256];

  if (modelType==3)
    sprintf(filename,"%s-%0.2f-%0.1f-%0.1f", modelName, thickness, width, height);
  else if (modelType==4)
    sprintf(filename,"%s-%0.1f-%0.1f", modelName, width, height);
  else
    sprintf(filename,"%s", modelName);

  sprintf(model_filename,"%s.txt", filename);
  save_model(model_filename);

#ifdef PPPL
  // for PPPL
  SimLicense_start("geomsim_core,geomsim_adv,meshsim_surface,meshsim_adapt,meshsim_adv",simLic);
#else
  // for SCOREC
  Sim_readLicenseFile(simLic);
#endif

  Sim_logOn("m3dc1_meshgen.log");
  SimModel_start();

  // Tessellation of GeomSim geometry requires Meshing to have started
  MS_init();

  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel sim_model = GM_new(1);
  if (useOld)
    make_sim_model_old (sim_model, rgn_bdry);
  else
    make_sim_model(sim_model, rgn_bdry);
  sprintf(model_filename,"%s.smd",filename);
  GM_write(sim_model,model_filename,0,0);


  //pParMesh sim_mesh = PM_new(0,sim_model,PMU_size());
  pMesh sim_mesh = M_new(0,sim_model);
  pACase meshCase = MS_newMeshCase(sim_model);

  int iface=0;
  GFIter faces = GM_faceIter(sim_model);
  int nGFace = GM_numFaces(sim_model);
  while( pGFace face= GFIter_next(faces))
  {
    MS_setMeshSize(meshCase,face, 1, meshSizes[iface++], NULL);
    cout<<"set mesh size of face "<<iface-1<<" to "<<meshSizes[iface-1]<<"\n";
  }
  GFIter_delete(faces);

  pSurfaceMesher surfMesh = SurfaceMesher_new(meshCase,sim_mesh);
  MS_setGlobalSizeGradationRate(meshCase, meshGradationRate);
  // for mesh matching
  int pbc=0;
  if (pbc)
  {
    pGEntity srcEdge = GM_entityByTag(sim_model,1,2);
    pGEntity destEdge  = GM_entityByTag(sim_model,1,4);
    //MS_setMeshMatchDegrees(meshCase,srcEdge,destEdge,NULL,NULL,NULL,0.0,NULL);
  }

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
  gmi_register_mesh();
  gmi_register_sim();

  sprintf(model_filename,"%s.smd",filename);
  gmi_model* mdl = gmi_load(model_filename);
 
  sprintf(model_filename,"%s.dmg",filename);
  gmi_write_dmg(mdl, model_filename);
  
  sprintf(model_filename,"%s.txt",filename);
  cout<<"m3dc1_meshgen: model is saved in M3D-C1-readable \""<<model_filename<<"\"\n";

  sprintf(model_filename,"%s.smd",filename);
  cout<<"m3dc1_meshgen: model is saved in Simmetrix-readable \""<<model_filename<<"\"\n";

  sprintf(model_filename,"%s.dmg",filename);
  cout<<"m3dc1_meshgen: model is saved in PUMI-readable \""<<model_filename<<"\"\n";

  int num_mfaces=M_numFaces(sim_mesh);

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.sms",filename, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.sms",filename, num_mfaces);
  M_write(sim_mesh,mesh_filename,0,0);
  cout<<"\nm3dc1_meshgen: mesh is saved in Simmetrix-readable \""<<mesh_filename<<"\"\n";

  MS_deleteMeshCase(meshCase);
  M_release(sim_mesh);
  Progress_delete(progress);

  progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pParMesh sim_pmesh = PM_load(mesh_filename, sim_model, progress);

  apf::Mesh* simApfMesh = apf::createMesh(sim_pmesh);
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh); //, reorder);

  if (reorder)
  { 
    std::cout<<"\n[INFO] adjacency-based mesh reordering in PUMI is ON\n\n";
    std::cout <<"m3dc1_meshgen: initial mesh before reordering is saved in \"initial-mesh.vtk\" for Paraview\n";
    
    apf::Numbering* nn=apf::numberOwnedDimension(simApfMesh, "elem", 2);
    writeVtkFiles("initial-mesh.vtk", simApfMesh);
    destroyNumbering(nn);
  }
  else 
    std::cout<<"\n[INFO] adjacency-based mesh reordering in PUMI is OFF (\"reorder\" in \""<<argv[1]<<"\"=0)\n\n";

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.smb",filename, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.smb",filename, num_mfaces);
  std::cout <<"m3dc1_meshgen: mesh is saved in PUMI-readable \""<<mesh_filename<<"\"\n";
  mesh->writeNative(mesh_filename);

  if (num_mfaces>1000)
    sprintf(mesh_filename,"%s-%dK.vtk",filename, num_mfaces/1000);
  else
    sprintf(mesh_filename,"%s-%d.vtk",filename, num_mfaces);
  std::cout <<"m3dc1_meshgen: mesh is saved in \""<<mesh_filename<<"\" for Paraview\n";


  std::cout<<"\n<< Check the attribute \"elem\" in Paraview for the element order! >>\n";

 // attach element id
  apf::Numbering* nr=apf::numberOwnedDimension(mesh, "elem", 2);

 // attach model face ID
  apf::Numbering* ng=apf::createNumbering(mesh,"gface",apf::getConstant(2),1);
  apf::MeshIterator* it = mesh->begin(2);
  apf::MeshEntity* e;
  while ((e = mesh->iterate(it)))
  {
    assert(gmi_dim(mdl, (gmi_ent*)(mesh->toModel(e)))==2);
    number(ng,e,0,0,gmi_tag(mdl, (gmi_ent*)(mesh->toModel(e))));
  }
  mesh->end(it);

  writeVtkFiles(mesh_filename, mesh);
  destroyNumbering(nr);
  destroyNumbering(ng);

  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);

  cout<<"\n<< Continue with \"simmodeler\" for more meshing control! >>\n\n";

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  gmi_sim_stop();

 // clean-up
  delete [] meshSizes;

  GM_release(sim_model);
  M_release(sim_pmesh);
  Progress_delete(progress);

  MS_exit();
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_logOff();
  Sim_unregisterAllKeys();

#ifdef PPPL
  // for PPPL
  SimLicense_stop();
#endif

  PCU_Comm_Free();
  MPI_Finalize();
}

void get_vacuum_rgn()
{
  bdbox[0]=a_param-b_param;
  bdbox[2]=a_param+b_param;
  bdbox[1]=d_param-e_param;
  bdbox[3]=d_param+e_param;

  num_pts=numInterPts;
  interpolate_points.resize(2*num_pts);
  for (int i=0; i<num_pts; i++)
  {
    double para =2*M3DC1_PI/(num_pts-1)*i;
    double xyz[3];
    aexp(para, xyz);
    interpolate_points.at(2*i)=xyz[0];
    interpolate_points.at(2*i+1)=xyz[1];
  }

  create_vtx(&gv1_id,&(interpolate_points.at(0)));
  cout<<"create_vtx "<<gv1_id<<"\n";
  create_edge(&ge1_id, &gv1_id, &gv1_id);
  cout<<"create_edge "<<ge1_id<<" ("<<gv1_id<<", "<<gv1_id<<")\n";
  attach_periodic_cubic_curve(&ge1_id, &num_pts, &(interpolate_points.at(0)));
  int innerWallEdges[]={ge1_id};
  int num_ge=1;
  create_loop(&loop_id,&num_ge,innerWallEdges);
  cout<<"create_loop "<<loop_id<<" ("<<ge1_id<<")\n";
  set_inner_wall_boundary (&loop_id);
  ++loop_id;
}

void get_multi_rgn()
{
  FILE* fp=fopen(bdry_files[0].c_str(),"r");

  if (!fp)
  {
    cout<<" unable to find file "<<bdry_files[0]<<endl;
    return;
  }


  fscanf(fp,"%d\n",&num_pts);
  interpolate_points.resize(2*(num_pts+1));
  for(int i=0; i<num_pts; i++)
  {
      fscanf(fp,"%lf %lf\n",&x,&y);
      interpolate_points.at(2*i)=x;
      interpolate_points.at(2*i+1)=y;
      if (x<bdbox[0]) bdbox[0]=x;
      if (x>bdbox[2]) bdbox[2]=x;
      if (y<bdbox[1]) bdbox[1]=y;
      if (y>bdbox[3]) bdbox[3]=y;
  }
  fclose(fp);
  interpolate_points.at(2*num_pts) = interpolate_points[0];
  interpolate_points.at(2*num_pts+1) = interpolate_points[1];
  num_pts++;
  mid[0]=(bdbox[0]+bdbox[2])/2.+offsetX;
  mid[1]=(bdbox[1]+bdbox[3])/2.+offsetY;


    //create inner wall bdry
    create_vtx(&gv1_id,&(interpolate_points.at(0)));
    cout<<"create_vtx "<<gv1_id<<"\n";
    create_edge(&ge1_id, &gv1_id, &gv1_id);
    cout<<"create_edge "<<ge1_id<<" ("<<gv1_id<<", "<<gv1_id<<")\n";
    attach_natural_cubic_curve(&ge1_id,&num_pts,&(interpolate_points.at(0)));
    // now set the loop closed
    int innerWallEdges[]={ge1_id};
    num_ge = 1;
    create_loop(&loop_id,&num_ge,innerWallEdges);
    cout<<"create_loop "<<loop_id<<" ("<<ge1_id<<")\n";
    set_inner_wall_boundary (&loop_id);
    ++loop_id;

    inner_outer_wall_pts();
}

void inner_outer_wall_pts()
{
    for (int i=1; i<numBdry; ++i)
    {
      gv1_id++;
      ge1_id++;

      vector<double> out_pts;
      int num_out_pts;

      FILE* fp=fopen(bdry_files[i].c_str(),"r");
      assert (fp);
      fscanf(fp,"%d\n",&num_out_pts);
      out_pts.resize(2*num_out_pts);
      for (int j=0; j<num_out_pts; j++)
      {
        fscanf(fp,"%lf %lf\n",&x,&y);
        out_pts.at(2*j)=x;
        out_pts.at(2*j+1)=y;
      }
      fclose(fp);
      if (fabs(out_pts.at(0)-out_pts.at(2*(num_out_pts-1))) >1e-16 || fabs(out_pts.at(1)-out_pts.at(2*(num_out_pts-1)+1)) >1e-16)
      {
	std::cout << "First and Last Points are not same" << "\n";
	out_pts.resize(2*(num_out_pts+1));
	out_pts.at(2*num_out_pts) = out_pts[0];
	out_pts.at(2*num_out_pts+1) = out_pts[1];
	num_out_pts++; 
      }
      create_vtx(&gv1_id,&(out_pts.at(0)));
      cout<<"create_vtx "<<gv1_id<<"\n";
      create_edge(&ge1_id, &gv1_id, &gv1_id);
      cout<<"create_edge "<<ge1_id<<" ("<<gv1_id<<", "<<gv1_id<<")\n";
       
      if (num_out_pts <= 10)
      	attach_natural_cubic_curve(&ge1_id,&num_out_pts,&(out_pts.at(0)));
      else
	attach_piecewise_linear_curve(&ge1_id,&num_out_pts,&(out_pts.at(0)));

      int outerWallEdges[]={ge1_id};
      num_ge = 1;
      create_loop(&loop_id,&num_ge,outerWallEdges);
      cout<<"create_loop "<<loop_id<<" ("<<ge1_id<<")\n";
      ++loop_id;
      set_outer_wall_boundary (&loop_id);
    }  // for (int i=0; i<nout_bdryFile; ++i)
}

int make_sim_model_old (pGModel& sim_model, vector< vector<int> >& rgn_bdry)
{
  std::map< int, std::vector<int> > loopContainer;
  std::map<int, std::pair<int, int> > edgeContainer;
  std::map<int, int> edgeType;
  std::map<int, std::vector<double> > vtxContainer;
  export_model_data(vtxContainer, edgeType, edgeContainer, loopContainer);

  std::map<int, pGVertex> vertices;
  std::map<int, pGEdge> edges;
  gmi_model* model = m3dc1_model::instance()->model;
  int numL=loopContainer.size();
#ifdef PPPL // SIMMODSUITE_MAJOR_VERSION >= 15
  pGIPart part = GM_rootPart(sim_model);
#else
  pGIPart part = GM_part(sim_model);
#endif
  pGRegion outerRegion = GIP_outerRegion(part);
  vector<pGEdge> innerLoop;
  vector<pGEdge> currentLoop;
  for( std::map<int, vector<int> >:: iterator it=loopContainer.begin(); it!=loopContainer.end(); it++)
  {
    int numE=it->second.size();
    for( int i=0; i<numE; i++)
    {
      int edge=it->second[i];
      std::pair<int, int> vtx=edgeContainer[edge];
      std::vector<double>& xyz= vtxContainer.at(vtx.first);
      vertices[vtx.first] = GR_createVertex(GIP_outerRegion(part), &xyz[0]);
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
      switch (curvType)
      {
        case LIN:
          {
            double xyz1[3], xyz2[3];
            GV_point(startVert,xyz1);
            GV_point(endVert,xyz2);
            curve = SCurve_createLine(xyz1, xyz2);
          }
          break;
        case BSPLINE:
          {
            M3DC1::BSpline** data= (M3DC1::BSpline**) gmi_analytic_data(model, ae);
            vector<double> ctrlPtsX,ctrlPtsY, knots,weight;
            int order;
            data[0]->getpara(order, ctrlPtsX, knots, weight);
            data[1]->getpara(order, ctrlPtsY, knots, weight);
            int numPts=ctrlPtsX.size();
            vector<double> ctrlPts3D (3*numPts);
            for( int k=0; k<numPts; k++)
            {
              ctrlPts3D.at(3*k)=ctrlPtsX.at(k);
              ctrlPts3D.at(3*k+1)=ctrlPtsY.at(k);
              ctrlPts3D[3*k+2]=0.0;
            }
            curve = SCurve_createBSpline(order,numPts,&ctrlPts3D[0],&knots[0],NULL);
          }
          break;
        default:
          std::cout<<" curve type not support by simmetrix "<<std::endl;
          throw 1;
          break;
      }
      pGEdge pe = GR_createEdge(GIP_outerRegion(part), startVert, endVert, curve, 1);
      std::cout<<"GR_createEdge ("<<GEN_tag(startVert)<<", "<<GEN_tag(endVert)<<")\n";
      edges[edge]=pe;
      currentLoop.push_back( pe);
    }

    double corner[3]={0,0,0}, xPt[3]={1,0,0}, yPt[3]={0,1,0};
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
    GR_createFace(GIP_outerRegion(part), faceEdges.size(),&(faceEdges[0]),&(faceDirs[0]),numloops,loopDef,planarSurface,1);
        std::cout<<"GR_createFace (";
    for (int i=0; i<faceEdges.size(); ++i)
      cout<<GEN_tag(faceEdges[i])<<" ";
    cout<<")\n";
  }
  printf("Number of vertices in Simmetrix model: %d\n",GM_numVertices(sim_model));
  printf("Number of edges in Simmetrix model: %d\n",GM_numEdges(sim_model));
  printf("Number of faces in Simmetrix model: %d\n",GM_numFaces(sim_model));
  printf("Number of regions in Simmetrix model: %d\n",GM_numRegions(sim_model));
}

int make_sim_model (pGModel& sim_model, vector< vector<int> >& rgn_bdry)
{
  std::map< int, std::vector<int> > loopContainer;
  std::map<int, std::pair<int, int> > edgeContainer;
  std::map<int, int> edgeType;
  std::map<int, std::vector<double> > vtxContainer;
  export_model_data(vtxContainer, edgeType, edgeContainer, loopContainer);
  std::map<int, pGVertex> vertices;
  std::map<int, pGEdge> edges;
  std::map< int, std::vector<pGEdge> > edgesOnLoop;
  gmi_model* model = m3dc1_model::instance()->model;  
  int numL=loopContainer.size();

#ifdef PPPL // SIMMODSUITE_MAJOR_VERSION >= 15
  pGIPart part = GM_rootPart(sim_model);
#else
  pGIPart part = GM_part(sim_model);
#endif
  pGRegion outerRegion = GIP_outerRegion(part);
  // Now we'll add loops
  for( std::map<int, vector<int> >:: iterator it=loopContainer.begin(); it!=loopContainer.end(); it++)
  {
    int numE=it->second.size();
    // First we'll add the vertices on the loop
    for( int i=0; i<numE; i++)
    {
      int edge=it->second[i];
      std::pair<int, int> vtx=edgeContainer[edge];
      std::vector<double>& xyz= vtxContainer.at(vtx.first);
      vertices[vtx.first] = GR_createVertex(GIP_outerRegion(part), &xyz[0]); 
    }
    for( int j=0; j<numE; j++)
    {
      pGVertex startVert, endVert;
      pGEdge pe;
      int edge = it->second[j];
      gmi_ent* ae = gmi_find(model, 1, edge);
      int curvType=edgeType.at(edge);
      std::pair<int, int> vtx=edgeContainer[edge];
      pCurve curve=NULL;
      switch (curvType)
      {
        case LIN: //linear
          {
	    startVert=vertices[vtx.first];
            endVert=vertices[vtx.second];
            double xyz1[3], xyz2[3];
            GV_point(startVert,xyz1);
            GV_point(endVert,xyz2);
            curve = SCurve_createLine(xyz1, xyz2);
            pe = GR_createEdge(GIP_outerRegion(part), startVert, endVert, curve, 1);
            std::cout <<"LIN: GR_createEdge "<<GEN_tag(startVert)<<", "<<GEN_tag(endVert)<<"\n";
          }
          break;
        case BSPLINE: // bspline
          {
            M3DC1::BSpline** data= (M3DC1::BSpline**) gmi_analytic_data(model, ae);
            vector<double> ctrlPtsX,ctrlPtsY, knots,weight;
            int order;
            startVert=vertices[vtx.first];
            endVert=vertices[vtx.second];
            data[0]->getpara(order, ctrlPtsX, knots, weight);
            data[1]->getpara(order, ctrlPtsY, knots, weight);
            int numPts=ctrlPtsX.size();
            vector<double> ctrlPts3D (3*(numPts));
            for( int k=0; k<numPts; k++)
            {
              ctrlPts3D.at(3*k)=ctrlPtsX.at(k);
              ctrlPts3D.at(3*k+1)=ctrlPtsY.at(k);
              ctrlPts3D[3*k+2]=0.0;
            }

	    // To make it consistent, we will define every edge in counter-clockwise direction.
	    // If curve is clockwise, set edge dir to 0, otherwise 1 to follow the above convention.
	    int edgeDir = 1;
	    bool clockwise = curveOrientation(ctrlPts3D);
            if (clockwise)
		edgeDir = 0;
            curve = SCurve_createBSpline(order,numPts,&ctrlPts3D[0],&knots[0],NULL);
           if (numE == 1)
           {
           	pe = GR_createEdge(GIP_outerRegion(part), startVert, startVert, curve, edgeDir);
                std::cout <<"BSPLINE: GR_createEdge "<<GEN_tag(startVert)<<", "<<GEN_tag(startVert)<<"\n";
           }
	   else if (numE == 2)
           {
		pe = GR_createEdge(GIP_outerRegion(part), startVert, endVert, curve, edgeDir);
                std::cout <<"BSPLINE: GR_createEdge "<<GEN_tag(startVert)<<", "<<GEN_tag(endVert)<<"\n";
           }
           GM_write(sim_model,"Debug_Edges.smd",0,0);
          }
          break;
        default:
          std::cout<<" curve type not support by simmetrix "<<std::endl;
          throw 1;
          break;
      }
      edges[edge]=pe;
      edgesOnLoop[it->first].push_back(pe);
    }
  }

    // Now add the faces
    double corner[3]={0,0,0}, xPt[3]={1,0,0}, yPt[3]={0,1,0};  // the points defining the surface of the face
    vector<pGEdge> faceEdges;
    pGEdge ge;
    vector<int> faceDirs;
    int numLoops;

  for (int i=0; i<rgn_bdry.size(); ++i)
  {
    std::cout << " =========================== FACE STATS for Face # " << i+1 << "=======================================\n";
    faceEdges.clear();
    faceDirs.clear();
    numLoops = rgn_bdry[i].size();
    int loopDef[numLoops];
    std::vector <int> loopsOnFace;
    for (int idx = 0; idx < numLoops; ++idx)
	loopsOnFace.push_back(rgn_bdry[i][idx]);
    int outerMostLoop = outerLoop(loopsOnFace, edgesOnLoop);
    loopsOnFace.clear();
    std::cout << "OuterMost Loop on this face = " << outerMostLoop << "\n";
    for (int idx=0; idx<rgn_bdry[i].size(); ++idx)
    {
	
        int numEdgesOnLoop = edgesOnLoop.at(rgn_bdry[i][idx]).size();
        loopDef[idx] = faceEdges.size();
        double dirs[numEdgesOnLoop];
	int edgeDir = 0;
	if (rgn_bdry[i][idx] == outerMostLoop)
		edgeDir = 1;
	for (int nEdge = 0; nEdge < numEdgesOnLoop; ++nEdge)
	{
     		ge = edgesOnLoop.at(rgn_bdry[i][idx])[nEdge];
        	faceEdges.push_back(ge);
		faceDirs.push_back(edgeDir);
	}
    } 
    pSurface planarSurface;
    planarSurface = SSurface_createPlane(corner,xPt,yPt);
    GM_write(sim_model,"Debug_Faces.smd",0,0);
    std::cout << "Number of Edges on Face  = " <<   faceEdges.size() << "\n";
    for (int f =0; f < faceEdges.size(); ++f)
 	std::cout << "Edge # " << f+1 << " length = " << GE_length(faceEdges[f],0,1) << "\n";
    std::cout << "Number of Loops = " << numLoops << "\n";
    for (int mk = 0; mk < numLoops; ++mk)
    	std::cout << "Dir = " << faceDirs[mk] << "\n";
    for (int mk = 0; mk < numLoops; ++mk)
    	std::cout << "Loop Def = " << loopDef[mk] << "\n"; 
    pGFace gf = GR_createFace(GIP_outerRegion(part), faceEdges.size(),
                  &(faceEdges[0]),&(faceDirs[0]),numLoops,&(loopDef[0]),planarSurface,1);
    std::cout<<"GR_createFace (";
    for (int j=0; j<faceEdges.size(); ++j)
      cout<<GEN_tag(faceEdges[j])<<" ";
    cout<<")\n";
    GEN_setNativeIntAttribute(gf, i+1, "faceType");
  }
  
  // The starting point of each curves may not be perfectly aligned. To align the starting points
  // on two consective loops can be done using the routine "setParallelVertices". This will be needed to set layered
  // mesh on finite thickness wall face.

  int innerLoop = 4;
  int outerLoop = 5;
  //setParallelVertices(innerLoop, outerLoop, edgesOnLoop);
   
  printf("Number of vertices in Simmetrix model: %d\n",GM_numVertices(sim_model));
  printf("Number of edges in Simmetrix model: %d\n",GM_numEdges(sim_model));
  printf("Number of faces in Simmetrix model: %d\n",GM_numFaces(sim_model));
  printf("Number of regions in Simmetrix model: %d\n",GM_numRegions(sim_model));
}

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

// Check the orientation of a curve. The method is applicable 
// to non-convex polygons too. 
// Returns true if the curve is clockwise, false if counter-clockwise.
bool curveOrientation(std::vector <double> &curvePts)
{
	int numPts = curvePts.size()/3;
	double sumEdges = 0.0;
	for (int i = 1; i < numPts; ++i)
	{
		double x1 = curvePts[(i-1)*3];
		double y1 = curvePts[(i-1)*3+1];
		double x2 = curvePts[i*3];
		double y2 = curvePts[i*3+1];	

		double areaUnderEdge = (x2 - x1)* (y2 + y1);
		sumEdges = sumEdges + areaUnderEdge; 
	}
	
	if (sumEdges > 0)
		return true;
	
	return false;
}

int outerLoop(std::vector <int> loops, std::map< int, std::vector<pGEdge> > edgesOnLoop)
{
	if (loops.size() == 1)
		return loops[0];

	int index = -1;
	double minDist = 1e16;
	double pt[3] = 	{10000000.0,0.0,0.0};
	for (int i = 0; i < loops.size(); ++i)
	{
		int indexTemp = loops[i];
		int numEdgesOnLoop = edgesOnLoop.at(indexTemp).size();
		for (int j =0; j < numEdgesOnLoop; ++j)
		{
			pGEdge ge = edgesOnLoop.at(indexTemp)[j];
			double outPt[3], outPar;
			GE_closestPoint(ge, pt, outPt, &outPar);
			double dist = getDist2D (pt, outPt);
			assert (outPt[0] < pt[0]);
			if (dist < minDist)
			{
				minDist = dist;
				index = indexTemp;
			}
		}		
	}	
	return index;
}


// The starting point of each curves may not be perfectly aligned. To align the starting points
// on two consective loops can be done using the given routine. This will be needed to set layered
// mesh on finite thickness wall face.
void setParallelVertices(int loopInner, int loopOuter, std::map< int, std::vector<pGEdge> > edgesOnLoop)
{
  	pGEdge ge1 = edgesOnLoop.at(loopInner)[0];
  	pGEdge ge2 = edgesOnLoop.at(loopOuter)[0];
  	pGVertex v1 = GE_vertex(ge1,0);
  	double vPt1[3];
  	GV_point(v1, vPt1);
  	double par;
  	GE_closestPoint(ge2, vPt1, 0, &par);
  	pGVertex newV = GE_insertVertex(ge2, par);
  	pPList gves = GV_edges(newV);
  	pGEdge geNew0 = static_cast<pGEdge>(PList_item(gves, 0));
  	pGEdge geNew1 = static_cast<pGEdge>(PList_item(gves, 1));
  	if (GE_vertex(geNew0,1) == newV)
        	pGEdge finalEdge = GM_combineEdges (geNew0, geNew1);
  	else	
        	pGEdge finalEdge = GM_combineEdges (geNew1, geNew0);
}
