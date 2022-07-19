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

char meshSizeFnName[]="meshSizeFn";
double meshSizeFn(const double gpt[3], void* data);

int make_sim_model (pGModel& sim_model, vector< vector<int> >& rgn_bdry);

void messageHandler(int type, const char *msg);
using namespace std;
double thickness=0.02, width=2.5, height=4.;
double offsetX=0, offsetY=0;
char modelName[256];
char in_bdryFile[256];
int numBdry=0, nout_bdryFile;
int numRgn=1;
char bdryFile[128];
vector<string> bdry_files;
char out_bdryFile[256];
int modelType=0;
int reorder=0;
int useVacuumParams=0;
int useOriginal = 1;
int numInterPts =20;
double meshSizes[] = {0.05, 0.05, 0.05};
int useAnaltyicSize=0;
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
int one=1, two=2, three=3, four=4;
int num_pts;
double mid[]={0.,0.};
double pt_pre[2], mergeTol;
    double dir1[2], dir2[2], normal1[2],normal2[2];
    double paras[]={0.,1.};
    double ctr_pts_o[12];


vector<double> interpolate_points;

void inner_wall_pts(  int num_pts);

void inner_outer_wall_pts();

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();

  if (argc!=2) 
  {    
    cout<<"[m3dc1_meshgen ERROR] usage: ./m3dc1_meshgen input_filename\n"; 
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
  vector< vector<int> > rgn_bdry;
  int index=0, numLoop, loopID;  

  while(EOF!=fscanf(infp,"%s",namebuff))
  {
    if (strcmp(namebuff,"numBdry")==0) 
    {
      fscanf(infp,"%d",&numBdry);
      cout<<"# bdry files = "<<numBdry<<"\n";
    }
    if (strcmp(namebuff,"numRgn")==0) 
    {
      fscanf(infp,"%d",&numRgn); 
      rgn_bdry.resize(numRgn);
      cout<<"# rgn = "<<rgn_bdry.size()<<"\n";
    }
    if (strcmp(namebuff,"bdryFile")==0) 
    {
      fscanf(infp, "%s", bdryFile);
      bdry_files.push_back(bdryFile);
    }

   if (strcmp(namebuff,"rgnBdry")==0)
   {
     fscanf(infp,"%d",&numLoop);  
     cout<<"index "<<index<<", numLoop "<<numLoop<<"\n";
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
      for(int i=0; i<3; i++)
        fscanf(infp, "%lf",meshSizes+i);
    }
    if (strcmp(namebuff,"useVacuumParams")==0) fscanf(infp, "%d",&useVacuumParams);
    if (strcmp(namebuff,"vacuumParams")==0) 
      fscanf(infp, "%lf %lf %lf %lf %lf\n", &a_param, &b_param, &c_param, &d_param, &e_param);
    if (strcmp(namebuff,"meshGradationRate")==0) fscanf(infp, "%lf",&meshGradationRate);
  }
  fclose(infp);
  if (numBdry==0)
    numBdry=bdry_files.size();
  else
    assert(numBdry==bdry_files.size());

  for (vector<string>::iterator sit=bdry_files.begin(); sit!=bdry_files.end(); ++sit)
    cout<<"bdry file "<<*sit<<"\n";

  for (int i=0; i<rgn_bdry.size(); ++i)
  {
    cout<<"Loop ID's for face "<<i<<": ";
    for (int j=0; j<rgn_bdry[i].size(); ++j)
    {
      cout<<rgn_bdry[i][j]<<" "; 
    }
    cout<<"\n";
  }

  mergeTol=0.3*thickness;

  // *************************
  //     interpolate points
  // *************************
  FILE* fp=fopen(bdry_files[0].c_str(),"r");

  if (!fp) 
  {
    cout<<" unable to find file "<<bdry_files[0]<<endl;
    return 1;
  }

  fscanf(fp,"%d\n",&num_pts);
  interpolate_points.resize(2*num_pts);
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
  mid[0]=(bdbox[0]+bdbox[2])/2.+offsetX;
  mid[1]=(bdbox[1]+bdbox[3])/2.+offsetY;

//  create inner wall bdry
    create_vtx(&gv1_id,&(interpolate_points.at(0)));
    create_vtx(&gv2_id,&(*(interpolate_points.rbegin()+1)));
    create_edge(&ge1_id, &gv1_id, &gv2_id);
    create_edge(&ge2_id, &gv2_id, &gv1_id);
    attach_natural_cubic_curve(&ge1_id,&num_pts,&(interpolate_points.at(0)));
    // now set the loop closed
    ctr_pts[0]=*(interpolate_points.rbegin()+1);
    ctr_pts[1]=*interpolate_points.rbegin();
    ctr_pts[10]=interpolate_points.at(0);
    ctr_pts[11]=interpolate_points.at(1);
    eval_normal(&ge1_id, paras+1, normal1); 
    eval_normal(&ge1_id, paras, normal2);
    dir1[0]=-normal1[1];
    dir1[1]=normal1[0];
    dir2[0]=-normal2[1];
    dir2[1]=normal2[0];
    double len=getDist2D(&(interpolate_points.at(0)),&(*(interpolate_points.rbegin()+1)));
    double offset_distance[4]={0.2*len,0.4*len,-0.2*len,-0.2*len};
    offset_point(ctr_pts, dir1, offset_distance, ctr_pts+2);
    offset_point(ctr_pts, dir1, offset_distance+1, ctr_pts+4);
    offset_point(ctr_pts+10, dir2, offset_distance+2, ctr_pts+6);
    offset_point(ctr_pts+10, dir2, offset_distance+3, ctr_pts+8);
    attach_b_spline_curve(&ge2_id,&order_p,&numCtrPts,ctr_pts,knots,NULL);
    int innerWallEdges[]={ge1_id,ge2_id};
    create_loop(&loop_id,&num_ge,innerWallEdges);
    set_inner_wall_boundary (&loop_id);
    cout<<"Loop "<<loop_id<<" created with GE "<<ge1_id<<" "<< ge2_id<<"\n";
    ++loop_id; 

  if (modelType==3) // &&  && thickness>0)
    inner_wall_pts(num_pts);
  else
    inner_outer_wall_pts();
  
    // increase ID's for geometric vertices and edges
    gv1_id=gv2_id+1;
    gv2_id=gv1_id+1;
    ge1_id=ge2_id+1;
    ge2_id=ge1_id+1;

    if (!useVacuumParams) // #edges=6
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
        create_vtx(&gv2_id,right);
        create_edge(&ge1_id, &gv1_id, &gv2_id);
        create_edge(&ge2_id, &gv2_id, &gv1_id);
        attach_b_spline_curve(&ge1_id,&order_p,&numCtrPts,ctrPtsBottom,knots,NULL);
        attach_b_spline_curve(&ge2_id,&order_p,&numCtrPts,ctrPtsTop,knots,NULL);
        int vacuumEdges[]={ge1_id,ge2_id};
        create_loop(&loop_id,&num_ge,vacuumEdges);
        cout<<"Loop "<<loop_id<<" created with GE "<<ge1_id<<" "<< ge2_id<<"\n";
        set_vacuum_boundary (&loop_id);
         ++loop_id;
    }
    else // #edges=5
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
      create_edge(&ge1_id, &gv1_id, &gv1_id);
      attach_periodic_cubic_curve(&ge1_id,&num_pts,&(interpolate_points_vacuum.at(0)));
      int vacuumEdges[]={ge1_id};
      num_ge=1;
      create_loop(&loop_id,&num_ge,vacuumEdges);
      cout<<"Loop "<<loop_id<<" created with GE "<<ge1_id<<"\n";
      set_vacuum_boundary (&loop_id);
      ++loop_id;
    }

  int numplot=5000;

  int numEdge;
  if (useVacuumParams)
    numEdge=numBdry*2+1;
  else
    numEdge=numBdry*2+2;
 
  cout<<__func__<<": simLic = "<<simLic<<", #edge = numEdge\n";
  cout<<"model #v "<<m3dc1_model::instance()->model->n[0]
      <<" #e "<<m3dc1_model::instance()->model->n[1] 
      <<" #f "<<m3dc1_model::instance()->model->n[2]<<"\n";

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
#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_sim.h>
#include "gmi_mesh.h"
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <ma.h>
#include <cstdlib>
#include <iostream>
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

  make_sim_model (sim_model ,rgn_bdry);
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
    if (!iface&&useAnaltyicSize)
    {
      MS_registerSizeExprFunc(meshSizeFnName, meshSizeFn, NULL);
      MS_setMeshSize(meshCase,face,1,0,"meshSizeFn($x,$y,$z)");
      iface++;
    }
    else
    {
      MS_setMeshSize(meshCase,face,1, meshSizes[iface++],NULL);
      cout<<"set mesh size of face "<<iface-1<<" to "<<meshSizes[iface-1]<<"\n";
    }
  }
  cout<<"# gface "<<iface<<"\n";
//assert(iface<=3);
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

  apf::Numbering* nr=apf::numberOwnedDimension(mesh, "elem", 2);
  writeVtkFiles(mesh_filename, mesh);
  destroyNumbering(nr);

  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);

  cout<<"\n<< Continue with \"simmodeler\" for more meshing control! >>\n\n";

  mesh->destroyNative();
  apf::destroyMesh(mesh);
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

#ifdef PPPL
  // for PPPL
  SimLicense_stop();
#endif

  PCU_Comm_Free();
  MPI_Finalize();
}

void inner_wall_pts(  int num_pts)
{
/*
    int vtxMerged=0, vtxInserted=0;
      vector<double> interpolate_points_o;
      double increment=1.0/double(num_pts-1);
      double dir[2], dir_pre[2];
      for(int i=0; i<interpolate_points.size()/2; i++)
      { 
        double para=i*increment; 
        eval_normal(&one, &para, dir);
        double pt_new[2];
        offset_point(&(interpolate_points.at(2*i)), dir, &thickness, pt_new);
        if (i!=0)
        { 
          double crossProduct=dir_pre[0]*dir[1]-dir_pre[1]*dir[0];
          double dotProduct=dir_pre[0]*dir[0]+dir_pre[1]*dir[1];
          if (crossProduct<0&&fabs(pt_new[0]-pt_pre[0])+fabs(pt_new[1]-pt_pre[1])<mergeTol)
          { 
            vtxMerged++;
            double pt_mid[2], dir_mid[2];
            double para_mid=(i-0.5)*increment;
            eval_position(&one, &para_mid, pt_mid);
            eval_normal(&one, &para_mid, dir_mid);
            offset_point(pt_mid, dir_mid, &thickness, pt_new);
            interpolate_points_o.pop_back();
            interpolate_points_o.pop_back();
          }
          else
          { 
            double tol_curv=0.3;
            if (dotProduct<tol_curv)
            { 
              double para_mid=(i-0.5)*increment;
              double pt_mid[2], dir_mid[2], pt_mid_new[2];
              eval_position(&one, &para_mid, pt_mid);
              eval_normal(&one, &para_mid, dir_mid);
              offset_point(pt_mid, dir_mid, &thickness, pt_mid_new);
              interpolate_points_o.push_back(pt_mid_new[0]);
              interpolate_points_o.push_back(pt_mid_new[1]);
              vtxInserted++;
            }
          }
        }
        interpolate_points_o.push_back(pt_new[0]);
        interpolate_points_o.push_back(pt_new[1]);
        pt_pre[0]=pt_new[0];
        pt_pre[1]=pt_new[1];
        dir_pre[0]=dir[0];
        dir_pre[1]=dir[1];
      }
      cout<<"[info] num pt merged "<<vtxMerged<<" num pt inserted "<< vtxInserted<<endl;

      FILE * fp_t=fopen("offset.txt","w");
      for(int i=0; i<interpolate_points_o.size()/2; i++)
      {
        fprintf(fp_t,"%lf %lf\n",interpolate_points_o.at(2*i),interpolate_points_o.at(2*i+1));
      }
      fclose(fp_t);

      for(int i=0; i<3; i++)
      {
        offset_point(ctr_pts+2*i, normal1, &thickness, ctr_pts_o+2*i);
        offset_point(ctr_pts+2*(i+3), normal2, &thickness, ctr_pts_o+2*(i+3));
      }
      create_vtx(&three,&(interpolate_points_o.at(0)));
      create_vtx(&four,&(*(interpolate_points_o.rbegin()+1)));
      create_edge(&three, &three, &four);
      create_edge(&four, &four, &three);
      int numpt_new=interpolate_points_o.size()/2;
      attach_natural_cubic_curve(&three,&numpt_new,&(interpolate_points_o.at(0)));
      attach_b_spline_curve(&four,&order_p,&numCtrPts,ctr_pts_o,knots,NULL);
      int outerWallEdges[]={three,four};
      create_loop(&two,&two,outerWallEdges);
      set_outer_wall_boundary (&two);
    }
*/
}

void inner_outer_wall_pts()
{

    for (int i=1; i<numBdry; ++i)
    {
      gv1_id=gv2_id+1;
      gv2_id=gv1_id+1;
      ge1_id=ge2_id+1;
      ge2_id=ge1_id+1;

      vector<double> out_pts;
      int num_out_pts;

      FILE* fp=fopen(bdry_files[i].c_str(),"r");
      cout<<"Open "<<bdry_files[i].c_str()<<"\n";
      assert (fp);
      fscanf(fp,"%d\n",&num_out_pts);
      out_pts.resize(2*num_out_pts);
      for (int i=0; i<num_out_pts; i++)
      {
        fscanf(fp,"%lf %lf\n",&x,&y);
        out_pts.at(2*i)=x;
        out_pts.at(2*i+1)=y;
      }
      fclose(fp);

      create_vtx(&gv1_id,&(out_pts.at(0)));
      create_vtx(&gv2_id,&(*(out_pts.rbegin()+1)));
      create_edge(&ge1_id, &gv1_id, &gv2_id);
      create_edge(&ge2_id, &gv2_id, &gv1_id);

      attach_natural_cubic_curve(&ge1_id,&num_out_pts,&(out_pts.at(0)));

      ctr_pts_o[0]=*(out_pts.rbegin()+1);
      ctr_pts_o[1]=*out_pts.rbegin();
      ctr_pts_o[10]=out_pts.at(0);
      ctr_pts_o[11]=out_pts.at(1);
      double dir1[2], dir2[2], normal1[2], normal2[2];
      double paras[]={0.,1.};
      eval_normal(&ge1_id, paras+1, normal1);
      eval_normal(&ge1_id, paras, normal2);
      dir1[0]=-normal1[1];
      dir1[1]=normal1[0];
      dir2[0]=-normal2[1];
      dir2[1]=normal2[0];
      double len=getDist2D(&(out_pts.at(0)),&(*(out_pts.rbegin()+1)));
      double offset_distance[4]={0.2*len,0.4*len,-0.2*len,-0.2*len};
      offset_point(ctr_pts_o, dir1, offset_distance, ctr_pts_o+2);
      offset_point(ctr_pts_o, dir1, offset_distance+1, ctr_pts_o+4);
      offset_point(ctr_pts_o+10, dir2, offset_distance+2, ctr_pts_o+6);
      offset_point(ctr_pts_o+10, dir2, offset_distance+3, ctr_pts_o+8);

      attach_b_spline_curve(&ge2_id,&order_p,&numCtrPts,ctr_pts_o,knots,NULL);
      int outerWallEdges[]={ge1_id,ge2_id};
      create_loop(&loop_id,&num_ge,outerWallEdges);
      cout<<"Loop "<<loop_id<<" created with GE "<<ge1_id<<" "<< ge2_id<<"\n";
      ++loop_id;
      set_outer_wall_boundary (&loop_id);
    }  // for (int i=0; i<nout_bdryFile; ++i)
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
      // GIP_insertVertexInRegion(part,&xyz[0],outerRegion);
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
        case LIN: //linear
          {
            double xyz1[3], xyz2[3];
            GV_point(startVert,xyz1);
            GV_point(endVert,xyz2);
            curve = SCurve_createLine(xyz1, xyz2);
          }
          break;
        case BSPLINE: // bspline
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
      cout<<"Create GE ["<<edge<<"] "<<GEN_tag(pe)<<" vtx "<<GEN_tag(startVert)<<", "
          <<GEN_tag(endVert)<<"\n";
      // GIP_insertEdgeInRegion(part, startVert, endVert, curve, 1, outerRegion);
      edges[edge]=pe;
    }
  }

  for (std::map<int, pGEdge>::iterator mit=edges.begin(); mit!=edges.end(); ++mit)
    cout<<"edges["<<mit->first<<"] = "<<GEN_tag(mit->second)<<"\n";

    // Now add the faces
    double corner[3]={0,0,0}, xPt[3]={1,0,0}, yPt[3]={0,1,0};  // the points defining the surface of the face
    vector<pGEdge> faceEdges;
  pGEdge ge;
    vector<int> faceDirs;
    int loopDef[2] = {0,0};
    int numLoops;
  vector<int>::iterator max_loop_it;
  vector<int>::iterator min_loop_it;

  for (int i=0; i<rgn_bdry.size(); ++i)
  {
    faceEdges.clear();
    faceDirs.clear();
    numLoops = (rgn_bdry[i].size()>1)?2:1;
    loopDef[1] = (numLoops==1)?0:2; //(numLoops-1)*2;
    cout<<"#loop "<<numLoops<<" loop_id "<<rgn_bdry[i][0]<<"\n";
      // find outmost loop id
    max_loop_it = std::max_element(rgn_bdry[i].begin(), rgn_bdry[i].end());
    min_loop_it = std::min_element(rgn_bdry[i].begin(), rgn_bdry[i].end());
    cout<<"#loop "<<rgn_bdry[i].size()<<" loop_id "<<rgn_bdry[i][0]
        <<" min_loop_id "<<*min_loop_it<<", max_loop_id "<<*max_loop_it<<"\n";
    for (int idx=0; idx<rgn_bdry[i].size(); ++idx)
    {  
       cout<<"loop_id "<<rgn_bdry[i][idx]<<"\n";
      if (rgn_bdry[i][idx]<*max_loop_it && rgn_bdry[i][idx]==*min_loop_it)
      {
        for (int j=0; j<=1; ++j)
          {
            ge = edges[rgn_bdry[i][idx]*2-j];
            faceEdges.push_back(ge);
            faceDirs.push_back(0);
            cout<<"inner loop "<<rgn_bdry[i][idx]*2-j<<" faceEdges.push_back("<<GEN_tag(ge)<<" dir 0\n"; 
          }
      }
      else if (rgn_bdry[i][idx]==*max_loop_it) // outer wall
        for (int j=1; j>=0; --j) 
        {
          ge = edges[rgn_bdry[i][idx]*2-j];
          faceEdges.push_back(ge);
          faceDirs.push_back(1);
          cout<<"outer loop "<<rgn_bdry[i][idx]*2-j<<" faceEdges.push_back("<<GEN_tag(ge)<<" dir 1\n";
        }
     } // for (int idx=0; idx<numLoops; ++idx)

    pSurface planarSurface;
    planarSurface = SSurface_createPlane(corner,xPt,yPt);
    for (int i=0; i<faceEdges.size(); ++i)
      cout<<"faceEdges["<<i<<"] = "<<GEN_tag(faceEdges[i])<<", dir "<<faceDirs[i]
          <<", #loops "<<numLoops <<", loopDef "<<loopDef[0]<<", "<<loopDef[1]<<"\n";
      
    GR_createFace(GIP_outerRegion(part), faceEdges.size(),
                  &(faceEdges[0]),&(faceDirs[0]),numLoops,loopDef,planarSurface,1);
    cout<<"+++ create gf +++\n";
    //GIP_insertFaceInRegion(part,faceEdges.size(),&(faceEdges[0]),&(faceDirs[0]),numloops,loopDef,planarSurface,1,outerRegion);
  }
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

// copied from meshSizeFn.h

double meshSizeFn(const double gpt[3], void* data)
{
  double r0=(bdbox[0]+bdbox[2])/2.;
  double z0=(bdbox[1]+bdbox[3])/2.;
  double rLen=(bdbox[2]-bdbox[0])/2.5;
  double zLen=(bdbox[3]-bdbox[1])/2.5;
  double r=gpt[0];
  double z=gpt[1];
 
  double fine_size=meshSizes[0];
  double coarse_size_axis=meshSizes[0];
  double coarse_size_wall=2*meshSizes[0];
  double val = (r-r0)*(r-r0)/rLen/rLen+(z-z0)*(z-z0)/zLen/zLen;
 
  double meshSize;
  double alpha=0.15; // ~r^2 at psi_0
  //cout<<val<<endl;
  if (val<alpha) meshSize=fine_size+coarse_size_axis*(alpha-val);
  else meshSize=fine_size+coarse_size_wall*(val-alpha);  
  //cout<<"meshSize "<<meshSize<<endl;
  return meshSize;
}
