/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <ReducedQuintic.h>
#include <ReducedQuinticImplicit.h>
#include <m3dc1_model.h>
#include "gmi_null.h" 
#include "m3dc1_scorec.h"
#include <iostream>
#include <assert.h>
#include "apf.h"
#include "apfMDS.h"
#include "petscksp.h"
#include "PCU.h"
#include "m3dc1_mesh.h" 
#include "apfMesh.h" 
using namespace apf;
using namespace std;
#define zeroValue 1e-6
#define numInter 12
#define numPlot 500
double eta =.001;     
double mu = .005;
double dt = 0.1;
char model[512];
char fieldName[512];
char mesh1[512];
char mesh2[512];
int compareU = 0;
int comparePsi = 0;
const double plotLen=0.2;
double end[]={1.1, 3.2};
double beta_12[] =
 {0.249286745170910, 0.501426509658179, 0.249286745170910, 0.063089014491502, 
  0.873821971016996, 0.063089014491502, 0.310352351033784, 0.636502499121399, 
  0.053145049844817, 0.636502499121399, 0.053145049844817, 0.310352351033784}; 

double alpha_12[] =
  {0.249286745170910, 0.249286745170910, 0.501426509658179, 0.063089014491502, 
  0.063089014491502, 0.873821971016996, 0.636502499121399, 0.053145049844817, 
  0.310352351033784, 0.310352351033784, 0.636502499121399, 0.053145049844817}; 

double area_weight_12[] =
  {0.116786275726379, 0.116786275726379, 0.116786275726379, 0.050844906370207, 
  0.050844906370207, 0.050844906370207, 0.082851075618374, 0.082851075618374, 
  0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374};
double alpha_79 [] =
 {0.333333333333333,-0.001900928704400, 0.500950464352200, 0.500950464352200, 
  0.023574084130543, 0.488212957934729, 0.488212957934729, 0.089726636099435, 
  0.455136681950283, 0.455136681950283, 0.196007481363421, 0.401996259318289, 
  0.401996259318289, 0.488214180481157, 0.255892909759421, 0.255892909759421, 
  0.647023488009788, 0.176488255995106, 0.176488255995106, 0.791658289326483, 
  0.104170855336758, 0.104170855336758, 0.893862072318140, 0.053068963840930, 
  0.053068963840930, 0.916762569607942, 0.041618715196029, 0.041618715196029, 
  0.976836157186356, 0.011581921406822, 0.011581921406822, 0.048741583664839, 
  0.606402646106160, 0.344855770229001, 0.606402646106160, 0.048741583664839, 
  0.344855770229001, 0.006314115948605, 0.615842614456541, 0.377843269594854, 
  0.615842614456541, 0.006314115948605, 0.377843269594854, 0.134316520547348, 
  0.559048000390295, 0.306635479062357, 0.559048000390295, 0.134316520547348, 
  0.306635479062357, 0.013973893962392, 0.736606743262866, 0.249419362774742, 
  0.736606743262866, 0.013973893962392, 0.249419362774742, 0.075549132909764, 
  0.711675142287434, 0.212775724802802, 0.711675142287434, 0.075549132909764, 
  0.212775724802802,-0.008368153208227, 0.861402717154987, 0.146965436053239, 
  0.861402717154987,-0.008368153208227, 0.146965436053239, 0.026686063258714, 
  0.835586957912363, 0.137726978828923, 0.835586957912363, 0.026686063258714, 
  0.137726978828923, 0.010547719294141, 0.929756171556853, 0.059696109149007, 
  0.929756171556853, 0.010547719294141, 0.059696109149007}; 
double beta_79[] =
 {0.333333333333333, 0.500950464352200,-0.001900928704400, 0.500950464352200, 
  0.488212957934729, 0.023574084130543, 0.488212957934729, 0.455136681950283, 
  0.089726636099435, 0.455136681950283, 0.401996259318289, 0.196007481363421, 
  0.401996259318289, 0.255892909759421, 0.488214180481157, 0.255892909759421,
  0.176488255995106, 0.647023488009788, 0.176488255995106, 0.104170855336758, 
  0.791658289326483, 0.104170855336758, 0.053068963840930, 0.893862072318140, 
  0.053068963840930, 0.041618715196029, 0.916762569607942, 0.041618715196029, 
  0.011581921406822, 0.976836157186356, 0.011581921406822, 0.344855770229001, 
  0.048741583664839, 0.606402646106160, 0.344855770229001, 0.606402646106160, 
  0.048741583664839, 0.377843269594854, 0.006314115948605, 0.615842614456541, 
  0.377843269594854, 0.615842614456541, 0.006314115948605, 0.306635479062357, 
  0.134316520547348, 0.559048000390295, 0.306635479062357, 0.559048000390295, 
  0.134316520547348, 0.249419362774742, 0.013973893962392, 0.736606743262866, 
  0.249419362774742, 0.736606743262866, 0.013973893962392, 0.212775724802802, 
  0.075549132909764, 0.711675142287434, 0.212775724802802, 0.711675142287434, 
  0.075549132909764, 0.146965436053239,-0.008368153208227, 0.861402717154987, 
  0.146965436053239, 0.861402717154987,-0.008368153208227, 0.137726978828923, 
  0.026686063258714, 0.835586957912363, 0.137726978828923, 0.835586957912363, 
  0.026686063258714, 0.059696109149007, 0.010547719294141, 0.929756171556853, 
  0.059696109149007, 0.929756171556853, 0.010547719294141 };
double area_weight_79 [] =
 {0.033057055541624, 0.000867019185663, 0.000867019185663, 0.000867019185663, 
  0.011660052716448, 0.011660052716448, 0.011660052716448, 0.022876936356421, 
  0.022876936356421, 0.022876936356421, 0.030448982673938, 0.030448982673938, 
  0.030448982673938, 0.030624891725355, 0.030624891725355, 0.030624891725355, 
  0.024368057676800, 0.024368057676800, 0.024368057676800, 0.015997432032024, 
  0.015997432032024, 0.015997432032024, 0.007698301815602, 0.007698301815602, 
  0.007698301815602,-0.000632060497488,-0.000632060497488,-0.000632060497488, 
  0.001751134301193, 0.001751134301193, 0.001751134301193, 0.016465839189576, 
  0.016465839189576, 0.016465839189576, 0.016465839189576, 0.016465839189576, 
  0.016465839189576, 0.004839033540485, 0.004839033540485, 0.004839033540485, 
  0.004839033540485, 0.004839033540485, 0.004839033540485, 0.025804906534650, 
  0.025804906534650, 0.025804906534650, 0.025804906534650, 0.025804906534650, 
  0.025804906534650, 0.008471091054441, 0.008471091054441, 0.008471091054441, 
  0.008471091054441, 0.008471091054441, 0.008471091054441, 0.018354914106280, 
  0.018354914106280, 0.018354914106280, 0.018354914106280, 0.018354914106280, 
  0.018354914106280, 0.000704404677908, 0.000704404677908, 0.000704404677908, 
  0.000704404677908, 0.000704404677908, 0.000704404677908, 0.010112684927462, 
  0.010112684927462, 0.010112684927462, 0.010112684927462, 0.010112684927462, 
  0.010112684927462, 0.003573909385950, 0.003573909385950, 0.003573909385950, 
  0.003573909385950, 0.003573909385950, 0.003573909385950}; 

double * alpha, * beta, * area_weight;

int xy2area (double coord[3][2], double xy[2], double areaCoord[2])
{
  double jacobi = (coord[1][1]-coord[2][1])*(coord[0][0]-coord[2][0]) + (coord[2][0]-coord[1][0])*(coord[0][1]-coord[2][1]);
  areaCoord[0]=(coord[1][1]-coord[2][1])*(xy[0]-coord[2][0])+(coord[2][0]-coord[1][0])*(xy[1]-coord[2][1]);
  areaCoord[0] /= jacobi;
  areaCoord[1] = (coord[2][1]-coord[0][1])*(xy[0]-coord[2][0])+(coord[0][0]-coord[2][0])*(xy[1]-coord[2][1]);
  areaCoord[1] /= jacobi;
  if(areaCoord[0]>-zeroValue && areaCoord[1]>-zeroValue&& 1-areaCoord[0]-areaCoord[1]>-zeroValue)
  return 1;
  else return 0;
}

void area2xy (double coord[3][2], double areaCoord[2], double xy[2])
{
  for(int i=0; i<2; i++)
    xy[i]=areaCoord[0]*coord[0][i]+areaCoord[1]*coord[1][i]+(1-areaCoord[0]-areaCoord[1])*coord[2][i];
}
double dofsBuff[512];
int evalField (apf::Mesh2* mesh, apf::Field* field,double xy[2], double * U, double *psi)
{
  const int dim=2;
  apf::MeshEntity* e;
  apf::MeshIterator* it = mesh->begin(dim);
  Vector3 vec;
  double coord[3][2];
  Downward down_ent;
  while ((e = mesh->iterate(it)))
  {
    int num_down_ent =  mesh->getDownward(e, 0, down_ent);
    assert(num_down_ent==3);
    for (int i=0; i<num_down_ent; i++)
    {
      mesh->getPoint(down_ent[i], 0, vec);
      coord[i][0]=vec[0];
      coord[i][1]=vec[1];
    }
    double areaCoord[2];
    if(xy2area(coord,xy, areaCoord))
    {
      ReducedQuinticImplicit shape;
      shape.setCoord(coord);
      double psiDofs[18], UDofs[18];
      for(int i=0; i<num_down_ent; i++)
      {
        apf::getComponents(field,down_ent[i],0,dofsBuff);
        for(int j=0; j<6; j++)
        {
          UDofs[i*6+j] = dofsBuff[j];
          psiDofs[i*6+j] = dofsBuff[6+j]; 
        }
      }
      shape.setDofs(UDofs);
      shape.eval_g(xy,U);
      shape.setDofs(psiDofs);
      shape.eval_g(xy,psi);
      mesh->end(it);
      return 1;
    }
  }
  mesh->end(it);
  return 0;
}
void printMeshInfo(apf::Mesh2 * mesh)
{
  apf::DynamicArray<apf::MeshTag*> tags;
  mesh->getTags(tags);
  for(int i=0; i<tags.getSize(); i++)
  {
    cout<<" mesh ith tag "<<mesh->getTagName(tags[i])<<endl;
  }
}
int main(int argc, char* argv[])
{
  cout<<"\t Usage: ./main model mesh1 mesh2 fieldName1 fieldName2 compareU comparePsi"<<endl;
  MPI_Init(&argc,&argv);
  if(numInter==12)
  {
    alpha = alpha_12;
    beta = beta_12;
    area_weight = area_weight_12;
  }
  else
  {
    alpha = alpha_79; 
    beta = beta_79; 
    area_weight = area_weight_79;
  }

  FILE* infp= fopen("input","r");
  if(!infp)
  {
    cout<<"File input not found"<<endl;
    return 1;
  }
  char namebuff[512];
  while(EOF!=fscanf(infp,"%s",namebuff))
  {
    if(strcmp(namebuff,"mu")==0) fscanf(infp,"%lf",&mu);
    if(strcmp(namebuff,"eta")==0) fscanf(infp,"%lf",&eta);
    if(strcmp(namebuff,"dt")==0) fscanf(infp,"%lf",&dt);
    if(strcmp(namebuff,"model")==0) fscanf(infp,"%s",model);
    if(strcmp(namebuff,"fieldName")==0) fscanf(infp,"%s",fieldName);
    if(strcmp(namebuff,"mesh1")==0) fscanf(infp,"%s",mesh1);
    if(strcmp(namebuff,"mesh2")==0) fscanf(infp,"%s",mesh2);
    if(strcmp(namebuff,"compareU")==0) fscanf(infp,"%d",&compareU);
    if(strcmp(namebuff,"comparePsi")==0) fscanf(infp,"%d",&comparePsi);
  }
  cout<<" mu "<<mu<<endl;
  cout<<" eta "<<eta<<endl;
  cout<<" dt "<<dt<<endl;
  cout<<" model "<<model<<endl;
  cout<<" mesh1 "<<mesh1<<endl;
  cout<<" mesh2 "<<mesh2<<endl;
  cout<<" compareU "<<compareU<<endl;
  cout<<" comparePsi "<<comparePsi<<endl;

  m3dc1_scorec_init();
  apf::Mesh2* meshes[2];
  m3dc1_model_load(model);
  meshes[0] = apf::loadMdsMesh(m3dc1_model::instance()->model, mesh1);  
  meshes[1] = apf::loadMdsMesh(m3dc1_model::instance()->model, mesh2);  
  printMeshInfo(meshes[0]);
  apf::Field* fields[2];
  apf::MeshTag* tags[2];
  tags[0]=meshes[0]->findTag(fieldName);
  tags[1]=meshes[1]->findTag(fieldName);
  for(int i=0; i<2; i++)
  {
    fields[i]=createPackedField(meshes[i], "psi", 12);
    apf::MeshIterator* it = meshes[i]->begin(0);
    it = meshes[i]->begin(0);
    while (apf::MeshEntity * e = meshes[i]->iterate(it))
    {
      double dofBuff[1024];
      meshes[i]->getDoubleTag(e, tags[i],dofBuff);
      apf::setComponents(fields[i], e,0, dofBuff);
    }
    meshes[i]->end(it);
  }
  apf::MeshEntity* e;
  const int dim=2;
  MeshIterator* it = meshes[0]->begin(dim);
  Vector3 vec;
  double coord[3][2];
  Downward down_ent;
  apf::Field * errorField = createStepField(meshes[0], "error", apf::VECTOR);
  apf::Field * mapField = createPackedField(meshes[1], "solutionMap", 1);

  //apf::Field * errorFieldPsi = createStepField(meshes[0], "errorPsi", apf::SCALAR);

  double maxJ=0;
  double xyMaxJ[2];
  double gradMaxJ[2];
  double errorH1[2]={0.,0.}, errorH2[2]={0.,0.};
  while ((e = meshes[0]->iterate(it)))
  {
    int num_down_ent =  meshes[0]->getDownward(e, 0, down_ent);
    assert(num_down_ent==3);
    for (int i=0; i<num_down_ent; i++)
    {
      meshes[0]->getPoint(down_ent[i], 0, vec);
      coord[i][0]=vec[0];
      coord[i][1]=vec[1];
    }
    ReducedQuinticImplicit shape;
    shape.setCoord(coord);
    double psiDofs[18], UDofs[18];
    for(int i=0; i<num_down_ent; i++)
    {
      apf::getComponents(fields[0],down_ent[i],0,dofsBuff);
      for(int j=0; j<6; j++)
      {
        UDofs[i*6+j] = dofsBuff[j];
        psiDofs[i*6+j] = dofsBuff[6+j];
      }
    }
    vector<double> UOrg(numInter), psiOrg(numInter);
    vector<double> UOrgDX(numInter), psiOrgDX(numInter);
    vector<double> UOrgDY(numInter), psiOrgDY(numInter);
    vector<double> UOrgDXX(numInter), UOrgDXY(numInter), UOrgDYY(numInter);
    vector<double> psiOrgDXX(numInter), psiOrgDXY(numInter), psiOrgDYY(numInter);

    double evalRes[6];
    for(int i=0; i<numInter; i++)
    {
      double xy[2];
      double areaCoord[2]={beta[i], alpha[i]};
      area2xy(coord, areaCoord, xy);
      shape.setDofs(UDofs);
      shape.eval_g(xy, evalRes);
      UOrg[i]=evalRes[0];
      UOrgDX[i]=evalRes[1];
      UOrgDY[i]=evalRes[2];
      UOrgDXX[i]=evalRes[3];
      UOrgDXY[i]=evalRes[4];
      UOrgDYY[i]=evalRes[5];
      shape.setDofs(psiDofs);
      shape.eval_g(xy, evalRes);
      double J = evalRes[3]+evalRes[5];
      if(J>maxJ)
      {
        maxJ=J;
        xyMaxJ[0]=xy[0];
        xyMaxJ[1]=xy[1];
        gradMaxJ[0]=evalRes[1];
        gradMaxJ[1]=evalRes[2];
      }
      psiOrg[i]=evalRes[0];
      psiOrgDX[i]=evalRes[1];
      psiOrgDY[i]=evalRes[2];
      psiOrgDXX[i]=evalRes[3];
      psiOrgDXY[i]=evalRes[4];
      psiOrgDYY[i]=evalRes[5];
    }
    if(compareU==1 || comparePsi==1)
    {
      double error[3];
      error[0]=error[1]=error[2]=0.;
      for(int i=0; i<numInter; i++) 
      {
        double xy[2];
        double areaCoord[2]={beta[i], alpha[i]};
        area2xy(coord, areaCoord, xy);
        double U[6], psi[6];
        evalField (meshes[1], fields[1], xy, U, psi);
        if(compareU)
        {
          double errorU =0; // (UOrg[i]-U[0])*(UOrg[i]-U[0]);
          errorU += (UOrgDX[i]-U[1])*(UOrgDX[i]-U[1])/dt;    
          errorU += (UOrgDY[i]-U[2])*(UOrgDY[i]-U[2])/dt;
          errorH1[0] +=errorU;
          errorU += mu*(UOrgDXX[i]-U[3])*(UOrgDXX[i]-U[3]);
          //errorU += mu*2*(UOrgDXY[i]-U[4])*(UOrgDXY[i]-U[4]);
          errorU += 2*mu*(UOrgDXX[i]-U[3])*(UOrgDYY[i]-U[5]);
          errorU += mu*(UOrgDYY[i]-U[5])*(UOrgDYY[i]-U[5]);
          errorH2[0]+=errorU;
          error[0]+=errorU*area_weight[i];
        }
        if(comparePsi)
        {
          double errorPsi = 0; // (psiOrg[i]-psi[0])*(psiOrg[i]-psi[0]);
          errorPsi += (psiOrgDX[i]-psi[1])*(psiOrgDX[i]-psi[1])/dt;
          errorPsi += (psiOrgDY[i]-psi[2])*(psiOrgDY[i]-psi[2])/dt;
          errorH1[1]+=errorPsi;
          errorPsi += eta*(psiOrgDXX[i]-psi[3])*(psiOrgDXX[i]-psi[3]);
          //errorPsi += 2*eta*(psiOrgDXY[i]-psi[4])*(psiOrgDXY[i]-psi[4]);
          errorPsi += 2*eta*(psiOrgDXX[i]-psi[3])*(psiOrgDYY[i]-psi[5]);
          errorPsi += eta*(psiOrgDYY[i]-psi[5])*(psiOrgDYY[i]-psi[5]);
          errorH2[1] +=errorPsi;
          error[1]+=errorPsi*area_weight[i];
        }
      }
      double jacobi = (coord[1][1]-coord[2][1])*(coord[0][0]-coord[2][0]) + (coord[2][0]-coord[1][0])*(coord[0][1]-coord[2][1]);
      for(int i=0; i<2; i++)
      {
        error[i] *=jacobi;
        error[i]= sqrt(error[i]);
      }
      apf::setComponents(errorField, e,0, error);
      //apf::setComponents(errorFieldPsi, e,0, error+1);
    }
  }
  meshes[0]->end(it);
  cout<<"errorH1 "<<errorH1[0]<<" "<<errorH1[1]<<" errorH2 "<<errorH2[0]<<" " <<errorH2[1]<<endl;
  //cout<<" max current "<<maxJ<<" at "<<xyMaxJ[0]<<","<<xyMaxJ[1]<<" grad "<<gradMaxJ[0]<<","<<gradMaxJ[1]<<endl;
  gradMaxJ[0]=0.221417;
  gradMaxJ[1]=-0.975179;
  xyMaxJ[0]=1.09025;
  xyMaxJ[1]=2.69775;
  double len = gradMaxJ[0]*gradMaxJ[0]+gradMaxJ[1]*gradMaxJ[1];
  len = sqrt(len);
  len /= plotLen;
  gradMaxJ[0] /=len;
  gradMaxJ[1] /=len;

  FILE* fp = fopen("j.txt","w");
  for(int i=0; i<numPlot+1; i++)
  {
    double t=double(i)/numPlot-0.5;
    double xy[]={xyMaxJ[0]+t*gradMaxJ[0],xyMaxJ[1]+t*gradMaxJ[1]};
    double U[6], psi[6];
    evalField (meshes[0], fields[0], xy, U, psi);
    fprintf(fp, "%f %f %f %f\n",xy[0], xy[1], psi[3]+psi[5], psi[0]);
  }
  fclose(fp);
  apf::writeVtkFiles("error",meshes[0]);

  int vdim=0; // start map the solution from mesh 0 to mesh 1
  it = meshes[1]->begin(vdim);
  double cd[3];
  while ((e = meshes[1]->iterate(it)))
  {
    meshes[1]->getPoint(e, 0, vec);
    cd[0]=vec[0];
    cd[1]=vec[1];
    cd[2]=vec[2];
    double U[6], psi[6];
    evalField (meshes[0], fields[0], cd, U, psi);
    apf::setComponents(mapField, e,0, U); 
  }
  meshes[1]->end(it);
  apf::writeVtkFiles("map",meshes[1]);
  //m3dc1_scorec_finalize();
  MPI_Finalize();
}
 
