/* 
   Copyright (C) 2005
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#include "modeler.h"
#include "AOMD.h"
#include "AOMD_cint.h"
#include "DiscreteModel.h"
#include <iostream>

#include "MeshAdapt.h"
#include "AdaptUtil.h"
#include "visUtil.h"
#include "PWLinearSField.h"
#include "fromMeshTools.h"
#include <stdlib.h>
#include <stdio.h>
#include "MeshSize.h"
#include "templateUtil.h"


int sizefield(pMesh,pSField,void *);
int sizefield1(pMesh,pSField,void *);
void preprocess(pMesh mesh);

double center[3], fine_radius;
double fine_size, coarse_size, gradation;

double L, L0, x_size, y_size, z_size, c_size;
double x_size_0, y_size_0, z_size_0;


int main(int argc, char *argv[])
{
printf("starting...\n");
  if(argc != 2 && argc !=3) {
    
    std::cout << "The run arguments are: " << argv[0] << " <mesh-file> <sizefield-file>\n";
    return 1;
  }
  
  if(argc == 2) {
    // keep the old discrete model interface
    pGModel model = new meshModel::DiscreteModel(argv[1], 2, 45, 45, 0, "struct-dmg.sms", "struct.dmg");
  }
  else if(argc == 3) {
    // inteface to define the benchmark mesh
    FILE *fp = fopen(argv[2], "r");

    int version;
    fscanf(fp, "%d", &version);
    if(version == 1) {
      // isotropic refinement
      fscanf(fp,"%lf %lf %lf %lf", &center[0], &center[1], &center[2], &fine_radius);
      fscanf(fp, "%lf %lf %lf\n", &fine_size, &coarse_size, &gradation);
      fclose(fp);
      
      printf("The defined geometry for the isotropic refined mesh is:\n");
      printf("The center is: %f %f %f\n", center[0], center[1], center[2]);
      printf("The radius for the fine zone is: %f\n", fine_radius);
      printf("The sizes are: %f %f %f\n", fine_size, coarse_size, gradation); 
    }
    else if(version == 2) {
      fscanf(fp,"%lf %lf %lf", &center[0], &center[1], &center[2]);
      fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fine_radius, &L, &L0, &x_size, &y_size, &z_size, &x_size_0, &y_size_0, &z_size_0, &c_size);
      printf("The defined geometry for the isotropic refined mesh is:\n");
      printf("The center is: %f %f %f\n", center[0], center[1], center[2]);
      printf("The anisotropic size is: %f %f %f %f %f %f \n", fine_radius, L, x_size, y_size, z_size, c_size);
      fclose(fp);
    }
    // load the mesh
    pMesh mesh=MS_newMesh(0);
    M_load(mesh, argv[1]);
    cout<<"Finish loading"<<endl;
    preprocess(mesh);

    cout<<"Finish prreprose"<<endl;
    
    pSField field=new PWLsfield(mesh);
    meshAdapt *rdr = new meshAdapt(mesh,field,1,0);
    
    if(version == 1)
      rdr->run(200,1, sizefield);
    else
      rdr->run(200,1, sizefield1);
    
    delete rdr;

    pGModel model = new meshModel::DiscreteModel(mesh, 2, 45, 45, 0, "struct-dmg.sms", "struct.dmg");
    
    M_delete(mesh);
    MS_exit();

  }
  return 0;
}


void preprocess(pMesh mesh)
{
  pPList elist = PList_new();
  double e_length;
  void *ptr;
  pEdge e;

  EIter eiter = M_edgeIter(mesh);
  
  while(e = EIter_next(eiter)) {
    e_length  = sqrt(adaptUtil::E_lengthSq(e));
    if(e_length > fine_radius)
      PList_append(elist, e);
  }
  EIter_delete(eiter);


  while(PList_size(elist)) {
    pSField field=new PWLsfield(mesh);
    meshAdapt *rdr = new meshAdapt(mesh,0,1,1);

    ptr = 0;
    while(e = (pEdge)PList_next(elist, &ptr))
      rdr->setAdaptLevel(e, 1);
    PList_delete(elist);
 
    rdr->run(1, 0, 0);
    delete rdr;
    delete field;

    elist = PList_new();
    eiter = M_edgeIter(mesh);
  
    while(e = EIter_next(eiter)) {
      e_length  = sqrt(adaptUtil::E_lengthSq(e));
      if(e_length > fine_radius)
	PList_append(elist, e);
    }
    EIter_delete(eiter);
    
  }

  PList_delete(elist);

  return;
}


int sizefield(pMesh mesh, pSField field, void *)
{
  
  VIter viter = M_vertexIter(mesh);
  pVertex vert;
  double vxyz[3], norm[3], R;



  while(vert = VIter_next(viter)) {
    V_coord(vert, vxyz);
    diffVt(vxyz, center, norm);
    R = sqrt(dotProd(norm, norm));
    
    if(R < fine_radius)     
      ((PWLsfield *)field)->setSize((pEntity)vert, fine_size);
    else
      ((PWLsfield *)field)->setSize((pEntity)vert, coarse_size);
  }
  
  VIter_delete(viter);

  EIter eiter = M_edgeIter(mesh);
  pEdge edge;
  int adjusted = 1;
  while(adjusted) {
    adjusted = 0;
    while(edge = EIter_next(eiter)) {
      pVertex v[2];
      double h[2], e_length, LM, values[2], update_size;
      int index = 0;

      for(int j=0; j<2; j++) {
	v[j] = E_vertex(edge, j);
	pMSize pT = ((PWLsfield *)field)->getSize((pEntity)v[j]);
	h[j] = pT->size(0);
      }

      
      e_length  = sqrt(adaptUtil::E_lengthSq(edge));
      
      LM = (h[0] - h[1]) / (e_length * (log (h[0]/h[1]))); 

      values[0] = pow((h[0]/h[1]), LM);
      values[1] = pow((h[1]/h[0]), LM);
      
      if(values[0] < values[1]) {
	double tmp = values[0];
	values[0] = values[1];
	values[1] = tmp;
	index = 1;
      }

      if(values[0] > 2.0) {
	update_size = h[!index] * pow(2., 1./LM);
	((PWLsfield *)field)->setSize((pEntity)v[index], update_size);
	adjusted = 1;
      }
      
      

    }
    EIter_reset(eiter);
    
  }

  EIter_delete(eiter);
 
  return 1;


}


int sizefield1(pMesh mesh, pSField field, void *)
{
  
  VIter viter = M_vertexIter(mesh);
  pVertex vert;
  double xyz[3], norm, norm1[3], R;
  double tol=0.01;
  double h[3], dirs[3][3];
  

  while(vert = VIter_next(viter)) {
    V_coord(vert, norm1);
    diffVt(norm1, center, xyz);
    R = sqrt(dotProd(xyz, xyz));

   
    if(R < fine_radius){
      h[0] = x_size * fabs(1. - exp (-fabs(R-fine_radius)*L)) + c_size;
      h[1] = y_size * fabs(1. - exp (-fabs(R-fine_radius)*L)) + z_size;
      h[2] = h[1];
    }
    else{
      h[0] = x_size_0 * fabs(1. - exp (-fabs(R-fine_radius)*L0)) + c_size;
      h[1] = y_size_0 * fabs(1. - exp (-fabs(R-fine_radius)*L0)) + z_size;
      h[2] = h[1];
    }

    norm=sqrt(R);
    if( norm>tol )
      {
        dirs[0][0]=xyz[0]/norm;
        dirs[0][1]=xyz[1]/norm;
        dirs[0][2]=xyz[2]/norm;
        if( xyz[0]*xyz[0] + xyz[1]*xyz[1] > tol*tol ) {
          dirs[1][0]=-1.0*xyz[1]/norm;
          dirs[1][1]=xyz[0]/norm;
          dirs[1][2]=0;
        } else {
          dirs[1][0]=-1.0*xyz[2]/norm;
          dirs[1][1]=0;
          dirs[1][2]=xyz[0]/norm;
        }
        crossProd(dirs[0],dirs[1],dirs[2]);
      }
    else
      {
        dirs[0][0]=1.0;
        dirs[0][1]=0.0;
        dirs[0][2]=0;
        dirs[1][0]=0.0;
        dirs[1][1]=1.0;
        dirs[1][2]=0;
        dirs[2][0]=0;
        dirs[2][1]=0;
        dirs[2][2]=1.0;
      }

    

    ((PWLsfield *)field)->setSize((pEntity)vert,dirs,h);
  }
  double beta[]={2.5,2.5,2.5};
  ((PWLsfield *)field)->anisoSmooth(beta);

  VIter_delete(viter);

  return 1;

}
