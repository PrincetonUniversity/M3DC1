/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef SLNTRANSFERUTIL_H
#define SLNTRANSFERUTIL_H
const double PI=3.141592653589793238;
const double  tol =1.e-10;
// para[4]={a,b,c,sin_theta, cos_theta} , order is the index of the a b c vertex
void evalGeoPara(double vtxcoord[3][2], double origin[2], double para[5],int order[3] );
// given two points, get a triangle with a=b=c
// output newv,a,origin, theta
void getTriangle(double v1[2], double v2[2], double newv[2], double & a, double origin[2], double theta[2]);
void printmatrix(double* A,int M, int N, char* filename);
// if the pt contained by the triangle, return 1; else 0
//http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
int containPoint(double vtxcoord[3][2], double pt[2]);
#endif
