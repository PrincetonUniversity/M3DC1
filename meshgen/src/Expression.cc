/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "Expression.h"
#include <iostream>
#include <cmath>
using namespace M3DC1;
void M3DC1::dummyAnalyticExpression(double phi, double dummy, double *xyz, void* userdata){}
void M3DC1::evalCoord(double para, double *xyz, void* userdata)
{
  Expression* R_p=((Expression**)userdata)[0];
  Expression* Z_p=((Expression**)userdata)[1];
  xyz[0]=R_p->eval(para);
  xyz[1]=Z_p->eval(para);
  //xyz[2]=0.0;
}

void M3DC1::evalNormalVector(Expression* xp,Expression* yp, double para, double* normalvec)
{
  double dx=xp->evalFirstDeriv(para);
  double dy=yp->evalFirstDeriv(para);
  double arclen=sqrt(dx*dx+dy*dy);
  //std::cout<<" dx "<<dx<<" dy "<<dy<<std::endl;
  normalvec[0]=dy/arclen;
  normalvec[1]=-1.0*dx/arclen;
  //normalvec[2]=0.0;
}
void M3DC1::evalCurvature(Expression* xp,Expression* yp, double para, double* curv)
{
  double dx=xp->evalFirstDeriv(para);
  double ddx=xp->evalSecondDeriv(para);
  double dy=yp->evalFirstDeriv(para);
  double ddy=yp->evalSecondDeriv(para);
  //std::cout<<"dx dy "<<dx<<" "<<dy<<" ddx ddy "<<ddx<<" "<<ddy<<std::endl;
  *curv=(dx*ddy-dy*ddx)/((dx*dx+dy*dy)*sqrt(dx*dx+dy*dy));
}
int M3DC1::calcuBinomial( int j, int i)
{
  if(j<i) return 0;
  int res=1;
  for( int ii=0; ii<i; ii++)
    res*=j-ii;
  for( int ii=0; ii<i; ii++)
    res/=ii+1;
  return res;
}
