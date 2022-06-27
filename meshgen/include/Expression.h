/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef EXPRESSION_H
#define EXPRESSION_H
#include <vector>
using std::vector;
/**
Base class of 2D parametric curve
given parametric coordinate x, it evaluates
a. physical coordinate
b. first deriviative and normal vector
c. second deriviative and curvatur
*/
namespace M3DC1{
class Expression
{
public:
  virtual double eval(double x) const = 0;
  virtual double evalFirstDeriv(double x) const = 0;
  virtual double evalSecondDeriv(double x) const = 0;
};

void dummyAnalyticExpression(double phi, double dummy, double *xyz, void* userdata);
void evalCoord(double para, double *xyz, void* userdata);
void evalNormalVector(Expression* xp,Expression* yp, double para, double* normalvec) ;
void evalCurvature(Expression* xp,Expression* yp, double para, double* curv) ;
int calcuBinomial( int j, int i);
};
#endif
