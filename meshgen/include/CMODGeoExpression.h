/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef CMODGEOEXPRESSION_H
#define CMODGEOEXPRESSION_H
#include <vector>
using std::vector;
#include "Expression.h"
/** Implementation of curves defined by triangluer functions*/
namespace M3DC1
{
class CMODExpressionR: public Expression
{
public:
  explicit CMODExpressionR(double a_p, double b_p, double c_p): a(a_p), b(b_p), c(c_p){}; 
  ~CMODExpressionR () {};
  virtual double eval(double x) const;
  virtual double evalFirstDeriv(double x) const;
  virtual double evalSecondDeriv(double x) const;
private:
  double a, b, c;
};
class CMODExpressionZ: public Expression
{
public:
  explicit CMODExpressionZ(double d_p, double e_p): d(d_p), e(e_p){};
  ~CMODExpressionZ () {};
  virtual double eval(double x) const;
  virtual double evalFirstDeriv(double x) const;
  virtual double evalSecondDeriv(double x) const;
private:
  double d, e;
};
};
#endif
