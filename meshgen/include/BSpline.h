/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef BSPLINE_H
#define BSPLINE_H
#include <vector>
using std::vector;
#include "Expression.h"
#include "PolyNomial.h"
// for clamped b-spline, the same as sim geo advanced
/** BSpline implementation
It also includes conversions between BSpline and monic polynomial
*/
namespace M3DC1 
{
class BSpline: public Expression
{
public:
  BSpline(int order_p, vector<double> &ctrlPts_p, vector<double> & knots_p, vector<double> & weight_p); 
  BSpline():order(-1){}
  ~BSpline () {};
  virtual double eval(double x) const;
  virtual double evalFirstDeriv(double x) const;
  virtual double evalSecondDeriv(double x) const;
  void print();
  void getpara(int & order_p, vector<double> & ctrlPts_p, vector<double>&  knots_p, vector<double> & weight_p);
  BSpline & operator = ( const PolyNomial& pn);
private:
  void calcuDerivCoeff();
  int order;
  vector <double> ctrlPts;
  vector <double> knots;
  vector <double> weight;

  vector <double> ctrlPts_1st;
  vector <double> ctrlPts_2nd;
};
};
#endif
