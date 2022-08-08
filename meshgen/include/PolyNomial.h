/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
using std::vector;
#include "Expression.h"
/** Implementation of monomic polynomial*/
namespace M3DC1
{
class PolyNomial: public Expression
{
public:
  explicit PolyNomial(int degree_p, vector <double> &coffs_p); 
  ~PolyNomial () {};
  void getcoeffs(vector <double> &coffs_p) const;
  virtual double eval(double x) const;
  virtual double evalFirstDeriv(double x) const;
  virtual double evalSecondDeriv(double x) const;
  void print();
private:
  int degree;
  vector <double> coffs;
  vector <double> firstDerivCoffs;
  vector <double> secondDerivCoffs;
};

};
#endif
