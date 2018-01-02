/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "CMODGeoExpression.h"
#include <cmath>
using namespace M3DC1;
double CMODExpressionR :: eval(double x) const
{
  return a + b*(cos(x + c*sin(x)));
}
double CMODExpressionR :: evalFirstDeriv(double x) const
{
  return -b*sin(x + c*sin(x))*(1.+c*cos(x));
}
double CMODExpressionR :: evalSecondDeriv(double x) const
{
  return -b*cos(x+c*sin(x))*(1.+c*cos(x))//*(1.+c*cos(x))
         +b*sin(x + c*sin(x))*c*sin(x);
}
double CMODExpressionZ :: eval(double x) const
{
  return d + e*sin(x);
}
double CMODExpressionZ :: evalFirstDeriv(double x) const
{
  return e*cos(x);
}
double CMODExpressionZ :: evalSecondDeriv(double x) const
{
  return -e*sin(x);
}

