/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef REDUCEDQUINTICIMPLICIT_H
#define REDUCEDQUINTICIMPLICIT_H
#include "ReducedQuintic.h"
// see steve jardin 2004
class ReducedQuinticImplicit : public ReducedQuintic
{
public:
  ReducedQuinticImplicit(){};
  ~ReducedQuinticImplicit(){};
  virtual void setCoord(double coord[3][2]);
  // dofs in glb coord
  virtual void setDofs(double dofs_p[18]);
  virtual void print();
  virtual void eval_l( double coord[2], double res[6] );
  virtual void eval_g( double coord[2], double res[6] );
  virtual void loc2glb( double coord[2]);
  virtual void glb2loc( double coord[2]);
  virtual void loc2glbDofs(double dofs_p[18]);
  virtual void glb2locDofs(double dofs_p[18]);
  virtual double getJacobi();
private:
  void setGeoPara(double a_p, double b_p, double c_p, double theta_p[2], double origin_p[2]);
  void rotateCoord(double coord[2], double sin_theta_p, double cos_theta_p);
  void rotateDof(double dofs_p[6], double sin_theta_p, double cos_theta_p);
  void calcuCoeffs();
  double coeffs[18][20];
  double dofs[18];
  double a, b, c,sin_theta,cos_theta,origin[2];
  int order[3]; // order of vtx
}; 

#endif
