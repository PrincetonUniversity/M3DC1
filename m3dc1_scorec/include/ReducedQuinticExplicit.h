/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef REDUCEDQUINTICEXPLICIT_H
#define REDUCEDQUINTICEXPLICIT_H
#include "ReducedQuintic.h"
// see Barnhill 1981
class ReducedQuinticExplicit : public ReducedQuintic
{
public:
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
  double dofs[18];
}; 

#endif
