/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef REDUCEDQUINTIC_H
#define REDUCEDQUINTIC_H
// see Barnhill 1981
class ReducedQuintic
{
public:
  ReducedQuintic(){};
  virtual ~ReducedQuintic(){};
  virtual void setCoord(double coord[3][2])=0;
  // dofs in glb coord
  virtual void setDofs(double dofs_p[18])=0;
  virtual void print()=0;
  // coord in glb coord, res in glb coord
  virtual void eval_g( double coord[2], double res[6] )=0;
  // coord in loc coord, res in glb coord
  virtual void eval_l( double coord[2], double res[6] )=0;
  virtual void loc2glb( double coord[2])=0;
  virtual void glb2loc( double coord[2])=0;
  virtual void loc2glbDofs(double dofs_p[18])=0;
  virtual void glb2locDofs(double dofs_p[18])=0;
  virtual double getJacobi()=0;
}; 

#endif
