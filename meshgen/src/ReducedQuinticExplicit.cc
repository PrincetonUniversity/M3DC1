/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "ReducedQuinticExplicit.h"
extern "C"
{
  void set_element_coord_i_ ( int*, double*);
  void eval_field_i_( double *, double*, double*, double*);
  void loc_xyz_2_glb_i_ ( double*, double*);
  void glb_xyz_2_loc_i_ ( double*, double*);
  void glb_dof_2_loc_i_ ( double*, double*);
  void loc_dof_2_glb_i_ ( double*, double*);
  void get_jacobi_i_( double *);
}
double ReducedQuinticExplicit:: getJacobi()
{
  double jac;
  get_jacobi_i_(&jac);
  return jac;
}
void ReducedQuinticExplicit:: setCoord(double coord[3][2])
{
  double coord_f[2][3];
  // convert to fortran convention
  for(int i=0; i<3; i++)
  {
    for( int j=0; j<2; j++)
      coord_f[j][i]=coord[i][j];
  }
  int dummy=0;
  set_element_coord_i_(&dummy,(double*)coord_f);
}
void ReducedQuinticExplicit:: setDofs(double dofs_p[18])
{
  for(int i=0;i<18; i++)
    dofs[i]=dofs_p[i];
  glb2locDofs(dofs);
}
void ReducedQuinticExplicit:: print()
{
}
void ReducedQuinticExplicit:: eval_g( double coord[2], double res[6] )
{
  double coord_loc[2]={coord[0],coord[1]};
  glb2loc(coord_loc);
  eval_field_i_(coord_loc,coord_loc+1,dofs,res);
}
void ReducedQuinticExplicit:: eval_l( double coord[2], double res[6] )
{
  double coord_loc[2]={coord[0],coord[1]};
  eval_field_i_(coord_loc,coord_loc+1,dofs,res);
}
void ReducedQuinticExplicit:: loc2glb( double coord[2])
{
  double coord_glb[2];
  loc_xyz_2_glb_i_(coord,coord_glb);
  coord[0]=coord_glb[0];
  coord[1]=coord_glb[1];
}
void ReducedQuinticExplicit:: glb2loc( double coord[2])
{
  double coord_loc[2];
  glb_xyz_2_loc_i_(coord,coord_loc);
  coord[0]=coord_loc[0];
  coord[1]=coord_loc[1];
}
void ReducedQuinticExplicit:: loc2glbDofs(double dofs_p[18])
{
  double dofs_glb[18];
  loc_dof_2_glb_i_(dofs_p,dofs_glb);
  for(int i=0; i<18; i++)
    dofs_p[i]=dofs_glb[i];
}
void ReducedQuinticExplicit:: glb2locDofs(double dofs_p[18])
{
  double dofs_loc[18];
  glb_dof_2_loc_i_(dofs_p,dofs_loc);
  for(int i=0; i<18; i++)
    dofs_p[i]=dofs_loc[i];
}
