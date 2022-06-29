/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "ReducedQuinticImplicit.h"
#include <slntransferUtil.h>
#include <iostream>
#include <assert.h>
#include <math.h>
using namespace std;

extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                 double* b, int* ldb, int* info );

void ReducedQuinticImplicit::setCoord(double coord[3][2])
{
    double origin[2];
    double para[5];
    // set up reduced quintic shape fns in the face
    evalGeoPara(coord, origin, para, order );
    setGeoPara(para[0],para[1],para[2],para+3,origin);
}
double ReducedQuinticImplicit::getJacobi()
{
  throw 1;
}
void ReducedQuinticImplicit::setGeoPara(double a_p, double b_p, double c_p, double theta_p[2], double origin_p[2])
{
  a=a_p;
  b=b_p;
  c=c_p;
  sin_theta=theta_p[0];
  cos_theta=theta_p[1];
  origin[0]=origin_p[0];
  origin[1]=origin_p[1];
  calcuCoeffs();
}

void ReducedQuinticImplicit::setDofs(double dofs_p[18])
{
  for( int i=0; i<3; i++)
  {
    int inode=order[i];
    int idx=inode*6;
    for(int j=0; j<6; j++)
      dofs[i*6+j]=dofs_p[idx+j];
  }
  glb2locDofs(dofs); 
}

// use lapack dgesv to solve the coffes matrix
void ReducedQuinticImplicit::calcuCoeffs()
{
  double a_sq=a*a;
  double a_cubic=a_sq*a;
  double a_quar=a_sq*a_sq;
  double a_quintic=a_sq*a_cubic;
  double b_sq=b*b;
  double b_cubic=b_sq*b;
  double b_quar=b_sq*b_sq;
  double b_quintic=b_sq*b_cubic;
  double c_sq=c*c;
  double c_cubic=c*c_sq;
  double c_quar=c_sq*c_sq;
  double c_quintic=c_sq*c_cubic;
  double A[20][20]=
    { 1., -b, 0., b_sq, 0, 0, -b_cubic, 0, 0, 0, b_quar, 0, 0, 0, 0, -b_quintic, 0, 0, 0, 0,
      0., 1., 0., -2.*b, 0., 0., 3*b_sq, 0, 0, 0, -4*b_cubic, 0, 0, 0, 0, 5*b_quar, 0, 0, 0, 0,
      0, 0, 1, 0, -b, 0, 0, b_sq, 0, 0, 0, -b_cubic, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, -6*b, 0, 0, 0, 12*b_sq, 0, 0, 0, 0, -20*b_cubic, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, -2*b, 0, 0, 0, 3*b_sq, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 0, -2*b, 0, 0, 0, 2*b_sq, 0, 0, 0, -2*b_cubic, 0, 0, 0,
      1, a, 0, a_sq, 0, 0, a_cubic, 0, 0, 0, a_quar, 0, 0, 0, 0, a_quintic, 0, 0, 0, 0,
      0, 1, 0, 2*a, 0, 0, 3*a_sq, 0, 0, 0, 4*a_cubic, 0, 0, 0, 0, 5*a_quar, 0, 0, 0, 0,
      0, 0, 1, 0, a, 0, 0, a_sq, 0, 0, 0, a_cubic, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 6*a, 0, 0, 0, 12*a_sq, 0, 0, 0, 0, 20*a_cubic, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 2*a, 0, 0, 0, 3*a_sq, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 0, 2*a, 0, 0, 0, 2*a_sq, 0, 0, 0, 2*a_cubic, 0, 0, 0,
      1, 0, c, 0, 0, c_sq, 0, 0, 0, c_cubic, 0, 0, 0, 0, c_quar, 0, 0, 0, 0, c_quintic,
      0, 1, 0, 0, c, 0, 0, 0, c_sq, 0, 0, 0, 0, c_cubic, 0, 0, 0, 0, c_quar, 0,
      0, 0, 1, 0, 0, 2*c, 0, 0, 0, 3*c_sq, 0, 0, 0, 0, 4*c_cubic, 0, 0, 0, 0, 5*c_quar,
      0, 0, 0, 2, 0, 0, 0, 2*c, 0, 0, 0, 0, 2*c_sq, 0, 0, 0, 0, 2*c_cubic, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 2*c, 0, 0, 0, 0, 3*c_sq, 0, 0, 0, 0, 4*c_cubic, 0,
      0, 0, 0, 0, 0, 2, 0, 0, 0, 6*c, 0, 0, 0, 0, 12*c_sq, 0, 0, 0, 0, 20*c_cubic,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*a_quar*c, 3*a_sq*c_cubic-2*a_quar*c, -2*a*c_quar+3*a_cubic*c_sq, c_quintic-4*a_sq*c_cubic, 5*a*c_quar,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*b_quar*c, 3*b_sq*c_cubic-2*b_quar*c, 2*b*c_quar-3*b_cubic*c_sq, c_quintic-4*b_sq*c_cubic, -5*b*c_quar }; 

     //printmatrix(&(A[0][0]),20,20,"A");
    // transpose A
    for ( int i=0; i<20; i++)
    {
      for ( int j=0;j<i; j++)
      {
        double a_ij=A[i][j];
        A[i][j]=A[j][i];
        A[j][i]=a_ij;
      }
    }
    // init the transpose(rhs)
    for( int i=0; i<18; i++)
    {
      for ( int j=0; j<20; j++)
        coeffs[i][j]=0.0;
      coeffs[i][i]=1.0;
    }
   
    //printmatrix(&(coeffs[0][0]), 18, 20, "coeffs"); 
    //solve to get coeffs
    int ipiv[20];
    int info;
    int N=20;
    int NRHS=18;
    dgesv_(&N,&NRHS,&(A[0][0]),&N, ipiv, &(coeffs[0][0]),&N,&info);
    if(info)
    {
      cout<<" error in LAPACK_dgesv "<<endl;
      throw 1;
    }
}

void ReducedQuinticImplicit::print()
{
}
void ReducedQuinticImplicit::eval_g( double coord[2], double res[6] )
{
  double coord_loc[2]={coord[0],coord[1]};
  glb2loc(coord_loc);
  eval_l( coord_loc, res );
}
void ReducedQuinticImplicit::eval_l( double coord[2], double res[6] )
{
  double coord_loc[2]={coord[0],coord[1]};
  double beta[20],coord_poly[20];
  int poly_idx[20][2]={0,0,1,0,0,1,2,0,1,1,0,2,3,0,2,1,1,2,0,3,4,0,3,1,2,2,1,3,0,4,5,0,3,2,2,3,1,4,0,5};
  for( int i=0; i<20; i++)
  {
    beta[i]=0;
    for( int j=0; j<18; j++)
    {
      beta[i]+=coeffs[j][i]*dofs[j];
    }
    //cout<<"beta["<<i<<"] "<<beta[i]<<endl;
  }
  double xi_poly[6], eta_poly[6];
  for( int i=0; i<6; i++)
  {
    xi_poly[i]=1.0;
    eta_poly[i]=1.0;
    if(i>0)
    {
      xi_poly[i]=xi_poly[i-1]*coord_loc[0];
      eta_poly[i]*=eta_poly[i-1]*coord_loc[1];
    }
  }
  // interpolate the 6 dofs
  for( int i=0; i<6; i++)
    res[i]=0.0;
  for( int i=0; i<20; i++)
  { 
    int xi_idx=poly_idx[i][0];
    int eta_idx=poly_idx[i][1];
    int xi_coff=xi_idx;
    int eta_coff=eta_idx;
    int xixi_coff=xi_coff*(xi_coff-1);
    int xieta_coff=xi_coff*eta_coff;
    int etaeta_coff=eta_coff*(eta_coff-1);
    res[0]+=beta[i]*xi_poly[xi_idx]*eta_poly[eta_idx];
    if(xi_coff)
    {
      assert(xi_idx>0);
      res[1]+=beta[i]*xi_coff*xi_poly[xi_idx-1]*eta_poly[eta_idx];
      if(xixi_coff)
      {
        assert(xi_idx>1);
        res[3]+=beta[i]*xixi_coff*xi_poly[xi_idx-2]*eta_poly[eta_idx];
      }
      if(eta_coff)
      {
        assert(eta_idx>0);
        res[4]+=beta[i]*xieta_coff*xi_poly[xi_idx-1]*eta_poly[eta_idx-1];
      }
    }
    if(eta_coff)
    {
      assert(eta_idx>0);
      res[2]+=beta[i]*eta_coff*xi_poly[xi_idx]*eta_poly[eta_idx-1];
      if(etaeta_coff)
      {
        assert(eta_idx>1);
        res[5]+=beta[i]*etaeta_coff*xi_poly[xi_idx]*eta_poly[eta_idx-2];
      }
    }
  }
  rotateDof(res,-sin_theta,cos_theta);
}

void ReducedQuinticImplicit:: loc2glb( double coord[2])
{
  rotateCoord(coord,-sin_theta,cos_theta);
  coord[0]+=origin[0];
  coord[1]+=origin[1];
}

void ReducedQuinticImplicit:: glb2loc( double coord[2])
{
  coord[0]-=origin[0];
  coord[1]-=origin[1];
  rotateCoord(coord,sin_theta, cos_theta);
}
void ReducedQuinticImplicit:: loc2glbDofs(double dofs_p[18])
{
  for(int i=0; i<3; i++)
  {
    rotateDof(dofs_p+i*6,-sin_theta,cos_theta);
  }
}
void ReducedQuinticImplicit:: glb2locDofs(double dofs_p[18])
{
  for(int i=0; i<3; i++)
  {
    rotateDof(dofs_p+i*6,sin_theta,cos_theta);
  }
}
void ReducedQuinticImplicit:: rotateCoord(double coord[2], double sin_theta_p, double cos_theta_p)
{
  double buffer;
  buffer=cos_theta_p*coord[0]+sin_theta_p*coord[1];
  coord[1]=cos_theta_p*coord[1]-sin_theta_p*coord[0];
  coord[0]=buffer;
}
void ReducedQuinticImplicit:: rotateDof(double dofs_p[6], double sin_theta_p, double cos_theta_p)
{
  double dofs_r[6];
  dofs_r[0]=dofs_p[0];
  dofs_r[1]=cos_theta_p*dofs_p[1]+sin_theta_p*dofs_p[2];
  dofs_r[2]=cos_theta_p*dofs_p[2]-sin_theta_p*dofs_p[1];
  dofs_r[3]=cos_theta_p*cos_theta_p*dofs_p[3]+2*sin_theta_p*cos_theta_p*dofs_p[4]+sin_theta_p*sin_theta_p*dofs_p[5];
  dofs_r[4]=-sin_theta_p*cos_theta_p*dofs_p[3]+(cos_theta_p*cos_theta_p-sin_theta_p*sin_theta_p)*dofs_p[4]+sin_theta_p*cos_theta_p*dofs_p[5];
  dofs_r[5]=sin_theta_p*sin_theta_p*dofs_p[3]-2*sin_theta_p*cos_theta_p*dofs_p[4]+cos_theta_p*cos_theta_p*dofs_p[5];
  for( int i=0; i<6; i++)
    dofs_p[i]=dofs_r[i];
}
