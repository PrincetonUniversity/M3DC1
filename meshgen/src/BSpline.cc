/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "BSpline.h"
#include <assert.h>
#include <iostream>
#include <cmath>
using std::cout;
using std::endl;
using namespace M3DC1;
BSpline :: BSpline(int order_p, vector<double> &ctrlPts_p, vector<double> & knots_p, vector<double> & weight_p)
{
  assert(order_p>1);
  order=order_p;
  ctrlPts=ctrlPts_p;
  knots=knots_p;
  //assert(knots_p.size()-ctrlPts_p.size()==order_p);
  weight=weight_p;
  calcuDerivCoeff();
}
void BSpline :: calcuDerivCoeff()
{
 // caculate coeffs of first deriv
  for( int i=1; i<ctrlPts.size(); i++)
  {
    double delta=double((order-1))/(knots.at(i+order-1)-knots.at(i));
    ctrlPts_1st.push_back((ctrlPts.at(i)-ctrlPts.at(i-1))*delta);
  }
  assert(ctrlPts_1st.size()==ctrlPts.size()-1);
  // caculate coeffs of second deriv
  for( int i=1; i<ctrlPts_1st.size(); i++)
  { 
    double delta=double((order-2))/(knots.at(i+order-1)-knots.at(i+1));
    ctrlPts_2nd.push_back((ctrlPts_1st.at(i)-ctrlPts_1st.at(i-1))*delta);
  }
  if(order>1) assert(ctrlPts_2nd.size()==ctrlPts.size()-2);
}
double BSpline :: eval(double x) const
{
  // first find the interval of x in knots
  int leftKnot=order-1;
  int leftPt=0;
  while(knots.at(leftKnot+1)<x)
  {
    leftKnot++;
    leftPt++;
    if(leftKnot==knots.size()-1) break;
  }
  vector<double> pts(&(ctrlPts[leftPt]),&ctrlPts[leftPt+order]);
  vector<double> localKnots(&(knots[leftKnot-order+2]),&(knots[leftKnot+order]));
  for( int r=1; r<=order; r++)
  {
    // from bottom to top to save a buff
    for( int i=order-1; i>=r; i-- )
    {
      double a_left=localKnots.at(i-1);
      double a_right=localKnots.at(i+order-r-1);
      double alpha;
      if(a_right==a_left) alpha=0.; // not sure??
      else  alpha=(x-a_left)/(a_right-a_left);
      pts.at(i)=(1.-alpha)*pts.at(i-1)+alpha*pts.at(i); 
    }
  }
  return pts.at(order-1);
}
double BSpline :: evalFirstDeriv(double x) const
{
  // first find the interval of x in knots
  int leftKnot=order-1;
  int leftPt=0;
  while(knots.at(leftKnot+1)<x)
  {
    leftKnot++;
    leftPt++;
  }
  int order_t=order-1;
  vector<double> pts(&(ctrlPts_1st.at(leftPt)),&(ctrlPts_1st[leftPt+order_t]));
  vector<double> localKnots(&(knots.at(leftKnot-order_t+2)),&(knots[leftKnot+order_t]));
  for( int r=1; r<=order_t; r++)
  {
    // from bottom to top to save a buff
    for( int i=order_t-1; i>=r; i-- )
    {
      double a_left=localKnots.at(i-1);
      double a_right=localKnots.at(i+order_t-r-1);
      double alpha;
      if(a_right==a_left) alpha=0.; // not sure??
      else  alpha=(x-a_left)/(a_right-a_left);
      pts.at(i)=(1.-alpha)*pts.at(i-1)+alpha*pts.at(i);
    }
  }
  return pts.at(order_t-1);
}
double BSpline :: evalSecondDeriv(double x) const
{
  if(order==2) return 0;
 // first find the interval of x in knots
  int leftKnot=order-1;
  int leftPt=0;
  while(knots.at(leftKnot+1)<x)
  {
    leftKnot++;
    leftPt++;
  }
  int order_t=order-2;
  vector<double> pts(&(ctrlPts_2nd.at(leftPt)),&(ctrlPts_2nd[leftPt+order_t]));
  vector<double> localKnots(&(knots.at(leftKnot-order_t+2)),&(knots[leftKnot+order_t]));
  for( int r=1; r<=order_t; r++)
  {
    // from bottom to top to save a buff
    for( int i=order_t-1; i>=r; i-- )
    {
      double a_left=localKnots.at(i-1);
      double a_right=localKnots.at(i+order_t-r-1);
      double alpha;
      if(a_right==a_left) alpha=0.; // not sure??
      else  alpha=(x-a_left)/(a_right-a_left);
      pts.at(i)=(1.-alpha)*pts.at(i-1)+alpha*pts.at(i);
    }
  }
  return pts.at(order_t-1);
}
void BSpline :: print()
{
  cout<<" ctrlPts "<<ctrlPts.size()<<endl;
  for( int i=0; i<ctrlPts.size(); i++)
    cout<< ctrlPts.at(i)<<" ";
  cout<<endl;
  cout<<" knots "<<knots.size()<<endl;
  for( int i=0; i<knots.size(); i++)
    cout<< knots.at(i)<<" ";
  cout<<endl;
}
void BSpline :: getpara(int& order_p, vector<double>& ctrlPts_p, vector<double>&  knots_p, vector<double>&  weight_p)
{
  order_p=order;
  ctrlPts_p=ctrlPts;
  knots_p=knots;
  weight_p=weight;
}
//H.Prautzsch Springer,2002
BSpline & BSpline :: operator = ( const PolyNomial& pn)
{
  vector <double> coffs_p;
  pn.getcoeffs(coffs_p);
  order=coffs_p.size();
  ctrlPts.clear();
  knots.clear();
  weight.clear();
  ctrlPts_1st.clear();
  ctrlPts_2nd.clear();
  knots.resize(order+order);
  ctrlPts.resize( order);
  for(int i=0; i<order; i++)
  {
    knots.at(i)=0.;
    knots.at(order+i)=1.0;
    double ctr=0.0;
    for( int j=0; j<=i; j++)
    {
      ctr+=(double)calcuBinomial(order-1-j,i-j)*coffs_p.at(order-1-j);
    }
    ctr/=(double)calcuBinomial(order-1,i);
    ctrlPts.at(i)=ctr;
  }
  double tol=1e-6;
  assert(fabs(ctrlPts.at(0)-pn.eval(0.0))<tol);
  assert(fabs(ctrlPts.at(order-1)-pn.eval(1.0))<tol);
  calcuDerivCoeff();
  return *this;
}

