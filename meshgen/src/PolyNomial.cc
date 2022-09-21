/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "PolyNomial.h"
#include <assert.h>
#include <iostream>
#include <cmath>
using std::cout;
using std::endl;
using namespace M3DC1;
PolyNomial::PolyNomial(int degree_p, vector <double> &coffs_p)
{
  degree=degree_p;
  coffs=coffs_p;
  assert(degree>0);
  assert(coffs.size()==degree);
  for( int i=0;i<degree-1;i++)
  {
    firstDerivCoffs.push_back(coffs[i]*(degree-1-i));
  }
  for( int i=0;i<degree-2;i++)
  {
    secondDerivCoffs.push_back(firstDerivCoffs[i]*(degree-2-i));
  }
}

void PolyNomial:: getcoeffs(vector <double> &coffs_p) const
{
  coffs_p=coffs;
}
void PolyNomial::print()
{
  cout<<" print PolyNomial: degree "<<degree<<endl;
  cout<<" expression: ";
  for(int i=0;i<degree;i++)
  { 
    if(i>0) cout<<" + ";
    cout<<coffs[i]<<"x^"<<degree-1-i;
  }
  cout<<endl;
  cout<<" first Deriv: ";
  if(degree-1<=0) cout<<"0";
  for(int i=0;i<degree-1;i++)
  {
    if(i>0) cout<<" + ";
    cout<<firstDerivCoffs[i]<<"x^"<<degree-2-i;
  }
  cout<<endl;
  cout<<" second Deriv coffs: ";
  if(degree-2<=0) cout<<"0";
  for(int i=0;i<degree-2;i++)
  {
    if(i>0) cout<<" + ";
    cout<<secondDerivCoffs[i]<<"x^"<<degree-3-i;
  }
  cout<<endl;
}
double PolyNomial::eval(double x) const
{
  double result=0;
  vector<double> x_n(degree,1.0);
  for( int i=degree-2;i>=0;i--)
  {
    x_n[i]=x_n[i+1]*x;
  }

  for( int i=0;i<degree;i++)
  {
     result+=coffs[i]*x_n[i];
  }
  
  return result;
}

double PolyNomial:: evalFirstDeriv(double x) const
{
  double result=0;
  if(degree==1)
    result=0.0;
  else
  {
    int thedegree=degree-1;
    vector<double> x_n(thedegree,1.0);
    for( int i=thedegree-2;i>=0;i--)
    {
      x_n[i]=x_n[i+1]*x;
    }

    for( int i=0;i<thedegree;i++)
    {
       result+=firstDerivCoffs[i]*x_n[i];
    }
   }
   return result;
}
  
double PolyNomial:: evalSecondDeriv(double x) const
{
  double result=0;
  if(degree<=2)
    result=0.0;
  else
  {
    int thedegree=degree-2;
    vector<double> x_n(thedegree,1.0);
    for( int i=thedegree-2;i>=0;i--)
    {
      x_n[i]=x_n[i+1]*x;
    }

    for( int i=0;i<thedegree;i++)
    {
       result+=secondDerivCoffs[i]*x_n[i];
    }
   }
   return result;
}

