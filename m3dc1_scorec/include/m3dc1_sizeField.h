/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_SIZEFIELD_H
#define M3DC1_SIZEFIELD_H
#include <apf.h>
#include <ma.h>
#include <apfVector.h>
#include <apfNumbering.h>
#include <iostream>
#include <vector>
#include <assert.h>

using namespace std;
class SizeFieldError : public ma::IsotropicFunction
{
  public:
    SizeFieldError(ma::Mesh* m, apf::Field* f, double targetError_p): mesh(m), field(f), targetError(targetError_p) {}
    virtual double getValue( ma::Entity* v)
    {
      double value;
      getComponents(field, v, 0, &value);
      assert(value>0);
      return value;
    }
    double getSize(ma::Entity* v)
    {
      apf:: Up edges;
      mesh->getUp(v,edges);
      double size=0;
      for(int i=0; i<edges.n; i++)
        size+=len(edges.e[i]);
      size/=edges.n;
      return size;
    }
    double len (ma::Entity* e)
    {
      apf::MeshEntity*  vertices[2];
      mesh->getDownward(e, 0, vertices);
      apf::Vector3 xyz2[2];
      for( int i=0; i<2; i++)
      {
        mesh->getPoint(vertices[i], 0, xyz2[i]);
      }
      return sqrt((xyz2[0][0]-xyz2[1][0])*(xyz2[0][0]-xyz2[1][0])+(xyz2[0][1]-xyz2[1][1])*(xyz2[0][1]-xyz2[1][1]));
    }
  private:
    ma::Mesh* mesh;
    apf::Field* field;
    double targetError;
};
class SizeFieldPsi : public ma::AnisotropicFunction
{
  public:
    SizeFieldPsi(apf::Field* f, double psi0_p, double psil_p, double param_p[13], int complexType_p)
    {
      field=f;
      complexType = complexType_p;
      psi0=psi0_p;
      psil=psil_p;
      for(int i=0; i<13; i++)
        param[i]=param_p[i];
    }
    virtual void getValue(
        ma::Entity* v,
        ma::Matrix& R,
        ma::Vector& h);
  private:
    apf::Field* field;
    int complexType;
    double psi0, psil;
    double param[13];
};

class Vortex : public ma::AnisotropicFunction
{
  public:
    Vortex(ma::Mesh* m, double center[3], double len)
    {
      std::cout<<" Vortex AnisotropicFunction center len "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" "<<len<<std::endl;
      mesh = m;
      modelLen = len;
      average = ma::getAverageEdgeLength(m);
      ma::Vector lower,upper;
      for(int i=0; i<3; i++)
        centroid[i]=center[i];
    }
    virtual void getValue(
        ma::Entity* v,
        ma::Matrix& R,
        ma::Vector& h)
    {
      ma::Vector x = ma::getPosition(mesh,v);
      double dx=x[0]-centroid[0];
      double dy=x[1]-centroid[1];
      double r=sqrt(dx*dx+dy*dy);
      if(r>1e-6) // if the vertex near to origin
      {
        dx=dx/r;
        dy=dy/r;
      }
      else
      {
        dx=1.;
        dy=0;
      }
      double small=0.2*average;
      double large=average;
      h[0]=small+large*modelLen*fabs(r-modelLen/3.);
      h[1]=large+large*fabs(r-modelLen/3.);
      h[2]=1.;
      R[0][0]=dx;
      R[1][0]=dy;
      R[2][0]=0;
      R[0][1]=-1.*dy;
      R[1][1]=dx;
      R[2][1]=0;   
      R[0][2]=0;
      R[1][2]=0;
      R[2][2]=1.;
    }
  private:
    ma::Mesh* mesh;
    double average;
    double modelLen;
    ma::Vector centroid;
};
#endif
