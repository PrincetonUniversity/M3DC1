/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_sizeField.h"
#include "m3dc1_mesh.h"
#include "PCU.h"
void SizeFieldPsi :: getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& h)
{
  //cout<<" getValue "<<endl;
  double norm;
  double hbar[3];
  double fldval;
  double psibar;
  double toltmp=0.01;
  assert(apf::getDimension(m3dc1_mesh::instance()->get_mesh(),v)==0);
  apf::Vector3 vcd;
  m3dc1_mesh::instance()->get_mesh()->getPoint(v, 0, vcd);
  //std::cout<<PCU_Comm_Self()<<" "<<v<<" getValue "<<vcd[0]<<" "<<vcd[1]<<" "<<vcd[2]<<std::endl;

  double value[12];
  assert(apf::countComponents(field)==(1+complexType)*6);
  ma::Vector xi(0,0,0);
  if(!hasEntity(field,v)) // happen in collaps???
  {
    std::cout<<PCU_Comm_Self()<<" warning to ma: collaps?? "<<std::endl;
    R[0][0]=1.;
    R[1][0]=0.;
    R[2][0]=0.0;
    R[0][1]=0;
    R[1][1]=1;
    R[2][1]=0.0;

    R[0][2]=0;
    R[1][2]=0;
    R[2][2]=1.;
    return;
  }
  getComponents(field, v, 0, value);
  fldval = value[0];
  psibar = (fldval - psi0)/(psil - psi0);
  if(psibar < param[0])
  {
    hbar[0] = param[3]*(1.-exp(-pow(fabs(psibar/param[0] - 1.), param[1]))) + param[8];
    hbar[1] = param[5]*(1.-exp(-pow(fabs(psibar/param[0] - 1.), param[1]))) + param[7];
    hbar[2] = hbar[1];
  }
  else
  {
    hbar[0] = param[4]*(1.-exp(-pow(fabs(psibar/param[0] - 1.), param[2]))) + param[8];
    hbar[1] = param[6]*(1.-exp(-pow(fabs(psibar/param[0] - 1.), param[2]))) + param[7];
    hbar[2] = hbar[1];
  }
  h[0] = 1./((1./hbar[0]) + (1./param[9])*(1./(1.+pow((psibar - param[12])/param[11], 2))));
  h[1] = 1./((1./hbar[1]) + (1./param[10])*(1./(1.+pow((psibar - param[12])/param[11], 2))));
  h[2] = h[1];

  double dpsidr=value[1];
  double dpsidz=value[2];

  double normgrad=sqrt(dpsidr*dpsidr+dpsidz*dpsidz);
  norm=sqrt(normgrad);
  if(norm>toltmp)
  {
    dpsidr=dpsidr/normgrad;
    dpsidz=dpsidz/normgrad;
  }
  else
  {
    dpsidr=1.0;
    dpsidz=0.0;
  }

  // use d(psi)/dr d(psi)/dz  as normal dirction
  R[0][0]=dpsidr;
  R[1][0]=dpsidz;
  R[2][0]=0.0;
  // use -d(psi)/dz d(psi)/dr as tangent dirction
  R[0][1]=-1.0*dpsidz;
  R[1][1]=dpsidr;
  R[2][1]=0.0;

  R[0][2]=0;
  R[1][2]=0;
  R[2][2]=1.;
  //double beta[3];
  //for (int i=0; i<3; i++)
    //beta[i]=1.0/smoothfact;
  //double beta[]={1.5,1.5,1.5};
  //((PWLsfield *)tensorField)->anisoSmooth(beta);
}
