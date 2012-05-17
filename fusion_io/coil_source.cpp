#include <iostream>
#include <math.h>

#include "coil_source.h"

void coil_segment::biot_savart_integrand(
		       const double r0, const double phi0, const double z0,
		       const double r1, const double phi1, const double z1,
		       const double dr, const double dphi, const double dz,
		       double *ir, double* iphi, double* iz)
{
  double delta_phi = phi0 - phi1;
  double delta_z = z0 - z1;
  double co = cos(delta_phi);
  double sn = sin(delta_phi);

  double dr2 = r0*r0 + r1*r1 - 2.*r0*r1*co + delta_z*delta_z;
  double dr32 = 1./pow(dr2,1.5);

  double delta_r[3];
  delta_r[0] = r0 - r1*co;
  delta_r[1] = r1*sn;
  delta_r[2] = delta_z;
  double dl[3];
  dl[0] =  dr*co + r1*dphi*sn;
  dl[1] = -dr*sn + r1*dphi*co;
  dl[2] =  dz;

  *ir   = (dl[1]*delta_r[2] - dl[2]*delta_r[1])*dr32;
  *iphi = (dl[2]*delta_r[0] - dl[0]*delta_r[2])*dr32;
  *iz   = (dl[0]*delta_r[1] - dl[1]*delta_r[0])*dr32;
}

bool coil_segment::eval(const double r, const double phi, const double z,
			double* br, double* bphi, double* bz)
{
  const int steps = 10;
  const double dr   = (R[1] - R[0])/(double)(steps);
  const double dphi = (Phi[1] - Phi[0])/(double)(steps);
  const double dz   = (Z[1] - Z[0])/(double)(steps);

  const double fac = 1.e-7*current;

  double r1 = R[0] + dr/2.;
  double phi1 = Phi[0] + dphi/2.;
  double z1 = Z[0] + dz/2.;

  for(int i=0; i<steps; i++) {
    double ir, iphi, iz;
    biot_savart_integrand(r, phi, z, r1, phi1, z1, dr, dphi, dz,
			  &ir, &iphi, &iz);

    *br += fac*ir;
    *bphi += fac*iphi;
    *bz += fac*iz;

    r1 += dr;
    phi1 += dphi;
    z1 += dz;
  }

  return true;
}

coil_source::coil_source()
{
  //interpolate = true;
}

coil_source::~coil_source()
{
}

int coil_source::load()
{
  return 0;
}

bool coil_source::eval(const double r, const double phi, const double z,
		       double* br, double* bphi, double* bz)
{
  iterator i = begin();

  while(i != end()) {
    if(!(i->eval(r,phi,z,br,bphi,bz))) return false;
    i++;
  };

  return true;
}
