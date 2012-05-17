#include <math.h>
#include <iostream>

#include "diiid_coils.h"

void diiid_icoils(coil_source* s, const double current, const int n, 
		  const double top_phase, const double bottom_phase)
{
  double c;
  const double r0 = 2.372;
  const double r1 = 2.161;
  const double z0 = 0.508;
  const double z1 = 1.005;
  const double dphi = 2.*3.14159265/6.;

  std::cerr << "Creating DIII-D I-coil set with\n"
	    << " current = " << current << " A\n"
	    << " n = " << n << "\n"
	    << " phase (top, bottom) = " 
	    << top_phase << ", " << bottom_phase
	    << std::endl;

  for(int i=0; i<6; i++) {
    double phi0 = dphi*(double)i;
    double phi1 = phi0 + dphi;

    // top coil
    c = current*cos(phi0*(double)n + 3.14159265*top_phase/180.);
    s->push_back(coil_segment(c, r0, r0, phi0, phi1, z0, z0));
    s->push_back(coil_segment(c, r0, r1, phi1, phi1, z0, z1));
    s->push_back(coil_segment(c, r1, r1, phi1, phi0, z1, z1));
    s->push_back(coil_segment(c, r1, r0, phi0, phi0, z1, z0));

    // bottom coil
    c = current*cos(phi0*(double)n + 3.14159265*bottom_phase/180.);
    s->push_back(coil_segment(c, r1, r1, phi0, phi1, -z1, -z1));
    s->push_back(coil_segment(c, r1, r0, phi1, phi1, -z1, -z0));
    s->push_back(coil_segment(c, r0, r0, phi1, phi0, -z0, -z0));
    s->push_back(coil_segment(c, r0, r1, phi0, phi0, -z0, -z1));
  }
}
