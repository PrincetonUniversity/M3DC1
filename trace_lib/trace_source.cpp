#include "trace_source.h"

#include <iostream>
#include <math.h>

#define MAX_SURFACE_PTS 1024

bool trace_field_source::get_surface(const double r0, const double phi0, 
				     const double z0, double ds, 
				     double** r, double** z, int* n)
{
  double br0 = 0, bphi0 = 0, bz0 = 0;
  double br1 = 0, bphi1 = 0, bz1 = 0;
  double b;
  double rt[MAX_SURFACE_PTS], zt[MAX_SURFACE_PTS];
  int cycle = 0;
  bool closed = true;

  *n = 1;
  rt[0] = r0;
  zt[0] = z0;

  while(true) {
    // evaluate field at old point
    br0 = 0.;
    bphi0 = 0.;
    bz0 = 0.;
    if(!eval(rt[*n-1], phi0, zt[*n-1], &br0, &bphi0, &bz0)) {
      if(!closed) break;
      else {
	ds = -ds;
	rt[*n-1] = r0;
	zt[*n-1] = z0;
	closed = false;
	continue;
      }
    }

    // guess position of new point
    b = sqrt(br0*br0 + bz0*bz0);
    rt[*n] = rt[*n-1] + ds*br0/b;
    zt[*n] = zt[*n-1] + ds*bz0/b;

    //calculate field at new point
    br1 = 0.;
    bphi1 = 0.;
    bz1 = 0.;
    if(!eval(rt[*n], phi0, zt[*n], &br1, &bphi1, &bz1)) {
      if(!closed) break;
      else {
	ds = -ds;
	rt[*n-1] = r0;
	zt[*n-1] = z0;
	closed = false;
	continue;
      }
    }
  
    // use average of field at old point and guess point
    // to calculate new point
    br0 = (br0 + br1)/2.;
    bz0 = (bz0 + bz1)/2.;
    b = sqrt(br0*br0 + bz0*bz0);
  
    rt[*n] = rt[*n-1] + ds*br0/b;
    zt[*n] = zt[*n-1] + ds*bz0/b;

    // determine if we have completed the surface
    if(zt[*n-1] > z0 && zt[*n] <= z0) cycle++;
    else if(zt[*n-1] < z0 && zt[*n] >= z0) cycle++;
    if(cycle==2) {
      closed = true;
      break;
    }

    (*n)++;

    if(*n >= MAX_SURFACE_PTS) {
      std::cerr << "Warining: exceeded maximum number of surface points" 
		<< std::endl;
      closed = false;
      break;
    }
  }

  if(*n==0) {
    std::cerr << "Error: no points found on surface" << std::endl;
    return false;
  }

  // copy surface points to arrays
  *r = new double[*n];
  *z = new double[*n];
  for(int i=0; i<*n; i++) {
    (*r)[i] = rt[i];
    (*z)[i] = zt[i];
  }

  return closed;
}
