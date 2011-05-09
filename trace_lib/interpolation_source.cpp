#include <iostream>

#include "interpolation_source.h"

interpolation_source::~interpolation_source()
{
}

void interpolation_source::free_grid()
{
  if(Br) delete[] Br;
  if(Bphi) delete[] Bphi;
  if(Bz) delete[] Bz;
  Br = Bphi = Bz = 0;
}

void interpolation_source::set_extent(const double r0, const double r1,
				      const double z0, const double z1,
				      const int nr, const int nphi, 
				      const int nz)
{
  r_extent[0] = r0;
  r_extent[1] = r1;
  z_extent[0] = z0;
  z_extent[1] = z1;
  n[0] = nr;
  n[1] = nphi;
  n[2] = nz;
}

bool interpolation_source::load()
{
  if(sources.size()==0)
    return true;

  if(r_extent[0] >= r_extent[1] || 
     z_extent[0] >= z_extent[1])  {
    std::cerr << "Error: interpolation grid has improper bounds" << std::endl;
    return false;
  }
  if(n[0] < 2 || n[1] < 1 || n[2] < 2) {
    std::cerr << "Error: interpolation grid has to few nodes" << std::endl;
    return false;
  }

  free_grid();

  Br = new double[n[0]*n[1]*n[2]];
  Bphi = new double[n[0]*n[1]*n[2]];
  Bz = new double[n[0]*n[1]*n[2]];

  if(!Br || !Bphi || !Bz) {
    std::cerr << "Error: interpolation grid allocation failed." << std::endl;
    return false;
  }

  double dr   = (r_extent[1] - r_extent[0])/(double)(n[0] - 1);
  double dphi = 2.*3.14159265/(double)n[1];
  double dz   = (z_extent[1] - z_extent[0])/(double)(n[2] - 1);
   
  std::cerr << "Calculating interpolation source..." << std::endl;

  int l=0;
  double r = r_extent[0];
  for(int i=0; i<n[0]; i++) {
    r += dr;
    double phi = 0.;
    for(int j=0; j<n[1]; j++) {
      phi += dphi;
      double z = z_extent[0];
      for(int k=0; k<n[2]; k++) {
	z += dz;
	Br[l] = 0.;
	Bphi[l] = 0.;
	Bz[l] = 0.;
	trace_source_list::iterator p = sources.begin();
	while(p != sources.end()) {
	  if(!(*p)->eval(r, phi, z, &Br[l], &Bphi[l], &Bz[l]))
	    return false;
	  p++;
	}
	l++;
      }
    }
  }

  return true;
}

bool interpolation_source::eval
(const double r, const double phi, const double z,
 double *br, double *bphi, double* bz)
{
  double x[3];
  x[0] = (r - r_extent[0])/(r_extent[1] - r_extent[0]) * (double)n[0];
  x[1] = phi/(2.*3.14159265) * (double)n[1];
  x[2] = (z - z_extent[0])/(z_extent[1] - z_extent[0]) * (double)n[2];
  double new_br, new_bphi, new_bz;
  // if(!trilinear_interpolation(x, n, &new_br, &new_bphi, &new_bz))
  if(!triquadratic_interpolation(x, n, &new_br, &new_bphi, &new_bz))
  //  if(!tricubic_interpolation(x, n, &new_br, &new_bphi, &new_bz))
    return false;

  *br += new_br;
  *bphi += new_bphi;
  *bz += new_bz;

  return true;
}

bool interpolation_source::center(double*, double*) const
{
  return false;
}

bool interpolation_source::extent(double*, double*, double*, double*) const
{
  return false;
}

bool interpolation_source::trilinear_interpolation
(const double x[3], const int n[3], double* br, double* bphi, double* bz)
  const
{
  int i[3];

  for(int j=0; j<3; j++)
    i[j] = (int)x[j];

  double dx[3];
  for(int j=0; j<3; j++) 
    dx[j] = x[j] - i[j];

  // check bounds
  if(i[0] < 0) return false;
  else if(i[0]+1 >= n[0]) return false;
  i[1] = (i[1] + n[1]) % n[1];
  if(i[2] < 0) return false;
  else if(i[2]+1 >= n[2]) return false;

  // calculate indices
  int index[2][2][2];
  int d[3];

  for(int j=0; j<2; j++) {
    d[0] = (i[0] + j)*n[1]*n[2];
    for(int k=0; k<2; k++) {
      if(i[1] + k >= n[1]) d[1] = (i[1] + k - n[1])*n[2];
      else                 d[1] = (i[1] + k       )*n[2];
      for(int l=0; l<2; l++) {
	d[2] = i[2] + l;
	index[j][k][l] = d[0] + d[1] + d[2];
      }
    }
  }

  // do interpolation
  double *f, *o;

  for(int j=0; j<3; j++) {
    if(j==0)      { f = Br;   o = br;   }
    else if(j==1) { f = Bphi; o = bphi; }
    else if(j==2) { f = Bz;   o = bz;   }

    *o = 
      (1.-dx[2]) * 
      ((1.-dx[1])*(f[index[0][0][0]]*(1.-dx[0]) + f[index[1][0][0]]*dx[0])
       +   dx[1] *(f[index[0][1][0]]*(1.-dx[0]) + f[index[1][1][0]]*dx[0]))
      + dx[2] *
      ((1.-dx[1])*(f[index[0][0][1]]*(1.-dx[0]) + f[index[1][0][1]]*dx[0])
       +   dx[1] *(f[index[0][1][1]]*(1.-dx[0]) + f[index[1][1][1]]*dx[0]));
  }

  return true;
}

bool interpolation_source::triquadratic_interpolation
(const double x[3], const int n[3], double* br, double* bphi, double* bz)
  const
{
  int i[3];

  for(int j=0; j<3; j++)
    i[j] = (int)(x[j] + 0.5);

  double dx[3];
  for(int j=0; j<3; j++) 
    dx[j] = x[j] - i[j];

  // check bounds
  if(i[0] < 1) return false;
  else if(i[0]+1 >= n[0]) return false;
  i[1] = (i[1] + n[1]) % n[1];
  if(i[2] < 1) return false;
  else if(i[2]+1 >= n[2]) return false;

  // determine indices 
  int index[3][3][3];
  int d[3];
  for(int j=0; j<3; j++) {
    d[0] = (i[0] + j - 1)*n[1]*n[2];
    for(int k=0; k<3; k++) {
      if(i[1] + k - 1 < 0)          d[1] = (i[1] + k - 1 + n[1])*n[2];
      else if(i[1] + k - 1 >= n[1]) d[1] = (i[1] + k - 1 - n[1])*n[2];
      else                          d[1] = (i[1] + k - 1       )*n[2];
      for(int l=0; l<3; l++) {
	d[2] = i[2] + l - 1;
	index[j][k][l] = d[0] + d[1] + d[2];
      }
    }
  }
  
  // do interpolation
  double *f, *o;
  double g[3], h[3];

  for(int j=0; j<3; j++) {
    if(j==0)      { f = Br;   o = br;   }
    else if(j==1) { f = Bphi; o = bphi; }
    else if(j==2) { f = Bz;   o = bz;   }

    for(int l=0; l<3; l++) {
      for(int m=0; m<3; m++) {
	g[m] = f[index[l][m][1]]
	  + dx[2]*((f[index[l][m][2]] - f[index[l][m][0]])/2.
		   + dx[2]*((f[index[l][m][2]] + f[index[l][m][0]])/2.
			    - f[index[l][m][1]]));
      }
      h[l] = g[1] 
	+ dx[1]*((g[2] - g[0])/2.
		 + dx[1]*((g[2] + g[0])/2. - g[1]));
    }


    *o = h[1] 
      + dx[0]*((h[2] - h[0])/2. 
	       + dx[0]*((h[2] + h[0])/2. - h[1]));
  }
  return true;
}

bool interpolation_source::tricubic_interpolation
(const double x[3], const int n[3], double* br, double* bphi, double* bz)
  const
{
  int i[3];

  for(int j=0; j<3; j++)
    i[j] = (int)x[j];

  double dx[3];
  for(int j=0; j<3; j++) 
    dx[j] = x[j] - i[j];

  // check bounds
  if(i[0] < 1) return false;
  else if(i[0]+2 >= n[0]) return false;
  i[1] = (i[1] + n[1]) % n[1];
  if(i[2] < 1) return false;
  else if(i[2]+2 >= n[2]) return false;

  // determine indices 
  int index[4][4][4];
  int d[3];
  for(int j=0; j<4; j++) {
    d[0] = (i[0] + j - 1)*n[1]*n[2];
    for(int k=0; k<4; k++) {
      if(i[1] + k - 1 < 0)          d[1] = (i[1] + k - 1 + n[1])*n[2];
      else if(i[1] + k - 1 >= n[1]) d[1] = (i[1] + k - 1 - n[1])*n[2];
      else                          d[1] = (i[1] + k - 1       )*n[2];
      for(int l=0; l<4; l++) {
	d[2] = i[2] + l - 1;
	index[j][k][l] = d[0] + d[1] + d[2];
      }
    }
  }

  // do interpolation
  double *s, *o;
  double f[4], g[4], h[4];

  for(int j=0; j<3; j++) {
    if(j==0)      { s = Br;   o = br;   }
    else if(j==1) { s = Bphi; o = bphi; }
    else if(j==2) { s = Bz;   o = bz;   }

    for(int l=0; l<4; l++) {
      for(int m=0; m<4; m++) {
	for(int n=0; n<4; n++) 
	  f[n] = s[index[l][m][n]];
	
	g[m] = f[1]
	  + dx[2]*((-2.*f[0] - 3.*f[1] + 6.*f[2] - f[3])/6.
		   + dx[2]*((f[0] - 2.*f[1] + f[2])/2.
			    + dx[2]*(-f[0] + 3.*f[1] - 3.*f[2] + f[3])/6.));
      }
      h[l] = g[1]
	+ dx[1]*((-2.*g[0] - 3.*g[1] + 6.*g[2] - g[3])/6.
		 + dx[1]*((g[0] - 2.*g[1] + g[2])/2.
			  + dx[1]*(-g[0] + 3.*g[1] - 3.*g[2] + g[3])/6.));
    }

    *o = h[1]
      + dx[0]*((-2.*h[0] - 3.*h[1] + 6.*h[2] - h[3])/6.
	       + dx[0]*((h[0] - 2.*h[1] + h[2])/2.
			+ dx[0]*(-h[0] + 3.*h[1] - 3.*h[2] + h[3])/6.));
  }
  return true;
}

