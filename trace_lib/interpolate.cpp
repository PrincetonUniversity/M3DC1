//=====================================================
// bicubic_interpolation_coeffs
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// calculates bicubic polynomial coefficients a of
// x, an array of dimension m x n,
// about index (i,j)
//=====================================================
bool bicubic_interpolation_coeffs(const double** x, const int m, const int n,
				  const int i, const int j, double a[4][4])
{
  if(i < 1 || i >= m-2) return false;
  if(j < 1 || j >= n-2) return false;

  a[0][0] = x[i][j];
  a[0][1] = (-2.*x[i][j-1]-3.*x[i][j]+6.*x[i][j+1]-x[i][j+2])/6.;
  a[0][2] = (    x[i][j-1]-2.*x[i][j]+   x[i][j+1]          )/2.;
  a[0][3] = (   -x[i][j-1]+3.*x[i][j]-3.*x[i][j+1]+x[i][j+2])/6.;
  
  a[1][0] = (-2.*x[i-1][j]- 3.*x[i][j]+6.*x[i+1][j]-x[i+2][j])/6.;
  a[1][1] = (4.*x[i-1][j-1] +6.*x[i-1][j]-12.*x[i-1][j+1]+2.*x[i-1][j+2]
	     +6.*x[i  ][j-1] +9.*x[i  ][j]-18.*x[i  ][j+1]+3.*x[i  ][j+2] 
	     -12.*x[i+1][j-1]-18.*x[i+1][j]+36.*x[i+1][j+1]-6.*x[i+1][j+2] 
	     +2.*x[i+2][j-1] +3.*x[i+2][j] -6.*x[i+2][j+1]   +x[i+2][j+2])/36.;
  a[1][2] = (-2.*x[i-1][j-1] +4.*x[i-1][j] -2.*x[i-1][j+1] 
	     -3.*x[i  ][j-1] +6.*x[i  ][j] -3.*x[i  ][j+1] 
	     +6.*x[i+1][j-1]-12.*x[i+1][j] +6.*x[i+1][j+1] 
	     -x[i+2][j-1]+ 2.*x[i+2][j]    -x[i+2][j+1])/12.;
  a[1][3] = (2.*x[i-1][j-1] -6.*x[i-1][j] +6.*x[i-1][j+1]-2.*x[i-1][j+2] 
	     +3.*x[i  ][j-1] -9.*x[i  ][j] +9.*x[i  ][j+1]-3.*x[i  ][j+2] 
	     -6.*x[i+1][j-1]+18.*x[i+1][j]-18.*x[i+1][j+1]+6.*x[i+1][j+2] 
	     +x[i+2][j-1] -3.*x[i+2][j] +3.*x[i+2][j+1]   -x[i+2][j+2])/36.;
  
  a[2][0] = (x[i-1][j]-2.*x[i][j]+x[i+1][j])/2.;
  a[2][1] = (-2.*x[i-1][j-1]-3.*x[i-1][j] +6.*x[i-1][j+1]   -x[i-1][j+2] 
	     +4.*x[i  ][j-1]+6.*x[i  ][j]-12.*x[i  ][j+1]+2.*x[i  ][j+2] 
	     -2.*x[i+1][j-1]-3.*x[i+1][j] +6.*x[i+1][j+1]   -x[i+1][j+2])/12.;
  a[2][2] = (   x[i-1][j-1]-2.*x[i-1][j]   +x[i-1][j+1] 
		-2.*x[i  ][j-1]+4.*x[i  ][j]-2.*x[i  ][j+1] 
		+x[i+1][j-1]-2.*x[i+1][j]   +x[i+1][j+1])/4.;
  a[2][3] = (  -x[i-1][j-1]+3.*x[i-1][j]-3.*x[i-1][j+1]+   x[i-1][j+2] 
	       +2.*x[i  ][j-1]-6.*x[i  ][j]+6.*x[i  ][j+1]-2.*x[i  ][j+2] 
	       -x[i+1][j-1]+3.*x[i+1][j]-3.*x[i+1][j+1]+   x[i+1][j+2])/12.;

  a[3][0] = (-x[i-1][j]+3.*x[i][j]-3.*x[i+1][j]+x[i+2][j])/6.;
  a[3][1] = (2.*x[i-1][j-1]+3.*x[i-1][j] -6.*x[i-1][j+1]   +x[i-1][j+2] 
	     -6.*x[i  ][j-1]-9.*x[i  ][j]+18.*x[i  ][j+1]-3.*x[i  ][j+2] 
	     +6.*x[i+1][j-1]+9.*x[i+1][j]-18.*x[i+1][j+1]+3.*x[i+1][j+2] 
	     -2.*x[i+2][j-1]-3.*x[i+2][j] +6.*x[i+2][j+1]   -x[i+2][j+2])/36.;
  a[3][2] = (  -x[i-1][j-1]+2.*x[i-1][j]   -x[i-1][j+1] 
	       +3.*x[i  ][j-1]-6.*x[i  ][j]+3.*x[i  ][j+1] 
	       -3.*x[i+1][j-1]+6.*x[i+1][j]-3.*x[i+1][j+1] 
	       +x[i+2][j-1]-2.*x[i+2][j]   +x[i+2][j+1])/12.;
  a[3][3] = (   x[i-1][j-1]-3.*x[i-1][j]+3.*x[i-1][j+1]   -x[i-1][j+2] 
		-3.*x[i  ][j-1]+9.*x[i  ][j]-9.*x[i  ][j+1]+3.*x[i  ][j+2] 
		+3.*x[i+1][j-1]-9.*x[i+1][j]+9.*x[i+1][j+1]-3.*x[i+1][j+2] 
		-x[i+2][j-1]+3.*x[i+2][j]-3.*x[i+2][j+1]   +x[i+2][j+2])/36.;
  return true;
}

//=====================================================
// cubic_interpolation_coeffs
// ~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// calculates cubic polynomial coefficients a of
// x, an array of dimension m,
// about index i
//=====================================================
bool cubic_interpolation_coeffs(const double* x, const int m, const int i,
				double* a)
{
  int ihermite;
  double xpi, xpip;

  ihermite = 0;
  switch(ihermite)
    {
    case 0:
    a[0] = x[i];
    if(i==0) {
      a[1] = (-3.*x[i] + 4.*x[i+1] - x[i+2])/2.;
      a[2] = (    x[i] - 2.*x[i+1] + x[i+2])/2.;
      a[3] = 0.;
    } else if(i == m-2) {
      a[1] = (-2.*x[i-1] - 3.*x[i] + 5.*x[i+1])/6.;
      a[2] = (    x[i-1] - 2.*x[i]    + x[i+1])/2.;
      a[3] = (   -x[i-1] + 3.*x[i] - 2.*x[i+1])/6.;
    } else if(i == m-1) {
      a[1] = (-x[i-1] + x[i])/3.;
      a[2] = ( x[i-1] - x[i])/2.;
      a[3] = (-x[i-1] + x[i])/6.;
      /*
        a[1] = (-x[i-1] + x[i])/4.;
        a[2] = ( x[i-1] - x[i])/2.;
        a[3] = (-x[i-1] + x[i])/4.;
      */
    } else {
      a[1] = (-2.*x[i-1] - 3.*x[i] + 6.*x[i+1] - x[i+2])/6.;
      a[2] = (    x[i-1] - 2.*x[i]    + x[i+1]         )/2.;
      a[3] = (   -x[i-1] + 3.*x[i] - 3.*x[i+1] + x[i+2])/6.;
    }
    break;

    case 1:
    a[0] = x[i];
    if(i==0) {
      xpip= .5*(x[i+2]-x[i]);
      a[1] = -2.*x[i] +2.*x[i+1] - xpip;
      a[2] =     x[i] -   x[i+1] + xpip;
      a[3] = 0.;
    } else if(i == m-1) {
      xpi = .5*(x[i+1]-x[i-1]);
      xpip=    (x[i+1]-x[i]);
      a[1] = xpi;
      a[2] = -.5*x[i+1] + 2.0*x[i] - 2.5*x[i-1] +     x[i-2];
      a[3] =  .5*x[i+1] - 1.5*x[i] + 1.5*x[i-1] - 0.5*x[i-2];
    } else if(i == m) {
      xpi = .5*(x[i+1]-x[i-1]);
      xpip=     x[i] - x[i-1];
      a[1] = xpip;
      a[2] = -.5*x[i] + 2.0*x[i-1] - 2.5*x[i-2] +     x[i-3];
      a[3] =  .5*x[i] - 1.5*x[i-1] + 1.5*x[i-2] - 0.5*x[i-3];
    } else {
      xpi = .5*(x[i+1]-x[i-1]);
      xpip= .5*(x[i+2]-x[i]);
      a[1] = xpi;
      a[2] = -3*x[i] +3*x[i+1] - 2*xpi - xpip;
      a[3] =  2*x[i] -2*x[i+1] +   xpi + xpip;
    }
    }
  return true;
}

//=====================================================
// cubic_interpolation
// ~~~~~~~~~~~~~~~~~~~
//
// interpolates function f(p) at point p0
// where f and p are arrays of length m
//=====================================================
bool cubic_interpolation(const int m, const double* p, const double p0, 
			 const double* f, double* f0)
{
  double a[4];
  int i, iguess;
  double dp;
  bool reverse;

  if(p[0] < p[m-1]) {
    reverse = false;
    if(p0 < p[0]) {
      *f0 = f[0];
      return false;
    }
    if(p0 > p[m-1]) {
      *f0 = f[m-1];
      return false;
    }
  } else {
    reverse = true;
    if(p0 > p[0]) {
      *f0 = f[0];
      return false;
    }
    if(p0 < p[m-1]) {
      *f0 = f[m-1];
      return false;
    }
  }

  iguess = m*(p0-p[0])/(p[m-1]-p[0]);

  if(p[iguess] < p0) {
    if(reverse) {
      // search backward
      for(i=iguess; i>=0; i--) {
	if(p[i] > p0) break;
      }
    } else {
      // search forward
      for(i=iguess; i<m-1; i++) {
	if(p[i+1] > p0) break;
      }
    }
  } else if(p[iguess] > p0) {
    if(reverse) {
      // search forward
      for(i=iguess; i<m-1; i++) {
	if(p[i+1] < p0) break;
      }
    } else { 
      // search backward
      for(i=iguess; i>=0; i--) {
	if(p[i] < p0) break;
      }
    }
  } else {
    i = iguess;
  }

  if(i < 0) i = 0;
  if(i > m-1) i = m-1;

  bool result = cubic_interpolation_coeffs(f,m,i,a);
  if(!result) return false;

  if(i >= m-1) {
    *f0 = a[0];
  } else {
    dp = (p0-p[i])/(p[i+1]-p[i]);
    *f0 = a[0] + (a[1] + (a[2] + a[3]*dp)*dp)*dp;
  }
  return true;
}
