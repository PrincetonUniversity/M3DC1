#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int nsplines = 100;
const int nout = 500;

static void find_center(const int n, const double* x, const double* z, 
			double* x0, double* z0);
static int parse_file(const char* filename, double** x, double** z);
static void scale_dn(const int n, const double *x, const double *z, 
		     const double x0, const double z0, 
		     const int nspl, const double* coeffs, 
		     double *t1, double* r1, const double scalefac);
static int test_data(const double x0, const double z0,
		      const int ntheta, const double* r0,
		      double** x, double**z);
static void write_data(const char* filename, const int n, 
		       const double* t, const double* r,
		       const double x0, const double z0);
static void write_gnuplot(const char* filename,
			  const char* f0, const char* f1,
			  const double scalefac);

static void hermite_spline_terms(const double t, const int nspl, 
				 double* c, const int nder);
static void solve_spline_coeffs(const int ndata, 
				const double* x, const double* z,
				const double x0, const double z0,
				const int nspl, double* c);
static void get_spline_values(const int nvals, const double* t, double* r, 
			      const int nspl, const double* c, const int nder);

static int least_squares(const int m, const int n, 
			  double* a, double* b, double* x);
static void transpose(const int m, const int n, double*);

int main(int argc, char* argv[])
{
  int i;
  double xmag, zmag, scalefac=1.;
  double *x0, *z0, *x1, *z1, *r0, *t0, *t, *r;
  int n0, n1;
  double* coeffs;

  if(argc < 3) {
    return 1;
  }

  // read equilibrium surface
  n0 = parse_file(argv[1], &x0, &z0);
  if(n0 == 0) return 1;

  // find geometric center of equilibrium surface data
  find_center(n0, x0, z0, &xmag, &zmag);
  printf("Center found at (%g, %g)\n", (float)(xmag), (float)(zmag));

  if(argc >=4) scalefac = atof(argv[3]);
  printf("Scalefac = %g\n", (float)scalefac);

  // find approximate equilibrium surface as r(theta)
  coeffs = malloc(2*nsplines*sizeof(double));
  solve_spline_coeffs(n0, x0, z0, xmag, zmag, nsplines, coeffs);
  free(x0);
  free(z0);

  // output approximate equilibrium surface
  t0 = malloc(nout*sizeof(double));
  r0 = malloc(nout*sizeof(double));
  for(i=0; i<nout; i++) t0[i] = 2.*M_PI*(double)i/(double)nout;
  get_spline_values(nout, t0, r0, nsplines, coeffs, 0);
  printf("Writing...\n");
  write_data("data0.out", nout, t0, r0, xmag, zmag);
  free(t0);
  free(r0);

  // read perturbed surface
  n1 = parse_file(argv[2], &x1, &z1);
  //n1 = test_data(xmag, zmag, ntheta, r0, &x1, &z1);
  if(n1 == 0) {
    printf("Error\n");
    return 1;
  }

  // scale the difference between the perturbed and equilibrium surface
  printf("Scaling data.\n");
  t = calloc(n1, sizeof(double));
  r = calloc(n1, sizeof(double));
  scale_dn(n1, x1, z1, xmag, zmag, nsplines, coeffs, t, r, scalefac);

  // output scaled perturbed surface
  printf("Writing data.\n");
  write_data("data1.out", n1, t, r, xmag, zmag);

  // write gnuplot script for visualizing surfaces
  printf("Writing gnuplot script.\n");
  write_gnuplot("gp", "data0.out", "data1.out", scalefac);

  printf("Done.\n");

  // free allocated memory
  free(x1);
  free(z1);
  free(t);
  free(r);
  free(coeffs);

  return 0;
}


int parse_file(const char* filename, double** x, double** z)
{
  FILE* file;
  char* line = 0;
  size_t len = 0;
  ssize_t read;
  int n = 0;
  int i;
  float dx, dz;

  file = fopen(filename, "r");

  if(!file) {
    printf("Error: could not open %s\n", filename);
    return 0;
  }

  // count number of lines in data file
  while((read = getline(&line, &len, file)) != -1) n++;
  if(line) free(line);

  *x = malloc(sizeof(double)*n);
  *z = malloc(sizeof(double)*n);

  rewind(file);
  for(i=0; i<n; i++) {
    fscanf(file, "%f\t%f\n", &dx, &dz);
    (*x)[i] = dx;
    (*z)[i] = dz;
  }

  fclose(file);
  return n;
}

void find_center(const int n, const double* x, const double* z, 
		 double*x0, double* z0)
{
  int i;
  *x0 = 0;
  *z0 = 0;

  for(i=0; i<n; i++) {
    *x0 += x[i];
    *z0 += z[i];
  }
  *x0 /= (double)n;
  *z0 /= (double)n;
}

static void hermite_spline_terms(const double t, const int nspl, double* c,
				 const int nder)
{
  int i, j;
  double dt, s;
  const double delta_t = 2.*M_PI/(double)nspl;

  for(i=0; i<2*nspl; i++)
    c[i] = 0.;

  i = (int)(t/delta_t) % nspl;
  j = (i+1) % nspl;

  dt = t - delta_t*(double)i;
  s = dt/delta_t;

  switch(nder) {
  case(0):
    c[2*i  ] = (1.     - 3.*s*s + 2.*s*s*s);
    c[2*i+1] = (     s - 2.*s*s +    s*s*s);
    c[2*j  ] = (         3.*s*s - 2.*s*s*s);
    c[2*j+1] = (       -    s*s +    s*s*s);
    break;
  case(1):
    c[2*i  ] = (   - 6.*s + 6.*s*s)/delta_t;
    c[2*i+1] = (1. - 4.*s + 3.*s*s)/delta_t;
    c[2*j  ] = (     6.*s - 6.*s*s)/delta_t;
    c[2*j+1] = (   - 2.*s + 3.*s*s)/delta_t;
    break;
  case(2):
    c[2*i  ] = (-6. + 12.*s)/(delta_t*delta_t);
    c[2*i+1] = (-4. +  6.*s)/(delta_t*delta_t);
    c[2*j  ] = ( 6. - 12.*s)/(delta_t*delta_t);
    c[2*j+1] = (-2. +  6.*s)/(delta_t*delta_t);
    break;
  default:
    printf("derivatives of order %d are not supported.\n", nder);
  }
}

static void solve_spline_coeffs(const int ndata, 
				const double* x, const double* z,
				const double x0, const double z0,
				const int nspl, double* coeffs)
{
  double r, t, dx, dz;
  double* cmat, *rhs;
  int i;
  
  cmat= malloc(ndata*nspl*2*sizeof(double));
  rhs = malloc(ndata*sizeof(double));
  
  for(i=0; i<ndata; i++) {
    dx = x[i]-x0;
    dz = z[i]-z0;
    r = sqrt(dx*dx + dz*dz);
    t = atan2(dz,dx);
    if(t<0.) t+=2.*M_PI;
    rhs[i] = r;
    hermite_spline_terms(t, nspl, &(cmat[2*i*nspl]), 0);
  }

  least_squares(ndata,2*nspl,cmat,rhs,coeffs);

  free(rhs);
  free(cmat);
}

static void get_spline_values(const int nvals, const double* t, double* r, 
			      const int nspl, const double* coeffs,
			      const int nder)
{
  int i, j;
  double* terms;
  
  terms = malloc(2*nspl*sizeof(double));

  for(i=0; i<nvals; i++) {
    hermite_spline_terms(t[i], nspl, terms, nder);
    r[i] = 0.;
    for(j=0; j<2*nspl; j++) {
      r[i] += terms[j]*coeffs[j];
    }
  }

  free(terms);
}

void scale_dn(const int n, const double *x, const double *z, 
	      const double x0, const double z0, 
	      const int nspl, const double *coeffs,
	      double* t1, double* r1,
	      const double scalefac)
{
  int i, k;
  double dx, dz, t, xp, zp, r, drdt;
  double nx, nz, xn, zn, l;
  double dxdt, dzdt;
  double last_t;
  const double tol=1e-4;
  const int max_its = 100;

  for(i=0; i<n; i++) {
    // initial guess for theta
    t = atan2(z[i]-z0,x[i]-x0);
    if(t < 0.) t += 2.*M_PI;

    for(k=0; k<max_its; k++) {
      // evaluate r, r', and r'' at theta = t
      get_spline_values(1, &t, &r, nspl, coeffs, 0);
      get_spline_values(1, &t, &drdt, nspl, coeffs, 1);
            
      // this is the point on the surface with theta = t
      xp = x0 + r*cos(t);
      zp = z0 + r*sin(t);

      // find direction normal to surface at theta = t
      dxdt = drdt*cos(t) - r*sin(t);
      dzdt = drdt*sin(t) + r*cos(t);
      nx =  dzdt;
      nz = -dxdt;
      l = sqrt(nx*nx + nz*nz);
      nx /= l;
      nz /= l;

      // distance from data point to surface point
      dx = x[i] - xp;
      dz = z[i] - zp;

      // new guess for a point on the surface
      xn = x[i] - (nx*dx + nz*dz)*nx;
      zn = z[i] - (nx*dx + nz*dz)*nz;
      last_t = t;
      t = atan2(zn-z0, xn-x0);
      if(t < 0.) t += 2.*M_PI;

      if(fabs(last_t-t) <= tol) break;
    }
    if(k==max_its) {
      printf("Warning: did not converge after %d iterations\n", max_its);
      printf("%g %g\n", last_t, t);
    }

    xp += scalefac*dx;
    zp += scalefac*dz;

    //    printf("old, new (%g,%g), (%g,%g)\n", x[i], z[i], xp, zp);

    // move back to r,theta coordinates
    r1[i] = sqrt((xp-x0)*(xp-x0) + (zp-z0)*(zp-z0));
    t1[i] = atan2((zp-z0), (xp-x0));
    if(t1[i]<0.) t1[i] += 2.*M_PI;
  }
}

void write_data(const char* filename, const int n, const double* t, 
		const double *r, const double x0, const double z0)
{
  FILE* file;
  int i;
  double x, z;

  file = fopen(filename, "w");

  for(i=0; i<n; i++) {
    x = x0 + r[i]*cos(t[i]);
    z = z0 + r[i]*sin(t[i]);
    fprintf(file, "%g\t%g\n", (float)x, (float)z);
  }

  fclose(file);
}

void write_gnuplot(const char* filename, const char*f0, const char*f1,
		   const double scalefac)
{
  FILE* file;

  file = fopen(filename, "w");

  fprintf(file, "set size ratio -1\n");
  fprintf(file, "set xlabel 'R (m)'\n");
  fprintf(file, "set ylabel 'Z (m)'\n");
  fprintf(file, "plot '%s' w lines title 'Equilibrium',\\\n \
 \'%s' w dots title 'Perturbed (x%g)\n", f0, f1, scalefac);

  fclose(file);
}

void transpose(const int m, const int n, double* a)
{
  int i, j;
  double* d = malloc(m*n*sizeof(double));

  for(i=0; i<m; i++) {
    for(j=0; j<n; j++) {
      d[j*m + i] = a[i*n + j];
    }
  }

  memcpy(a, d, m*n*sizeof(double));
  free(d);
}

int least_squares(const int m, const int n, double* a, double* b, double* x)
{
  int nrhs, lda, ldb, lwork, info, i;
  double opt;
  char trans;
  double* work;

  // put data into column-major form for fortran call
  transpose(m,n,a);

  trans = 'N';
  nrhs = 1;
  lda = m;
  ldb = m;

  printf("solving...\n");
  // determine optimal workspace size
  lwork = -1;
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &opt, &lwork, &info);
  lwork = (int)opt;
  printf("optimal size of work = %d\n", lwork);
  work = malloc(lwork*sizeof(double));
  // solve least squares problem
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  free(work);
  printf("done solving.  Info=%d\n", info);

  // copy solution
  for(i=0; i<n; i++)
    x[i] = b[i];

  return info;
}

int test_data(const double x0, const double z0,
	       const int ntheta, const double* r0,
	       double** x, double**z)
{
  int i;
  double r,t;

  (*x) = malloc(ntheta*sizeof(double));
  (*z) = malloc(ntheta*sizeof(double));

  for(i=0; i<ntheta; i++) {
    t = 2.*M_PI*(double)i/(double)ntheta;
    r = r0[i] + 0.01*sin(5.*t);
    (*x)[i] = x0 + r*cos(t);
    (*z)[i] = z0 + r*sin(t);
  }

  return ntheta;
}
