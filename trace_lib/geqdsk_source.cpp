#include "geqdsk_source.h"
#include "interpolate.h"

#include <iostream>
#include <fstream>
#include <iomanip>

/*
extern "C" {
  void bicubic_interpolation_coeffs_(double*, int*, int*, int*, int*, double*);
  void cubic_interpolation_(int*, double*, double*, double*, double*);
}
*/

static double pow(double x, int p)
{
  double r = 1.;
  if(p < 0) return 0.;
  for(int i=0; i<p; i++) r = r*x;
  return r;
}

geqdsk_source::geqdsk_source()
  : psi(0), psirz(0), fpol(0)
{
}

geqdsk_source::~geqdsk_source()
{
  if(psirz) {
    for(int i=0; i<nh; i++) delete[] psirz[i];
    delete[] psirz;
  }
  if(fpol) delete[] fpol;
  if(psi) delete[] psi;
}

bool geqdsk_source::load()
{
  std::fstream gfile;
  std::string dum;


  // Read eqdsk file
  gfile.open(filename.c_str(), std::fstream::in);

  if(!gfile) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return false;
  }

  double rdim, zdim, rcentr;
  double simag, sibry, bcentr, current;
    
  for(int i=0; i<5; i++) gfile >> dum;
  gfile >> nw >> nh;
  gfile >> rdim >> zdim >> rcentr >> rleft >> zmid;
  gfile >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
  gfile >> current >> simag >> dum >> rmaxis >> dum;
  gfile >> zmaxis >> dum >> sibry >> dum >> dum;

  std::cout << "nw = " << nw << ", nh = " << nh << std::endl;
  std::cout << "rmaxis = " << rmaxis << ", zmaxis = " << zmaxis << std::endl;
  std::cout << "rleft = " << rleft << ", zmid = " << zmid << std::endl;
  std::cout << "rdim = " << rdim << ", zdim = " << zdim << std::endl;
  std::cout << "simag = " << simag << std::endl;

  psirz = new double*[nh];
  for(int i=0; i<nh; i++) psirz[i] = new double[nw];
  fpol = new double[nw];  

  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> fpol[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // press
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // ffprime
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // pprime
  for(int i=0; i<nh; i++) 
    for(int j=0; j<nw; j++) gfile >> std::setw(16) >> psirz[i][j];

  gfile.close();

  // calculate spacing of grid points
  dx = rdim/(double)(nw - 1);
  dz = zdim/(double)(nh - 1);

  // set up array of flux values
  psi = new double[nw];
  for(int i=0; i<nw; i++) 
    psi[i] = (sibry - simag)*(double)i/(double)(nw - 1) + simag;

  return true;
}


bool geqdsk_source::eval(const double r, const double phi, const double z,
			 double* b_r, double* b_phi, double* b_z)
{
  double a[4][4];

  double p = (r-rleft)/dx;
  double q = (z-zmid)/dz + (double)nh/2.;
  int i = (int)p;
  int j = (int)q;

  /*
  // convert i, j to fortran indices
  i++; j++; p++; q++;
  bicubic_interpolation_coeffs_(psirz,&nw,&nh,&i,&j,a);
  */

  bool result=bicubic_interpolation_coeffs((const double**)psirz,nh,nw,j,i,a);
  if(!result) return false;

  double temp;
  double si = 0.;
  double dw = p-(double)i;
  double dh = q-(double)j;
  //  for(int n=0; n<4; n++) {
  //    for(int m=0; m<4; m++) {
  for(int n=0; n<4; n++) {
    for(int m=0; m<4; m++) {
      si += a[m][n]*pow(dw,n)*pow(dh,m);
      
      temp = a[m][n]*n*pow(dw,n-1)*pow(dh,m);
      *b_z -= temp/(r*dx);

      temp = a[m][n]*m*pow(dw,n)*pow(dh,m-1);
      *b_r += temp/(r*dz);
    }
  }

  double f;
  //  cubic_interpolation_(&nw, psi, &si, fpol, &f);
  cubic_interpolation(nw, psi, si, fpol, &f);

  *b_phi += f/r;
  
  return true;
}

bool geqdsk_source::center(double *r0, double *z0) const
{
  *r0 = rmaxis;
  *z0 = zmaxis;
  return true;
}

bool geqdsk_source::extent(double* r0, double* r1, double* z0, double* z1)
  const
{
  *r0 = rleft;
  *r1 = rleft + dx*(double)(nw-1);
  *z0 = zmid - dz*(double)((nh-1)/2);
  *z1 = zmid + dz*(double)((nh-1)/2);

  return true;
}
