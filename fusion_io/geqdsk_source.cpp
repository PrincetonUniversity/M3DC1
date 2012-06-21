#include "fusion_io.h"
#include "interpolate.h"

#include <iostream>
#include <fstream>
#include <iomanip>

static double pow(double x, int p)
{
  double r = 1.;
  if(p < 0) return 0.;
  for(int i=0; i<p; i++) r = r*x;
  return r;
}


geqdsk_source::geqdsk_source()
  : psi(0), psirz(0), fpol(0), ffprime(0), press(0)
{
}

geqdsk_source::~geqdsk_source()
{
  close();
}

int geqdsk_source::open(const char* filename)
{
  std::fstream gfile;
  std::string dum;

  // Read eqdsk file
  gfile.open(filename, std::fstream::in);

  if(!gfile) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return 1;
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

  psirz = new double[nw*nh];
  ffprime = new double[nw];
  fpol = new double[nw];
  press = new double[nw];

  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> fpol[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> press[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> ffprime[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // pprime
  for(int i=0; i<nw*nh; i++) gfile >> std::setw(16) >> psirz[i];

  gfile.close();

  // calculate spacing of grid points
  dx = rdim/(double)(nw - 1);
  dz = zdim/(double)(nh - 1);

  // set up array of flux values
  psi = new double[nw];
  for(int i=0; i<nw; i++) 
    psi[i] = (sibry - simag)*(double)i/(double)(nw - 1) + simag;

  return 0;
}

int geqdsk_source::close()
{
  if(psirz) delete[] psirz;
  if(fpol) delete[] fpol;
  if(ffprime) delete[] ffprime;
  if(press) delete[] press;
  if(psi) delete[] psi;

  psirz = 0;
  fpol = 0;
  ffprime = 0;
  press = 0;
  psi = 0;

  return FIO_SUCCESS;
}

int geqdsk_source::center(double *r0, double *z0) const
{
  *r0 = rmaxis;
  *z0 = zmaxis;
  return FIO_SUCCESS;
}

int geqdsk_source::extent(double* r0, double* r1, double* z0, double* z1)
  const
{
  *r0 = rleft;
  *r1 = rleft + dx*(double)(nw-1);
  *z0 = zmid - dz*(double)((nh-1)/2);
  *z1 = zmid + dz*(double)((nh-1)/2);

  return FIO_SUCCESS;
}

int geqdsk_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  return FIO_SUCCESS;
}

int geqdsk_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_CURRENT_DENSITY);
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_TOTAL_PRESSURE);

  return FIO_SUCCESS;
}

int geqdsk_source::get_field(const field_type t,fio_field** f,
			     const fio_option_list* opt)
{
  *f = 0;
  fio_field* mf;

  switch(t) {
  case(FIO_CURRENT_DENSITY):
    mf = new geqdsk_current_density(this);
    break;

  case(FIO_MAGNETIC_FIELD):
    mf = new geqdsk_magnetic_field(this);
    break;

  case(FIO_TOTAL_PRESSURE):
    mf = new geqdsk_pressure_field(this);
    break;

  default:
    return FIO_UNSUPPORTED;
  };
 
  *f = mf;
  return FIO_SUCCESS;
}

int geqdsk_source::interpolate_psi(const double r0, const double z0,
				   double* si) const
{
  double a[16];

  double p = (r0-rleft)/dx;
  double q = (z0-zmid)/dz + (double)nh/2.;
  int i = (int)p;
  int j = (int)q;
  int ierr;

  if(i < 1 || i > nw) return FIO_OUT_OF_BOUNDS;
  if(j < 1 || j > nh) return FIO_OUT_OF_BOUNDS;

  // convert i, j to fortran indices
  i++; j++; p++; q++;
  bicubic_interpolation_coeffs_(psirz,&nw,&nh,&p,&q,a,&ierr);
  if(ierr!=0) {
    std::cerr << "Interpolation error" << std::endl;
    return ierr;
  }

  for(int n=0; n<6; n++)
    si[i] = 0.;

  double temp;
  for(int n=0; n<4; n++) {
    for(int m=0; m<4; m++) {
      int index = m*4 + n;

      si[0] += a[index]      *pow(p-i,n  )*pow(q-j,m  );
      
      temp = a[index]      *n*pow(p-i,n-1)*pow(q-j,m  );
      si[1] += temp/dx;

      temp = a[index]      *m*pow(p-i,n  )*pow(q-j,m-1);
      si[2] += temp/dz;

      temp = a[index]*n*(n-1)*pow(p-i,n-2)*pow(q-j,m  );
      si[3] += temp/(dx*dx);

      temp = a[index]    *m*n*pow(p-i,n-1)*pow(q-j,m-1);
      si[4] += temp/(dx*dz);

      temp = a[index]*m*(m-1)*pow(p-i,n  )*pow(q-j,m-2);
      si[5] += temp/(dz*dz);
    }
  }
  
  return FIO_SUCCESS;
}
