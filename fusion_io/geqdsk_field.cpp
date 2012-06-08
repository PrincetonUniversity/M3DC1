#include "fusion_io.h"

extern "C" {
  void bicubic_interpolation_(const int*, const int*, const double*, const double*, 
			      const double*, const double*, const double*, 
			      double*, double*, double*, int*);
  void bicubic_interpolation_coeffs_(const double*, const int*, const int*, const int*, const int*, 
				     double*, int*);
  void cubic_interpolation_(const int*, const double*, const double*, const double*, double*);
}

static double pow(double x, int p)
{
  double r = 1.;
  if(p < 0) return 0.;
  for(int i=0; i<p; i++) r = r*x;
  return r;
}

int geqdsk_magnetic_field::eval(const double* x, double* b)
{
  double a[16];

  double p = (x[0]-source->rleft)/source->dx;
  double q = (x[2]-source->zmid)/source->dz + (double)source->nh/2.;
  int i = (int)p;
  int j = (int)q;
  int ierr;

  if(i < 1 || i > source->nw) return false;
  if(j < 1 || j > source->nh) return false;

  // convert i, j to fortran indices
  i++; j++; p++; q++;
  bicubic_interpolation_coeffs_(source->psirz,&source->nw,&source->nh,&i,&j,a,&ierr);
  if(ierr!=0) return ierr;

  b[0] = b[2] = 0;

  double temp;
  double si = 0.;
  for(int n=0; n<4; n++) {
    for(int m=0; m<4; m++) {
      int index = m*4 + n;

      si += a[index]*pow(p-i,n)*pow(q-j,m);
      
      temp = a[index]*n*pow(p-i,n-1)*pow(q-j,m  );
      b[2] -= temp/(x[0]*source->dx);

      temp = a[index]*m*pow(p-i,n  )*pow(q-j,m-1);
      b[0] += temp/(x[0]*source->dz);
    }
  }

  double f;
  cubic_interpolation_(&source->nw, source->psi, &si, source->fpol, &f);

  b[1] = f/x[0];
  
  return FIO_SUCCESS;
}
