#include "fusion_io.h"


int gato_magnetic_field::eval(const double* x, double* b)
{
  int e, ierr;
  double psi, f;

  ierr = source->interpolate_psi(x, &psi, &e);
  if(ierr != FIO_SUCCESS) return ierr;

  if(part != FIO_PERTURBED_ONLY) {
    b[0] = -source->dpsidz[e]/x[0];
    b[2] =  source->dpsidr[e]/x[0];
    
    ierr = source->interpolate_flux_function(source->ftor, psi, &f);
    if(ierr != FIO_SUCCESS) return ierr;
    
    b[1] =  f/x[0];
  }
  
  return FIO_SUCCESS;
}

int gato_pressure_field::eval(const double* x, double* p)
{
  double psi, z;
  int ierr;

  ierr = source->interpolate_psi(x, &psi);
  if(ierr != FIO_SUCCESS) return ierr;

  if(part != FIO_PERTURBED_ONLY) {
    ierr = source->interpolate_flux_function(source->pressure, psi, p);
    if(ierr != FIO_SUCCESS) return ierr;
  }
  return FIO_SUCCESS;
}

