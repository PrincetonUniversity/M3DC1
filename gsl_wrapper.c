#include <gsl/gsl_sf_bessel.h>

void gsl_bessel_i_(const int* n, const double* x, double *out)
{
  *out = gsl_sf_bessel_In(*n, *x);
}
