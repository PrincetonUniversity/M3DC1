#include <gsl/gsl_sf_bessel.h>

void gsl_bessel_i_(const int* n, const double* x, double *out)
{
  *out = gsl_sf_bessel_In(*n, *x);
}

void gsl_bessel_j_(const int* n, const double* x, double *out)
{
  *out = gsl_sf_bessel_Jn(*n, *x);
}

void gsl_bessel_y_(const int* n, const double* x, double *out)
{
  *out = gsl_sf_bessel_Yn(*n, *x);
}

void gsl_erf_(const double* x, double *out)
{
  *out = gsl_sf_erf(*x);
}

