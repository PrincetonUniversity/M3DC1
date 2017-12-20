#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_ellint.h>

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

void gsl_elliptic_e_(const double* k, double *out)
{
  *out = gsl_sf_ellint_Ecomp(*k, GSL_PREC_SINGLE);
}

void gsl_elliptic_k_(const double* k, double *out)
{
  *out = gsl_sf_ellint_Kcomp(*k, GSL_PREC_SINGLE);
}

void gsl_erf_(const double* x, double *out)
{
  *out = gsl_sf_erf(*x);
}

