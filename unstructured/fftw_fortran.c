#include <fftw3.h>

#include <stdio.h>

static fftw_plan p;

void fftw_init_(double* in, fftw_complex* out, int* N)
{
  p = fftw_plan_dft_r2c_1d(*N, in, out, FFTW_ESTIMATE);
}

void fftw_run_()
{
  fftw_execute(p);
}

void fftw_destroy_()
{
  fftw_destroy_plan(p);
}

void fftw_(double* in, fftw_complex* out, int* N)
{
  p = fftw_plan_dft_r2c_1d(*N, in, out, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}
