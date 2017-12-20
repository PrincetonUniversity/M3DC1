#include <fftw3.h>

#include <stdio.h>

void fftw_(double* in, fftw_complex* out, int* N)
{
  fftw_plan p;
  
  p = fftw_plan_dft_r2c_1d(*N, in, out, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}
