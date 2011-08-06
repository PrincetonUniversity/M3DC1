#include <fftw3.h>

#include <stdio.h>

void fftw_(double* in, fftw_complex* ou, int* N)
{
  fftw_plan p;
  int i;

  fftw_complex* out;

  out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(*N));
  
  p = fftw_plan_dft_r2c_1d(*N, in, out, FFTW_ESTIMATE);

  fftw_execute(p);

  for(i=0; i<*N/2+1; i++) {
    printf("FFTW result: out[%d] = (%g, %g)\n", i, out[i][0], out[i][1]);
  }

  fftw_destroy_plan(p);
  fftw_free(out);
}
