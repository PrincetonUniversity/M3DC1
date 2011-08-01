#include <m3dc1_source.h>

static m3dc1_source* source0;
static m3dc1_source* source1;

extern "C" void m3dc1_load_file_(int* time,  int* ierr)
{
  source0 = new m3dc1_source("C1.h5", -1);
  if(!source0->load()) {
    *ierr = 1;
    return;
  }

  if(*time >= 0) {
    source1 = new m3dc1_source("C1.h5", *time);
    if(!source1->load()) {
      delete(source0);
      *ierr = 1;
      return;
    }
  } else source1 = 0;

  *ierr = 0;
}

extern "C" void m3dc1_unload_file_()
{
  delete(source0);
  if(source1) delete(source1);
}

extern "C" void m3dc1_get_field_(const double* R, const double* Phi, const double* Z, 
				 double* Br, double* Bphi, double* Bz)
{
  *Br = 0.;
  *Bz = 0.;
  *Bphi = 0.;
  source0->eval(*R, *Phi, *Z, Br, Bphi, Bz);
  if(source1) source1->eval(*R, *Phi, *Z, Br, Bphi, Bz);
}

