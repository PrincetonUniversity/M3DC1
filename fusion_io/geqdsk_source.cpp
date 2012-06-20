#include "fusion_io.h"

#include <iostream>
#include <fstream>
#include <iomanip>

geqdsk_source::geqdsk_source()
  : psi(0), psirz(0), fpol(0)
{
}

geqdsk_source::~geqdsk_source()
{
  if(psirz) delete[] psirz;
  if(fpol) delete[] fpol;
  if(psi) delete[] psi;
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
  fpol = new double[nw];  

  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> fpol[i];
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // press
  for(int i=0; i<nw; i++) gfile >> std::setw(16) >> dum;      // ffprime
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
  fields->push_back(FIO_MAGNETIC_FIELD);

  return FIO_SUCCESS;
}

int geqdsk_source::get_field(const field_type t,fio_field** f,
			     const fio_option_list* opt)
{
  *f = 0;
  fio_field* mf;
  bool unneeded_species = false;
  int s, result;

  result = FIO_SUCCESS;

  opt->get_option(FIO_SPECIES, &s);

  switch(t) {
  case(FIO_MAGNETIC_FIELD):
    mf = new geqdsk_magnetic_field(this);
    if(s!=0) unneeded_species = true;
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  if(result==FIO_BAD_SPECIES) {
    std::cerr << "Unsupported species: " << fio_species(s).name() << std::endl;
    return result;
  }
  
  *f = mf;
  return result;
}
