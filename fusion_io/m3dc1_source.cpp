#include "fusion_io.h"
#include <iostream>

int m3dc1_source::open(const char* filename)
{
  if(!file.open(filename))
    return FIO_FILE_ERROR;

  // read some relevant data from the file
  file.read_parameter("extsubtract", &extsubtract);
  file.read_parameter("icomplex", &icomplex);
  file.read_parameter("3d", &i3d);
  file.read_parameter("eqsubtract", &eqsubtract);
  file.read_parameter("linear", &linear);
  file.read_parameter("zeff", &zeff);
  file.read_parameter("ion_mass", &ion_mass);
  file.read_parameter("n0_norm", &n0);
  file.read_parameter("l0_norm", &L0);
  file.read_parameter("b0_norm", &B0);

  // define some normalization quantities
  p0 = B0*B0/(4.*M_PI);

  // convert normalization quantities to mks
  n0 /= 1.e-6;
  B0 /= 1e4;
  L0 /= 100.;
  p0 /= 10.;

  // determine ion species (assume one proton and no electrons)
  ion_species = fio_species(ion_mass, 1, 0);

  return FIO_SUCCESS;
}

int m3dc1_source::close()
{
  if(!file.close())
    return FIO_FILE_ERROR;
  return FIO_SUCCESS;
}

int m3dc1_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  opt->add_option(FIO_TIMESLICE, 0);
  opt->add_option(FIO_LINEAR_SCALE, 1.);
  opt->add_option(FIO_PERTURBED_ONLY, 0);
  opt->add_option(FIO_SPECIES, 0);

  return FIO_SUCCESS;
}

int m3dc1_source::get_field(const field_type t,fio_field** f,
			    const fio_option_list* opt)
{
  *f = 0;
  m3dc1_fio_field* mf;
  bool unneeded_species = false;
  int s, result;

  opt->get_option(FIO_SPECIES, &s);
  if(s==FIO_MAIN_ION) s = ion_species;

  switch(t) {
  case(FIO_MAGNETIC_FIELD):
    mf = new m3dc1_magnetic_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_TOTAL_PRESSURE):
    mf = new m3dc1_scalar_field(this, "P", p0);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_PRESSURE):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "Pe", p0);
    } else if(s==ion_species) {
      mf = new m3dc1_scalar_field(this, "P", p0);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  case(FIO_DENSITY):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "den", n0*zeff);
    } else if(s==ion_species) {
      mf = new m3dc1_scalar_field(this, "den", n0);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  if(result==FIO_BAD_SPECIES) {
    std::cerr << "Unsupported species: " << fio_species(s).name() << std::endl;
    std::cerr << "Main ions: " << ion_species.name() << std::endl;
    return result;
  }

  result = mf->load(opt);
  if(result == FIO_SUCCESS) {
    *f = mf;
  } else {
    delete(mf);
  }
  return result;
}
