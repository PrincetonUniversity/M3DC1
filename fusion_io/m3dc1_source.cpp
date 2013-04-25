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
  file.read_parameter("bzero", &bzero);
  file.read_parameter("rzero", &rzero);
  file.read_parameter("n0_norm", &n0);
  file.read_parameter("l0_norm", &L0);
  file.read_parameter("b0_norm", &B0);

  const double c = 3.e10;
  const double m_p = 1.6726e-24;

  // define some normalization quantities
  p0 = B0*B0/(4.*M_PI);
  J0 = B0*c/(4.*M_PI*L0);
  v0 = B0/sqrt(4.*M_PI*ion_mass*m_p*n0);

  // convert normalization quantities to mks
  n0 /= 1.e-6;
  B0 /= 1e4;
  L0 /= 100.;
  p0 /= 10.;
  J0 /= (c*1e-5);
  v0 /= 100.;

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
  opt->add_option(FIO_PHASE, 0.);
  opt->add_option(FIO_PART, FIO_TOTAL);
  opt->add_option(FIO_SPECIES, 0);

  return FIO_SUCCESS;
}

int m3dc1_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_DENSITY);
  fields->push_back(FIO_CURRENT_DENSITY);
  fields->push_back(FIO_FLUID_VELOCITY);
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_PRESSURE);
  fields->push_back(FIO_TOTAL_PRESSURE);
  fields->push_back(FIO_VECTOR_POTENTIAL);

  return FIO_SUCCESS;
}

int m3dc1_source::get_field(const field_type t,fio_field** f,
			    const fio_option_list* opt)
{
  *f = 0;
  m3dc1_fio_field* mf;
  bool unneeded_species = false;
  int s, result;

  result = FIO_SUCCESS;

  opt->get_option(FIO_SPECIES, &s);
  if(s==FIO_MAIN_ION) s = ion_species;

  switch(t) {
  case(FIO_DENSITY):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "den", n0*zeff);
    } else if(s==ion_species) {
      mf = new m3dc1_scalar_field(this, "den", n0);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  case(FIO_CURRENT_DENSITY):
    mf = new m3dc1_current_density(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_ELECTRIC_FIELD):
    mf = new m3dc1_electric_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_FLUID_VELOCITY):
    mf = new m3dc1_fluid_velocity(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_GRAD_VECTOR_POTENTIAL):
    mf = new m3dc1_grad_vector_potential(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_MAGNETIC_FIELD):
    mf = new m3dc1_magnetic_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_PRESSURE):
    if(s==fio_electron) {
      mf = new m3dc1_scalar_field(this, "Pe", p0);
    } else if(s==ion_species) {
      mf = new m3dc1_pi_field(this);
    } else {
      result = FIO_BAD_SPECIES;
    }
    break;

  case(FIO_SCALAR_POTENTIAL):
    mf = new m3dc1_phi_field(this);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_TOTAL_PRESSURE):
    mf = new m3dc1_scalar_field(this, "P", p0);
    if(s!=0) unneeded_species = true;
    break;

  case(FIO_VECTOR_POTENTIAL):
    mf = new m3dc1_vector_potential(this);
    if(s!=0) unneeded_species = true;
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
