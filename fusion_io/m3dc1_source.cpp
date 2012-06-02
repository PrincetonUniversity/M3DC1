#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "fusion_io_defs.h"

int m3dc1_source::open(const char* filename)
{
  if(!file.open(filename))
    return FIO_FILE_ERROR;
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

  return FIO_SUCCESS;
}

int m3dc1_source::get_field(const field_type t,fio_field** f,const fio_option_list* opt)
{
  *f = 0;
  m3dc1_fio_field* mf;

  switch(t) {
  case(FIO_MAGNETIC_FIELD):
    mf = new m3dc1_magnetic_field(this);
    break;

  case(FIO_PRESSURE):
    mf = new m3dc1_scalar_field(this,"P");
    break;

  case(FIO_DENSITY):
    mf = new m3dc1_scalar_field(this,"den");
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  int result = mf->load(opt);
  if(result == FIO_SUCCESS) {
    *f = mf;
  } else {
    delete(mf);
  }
  return result;
}
