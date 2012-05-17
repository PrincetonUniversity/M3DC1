#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "fusion_io.h"

int m3dc1_source::open(const char* filename)
{
  if(!file.open(filename))
    return 1;
  return 0;
}

int m3dc1_source::close()
{
  if(!file.close())
    return 1;
  return 0;
}

int m3dc1_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  return 0;
}

int m3dc1_source::get_field(const field_type t,fio_field** f,const fio_option_list* opt)
{
  *f = 0;
  m3dc1_fio_field* mf;

  switch(t) {
  case(FIO_MAGNETIC_FIELD):
    mf = new m3dc1_magnetic_field();
    break;

  case(FIO_PRESSURE):
    mf = new m3dc1_scalar_field("P");
    break;

  case(FIO_DENSITY):
    mf = new m3dc1_scalar_field("den");
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  int result = mf->load(&file, opt);
  if(result == FIO_SUCCESS) {
    *f = mf;
  } else {
    delete(mf);
  }
  return result;
}
