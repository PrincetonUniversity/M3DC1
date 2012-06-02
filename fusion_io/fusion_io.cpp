#include "fusion_io.h"

int fio_open_source(fio_source** src, const int type, const char* filename)
{
  int ierr;
  *src = 0;

  switch(type) {
  case(FIO_M3DC1_SOURCE):
    *src = new m3dc1_source();
    ierr = (*src)->open(filename);
    break;

  default:
    std::cerr << "Source type " << type << " unsupported." << std::endl;
    return FIO_UNSUPPORTED;
  }

  if(ierr != FIO_SUCCESS) {
    delete(*src);
    return ierr;
  }
  
  return FIO_SUCCESS;
}

int fio_close_source(fio_source** source)
{
  if(*source) {
    int result = (*source)->close();
    delete(*source);
    *source = 0;
    return result;
  } else return 1;
}

int fio_close_field(fio_field** field)
{
  if(*field) {
    delete(*field);
    *field = 0;
    return FIO_SUCCESS;
  } else return 1;
}
