#include "fusion_io.h"
#include <iostream>

int fio_open_source(fio_source** src, const int type, const char* filename)
{
  int ierr;
  *src = 0;

  switch(type) {
  case(FIO_GEQDSK_SOURCE):
    *src = new geqdsk_source();
    ierr = (*src)->open(filename);
    break;

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

int fio_get_field_name(const field_type f, std::string* s)
{
  switch(f) {
  case(FIO_ALPHA):             *s = "alpha";            break;
  case(FIO_CURRENT_DENSITY):   *s = "current density";  break;
  case(FIO_DENSITY):           *s = "density";          break;
  case(FIO_ELECTRIC_FIELD):    *s = "electric field";   break;
  case(FIO_FLUID_VELOCITY):    *s = "fluid velocity";   break; 
  case(FIO_MAGNETIC_FIELD):    *s = "magnetic field";   break;
  case(FIO_PRESSURE):          *s = "pressure";         break;
  case(FIO_SCALAR_POTENTIAL):  *s = "scalar potential"; break;
  case(FIO_TOTAL_PRESSURE):    *s = "total pressure";   break;
  case(FIO_VECTOR_POTENTIAL):  *s = "vector potential"; break;
  case(FIO_VELOCITY):          *s = "velocity";         break; 
  default:
    *s = "Unnamed field";
    return 1;
  }

  return FIO_SUCCESS; 
}
