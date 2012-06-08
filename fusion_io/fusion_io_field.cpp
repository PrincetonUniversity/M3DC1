#include "fusion_io.h"

fio_field& fio_field::operator+(const fio_field& f)
{
  return *(new fio_field_sum(this, &f));
}

fio_field& fio_field::operator*(const fio_field& f)
{
  return *(new fio_field_product(this, &f));
}
