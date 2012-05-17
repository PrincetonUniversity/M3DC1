#include "compound_field.h"

int fio_compound_field::add_field(fio_field* f)
{
  fields.push_back(f);
  return 0;
}

int fio_compound_field::eval(const double* x, double* v)
{
  *v = 0;

  field_list::iterator i = fields.begin();
  while(i != fields.end()) {
    int result = (*i)->eval(x, v);
    if(result != 0) return result;

    i++;
  }

  return 0;
}
