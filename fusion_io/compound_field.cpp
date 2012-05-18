#include "compound_field.h"
#include "fusion_io_defs.h"

#include <iostream>

fio_compound_field::
component_field::component_field(fio_field* f, const int o, const int d)
{
  field = f;
  op = o;
  factor = d;
}

fio_compound_field::
component_field::component_field(const component_field& f)
{
  field = f.field;
  op = f.op;
  factor = f.factor;
}


fio_compound_field::fio_compound_field()
{
  dim = 0;
}

int fio_compound_field::dimension() const
{
  return dim;
}

int fio_compound_field::add_field(fio_field* f, const int op, const double d)
{
  if(fields.empty()) {
    dim = f->dimension();
  } else {
    if(f->dimension() != dimension()) {
      std::cerr << "Error: adding field of incompatible dimension" 
		<< std::endl;
      return 1;
    }
  }

  fields.push_back(component_field(f, op, d));
  return FIO_SUCCESS;
}

int fio_compound_field::eval(const double* x, double* v)
{
  double z[dimension()];
  for(int j=0; j<dimension(); j++)
    v[j] = 0;

  field_list::iterator i = fields.begin();
  while(i != fields.end()) {
    int result = i->field->eval(x, z);
    for(int j=0; j<dimension(); j++)
      switch(i->op) {
      case(FIO_ADD):       v[j] += i->factor*z[j]; break;
      case(FIO_MULTIPLY):  v[j] *= i->factor*z[j]; break;
      case(FIO_DIVIDE):    v[j] /= i->factor*z[j]; break;
      }
      
    if(result != FIO_SUCCESS) return result;

    i++;
  }

  return FIO_SUCCESS;
}
