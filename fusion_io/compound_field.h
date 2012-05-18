#ifndef COMPOUND_FIELD_H
#define COMPOUND_FIELD_H

#include <deque>
#include "fusion_io_field.h"
#include "fusion_io_defs.h"

class fio_compound_field : public fio_field {
  struct component_field {
    fio_field* field;
    int        op;
    double     factor;
    component_field(fio_field* f, const int o, const int d);
    component_field(const component_field& f);
  };
  typedef std::deque<component_field> field_list;
  field_list fields;
  int dim;

 public:
  fio_compound_field();
  int dimension() const;
  int eval(const double*, double*);
  int add_field(fio_field*, const int op=FIO_ADD, const double d=1.);
};

#endif
