#ifndef COMPOUND_FIELD_H
#define COMPOUND_FIELD_H

#include <deque>
#include "fusion_io_field.h"

class fio_compound_field {
  typedef std::deque<fio_field*> field_list;
  field_list fields;

 public:
  int eval(const double*, double*);
  int add_field(fio_field*);
};

#endif
