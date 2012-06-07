#ifndef FUSION_IO_SOURCE_H
#define FUSION_IO_SOURCE_H

#include "fusion_io_species.h"
#include "fusion_io_field.h"
#include "options.h"

typedef int field_type;

class fio_source {
 public:
  virtual ~fio_source()
    { }
  virtual int open(const char*) = 0;
  virtual int close() = 0;

  virtual int get_field_options(fio_option_list*) const = 0;
  virtual int get_field(const field_type, fio_field**, const fio_option_list*, const fio_species*) = 0;
};

#endif 
