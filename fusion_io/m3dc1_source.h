#ifndef M3DC1_SOURCE_H
#define M3DC1_SOURCE_H

#include <m3dc1_file.h>

#include "fusion_io_source.h"

class m3dc1_source : public fio_source {
  m3dc1_file file;

 public:
  int open(const char*);
  int close();

  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);
};


#endif
