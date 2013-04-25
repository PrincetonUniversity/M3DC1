#ifndef FUSION_IO_H
#define FUSION_IO_H

#include "fusion_io_defs.h"
#include "fusion_io_species.h"
#include "fusion_io_source.h"
#include "fusion_io_field.h"
#include "fio_operations.h"
#include "compound_field.h"
#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "geqdsk_source.h"
#include "geqdsk_field.h"
#include "gato_source.h"
#include "gato_field.h"
#include "c_interface.h"

#include <string>

int fio_close_source(fio_source** source);
int fio_close_field(fio_field** field);
int fio_get_field_name(field_type, std::string*);
int fio_get_option_name(const int, std::string*);
int fio_open_source(fio_source** src, const int type, const char* filename);

#endif
