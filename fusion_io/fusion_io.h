#ifndef FUSION_IO_H
#define FUSION_IO_H

#include "fusion_io_defs.h"
#include "fusion_io_source.h"
#include "fusion_io_field.h"
#include "fio_operations.h"
#include "compound_field.h"
#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "geqdsk_source.h"
#include "geqdsk_field.h"
#include "c_interface.h"

int fio_open_source(fio_source** src, const int type, const char* filename);
int fio_close_source(fio_source** source);
int fio_close_field(fio_field** field);

#endif
