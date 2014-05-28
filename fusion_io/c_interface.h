#ifndef FIO_C_INTERFACE_H
#define FIO_C_INTERFACE_H

extern "C" {
  int fio_add_field(const int, const int, const int, const double);
  int fio_close_field(const int);
  int fio_close_series(const int);
  int fio_close_source(const int);
  int fio_create_compound_field(int*);
  int fio_eval_field(const int, const double*, double*);
  int fio_eval_field_deriv(const int, const double*, double*);
  int fio_eval_series(const int, const double, double*);
  int fio_get_options(const int);
  int fio_get_available_fields(const int, int*, field_type**);
  int fio_get_field(const int, const int, int*);
  int fio_get_series(const int, const int, int*);
  int fio_open_source(const int, const char*, int*);
  int fio_set_int_option(const int, const int);
  int fio_set_str_option(const int, const char*);
  int fio_set_real_option(const int, const double);
}

#endif
