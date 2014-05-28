#include "fusion_io.h"

extern "C" {
  void fio_add_field_(const int*, const int*, const int*, const double*, int*);
  void fio_close_field_(const int*, int*);
  void fio_close_series_(const int*, int*);
  void fio_close_source_(const int*, int*);
  void fio_create_compound_field_(int*, int*);
  void fio_eval_field_(const int*, const double*, double*, int*);
  void fio_eval_field_deriv_(const int*, const double*, double*, int*);
  void fio_eval_series_(const int*, const double*, double*, int*);
  void fio_get_options_(const int*, int*);
  void fio_get_field_(const int*, const int*, int*, int*);
  void fio_get_series_(const int*, const int*, int*, int*);
  void fio_open_source_(const int*, const char*, int*, int*);
  void fio_set_int_option_(const int*, const int*, int*);
  void fio_set_str_option_(const int*, const char*, int*);
  void fio_set_real_option_(const int*, const double*, int*);
}

void fio_add_field_(const int* icfield, const int* ifield, 
		const int* op, const double* fac, int* ierr)
{
  *ierr = fio_add_field(*icfield, *ifield, *op, *fac);
}

void fio_close_field_(const int* ifield, int* ierr)
{
  *ierr = fio_close_field(*ifield);
}

void fio_close_series_(const int* iseries, int* ierr)
{
  *ierr = fio_close_series(*iseries);
}

void fio_close_source_(const int* isrc, int* ierr)
{
  *ierr = fio_close_source(*isrc);
}

void fio_create_compound_field_(int* ifield, int* ierr)
{
  *ierr = fio_create_compound_field(ifield);
}

void fio_eval_field_(const int* ifield, const double* x, double* v, int* ierr)
{
  *ierr = fio_eval_field(*ifield, x, v);
}

void fio_eval_field_deriv_(const int* ifield, const double* x, double* v, 
			   int* ierr)
{
  *ierr = fio_eval_field_deriv(*ifield, x, v);
}

void fio_eval_series_(const int* iseries, const double* x, double* v, int* ierr)
{
  *ierr = fio_eval_series(*iseries, *x, v);
}

void fio_get_field_(const int* isrc, const int* itype,
		    int* handle, int* ierr)
{
  *ierr = fio_get_field(*isrc, *itype, handle);
}

void fio_get_options_(const int* isrc, int* ierr)
{
  *ierr = fio_get_options(*isrc);
}

void fio_get_series_(const int* isrc, const int* itype,
		    int* handle, int* ierr)
{
  *ierr = fio_get_series(*isrc, *itype, handle);
}

void fio_open_source_(const int* type, const char* filename, 
		  int* handle, int* ierr)
{
  *ierr = fio_open_source(*type, filename, handle);
}

void fio_set_int_option_(const int* iopt, const int* v, int* ierr)
{
  *ierr = fio_set_int_option(*iopt, *v);
}
void fio_set_str_option_(const int* iopt, const char* v, int* ierr)
{
  *ierr = fio_set_str_option(*iopt, v);
}
void fio_set_real_option_(const int* iopt, const double* v, int* ierr)
{
  *ierr = fio_set_real_option(*iopt, *v);
}
