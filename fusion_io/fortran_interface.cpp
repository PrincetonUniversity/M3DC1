#include "fusion_io.h"

extern "C" {
  void fio_add_field_(const int*, const int*, const int*, const double*, int*);
  void fio_close_field_(const int*, int*);
  void fio_close_source_(const int*, int*);
  void fio_create_compound_field_(int*, int*);
  void fio_eval_field_(const int*, const double*, double*, int*);
  void fio_get_options_(const int*, int*);
  void fio_get_field_(const int*, const int*, int*, int*);
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

void fio_get_field_(const int* isrc, const int* itype, 
		int* handle, int* ierr)
{
  *ierr = fio_get_field(*isrc, *itype, handle);
}

void fio_get_options_(const int* isrc, int* ierr)
{
  *ierr = fio_get_options(*isrc);
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


/*
void fio_add_field_(const int* icfield, const int* ifield, 
		const int* op, const double* fac, int* ierr)
{
  fio_field* f = field_list[*icfield];
  if(typeid(*f) != typeid(fio_compound_field)) {
    std::cerr << "Error: can only add a field to a compound field"
	      << std::endl;
    *ierr = 1;
    return;
  }
  *ierr = ((fio_compound_field*)f)->add_field(field_list[*ifield], *op, *fac);
}

void fio_close_field_(const int* ifield, int* ierr)
{
  if(field_list[*ifield])
    delete(field_list[*ifield]);

  field_list[*ifield] = (fio_field*)0;

  *ierr = FIO_SUCCESS;
}

void fio_close_source_(const int* handle, int* ierr)
{
  *ierr = fio_close_source(&(source_list[*handle]));
}

void fio_create_compound_field_(int* ifield, int* ierr)
{
  fio_field* f = new fio_compound_field();

  *ifield = field_list.size();
  field_list.push_back(f);
  *ierr = FIO_SUCCESS;
}

void fio_eval_field_(const int* ifield, const double* x, double* v, int* ierr)
{
  *ierr = field_list[*ifield]->eval(x, v);
}

void fio_get_field_(const int* src, const int* type, 
		int* handle, int* ierr)
{
  fio_field* f;
  *ierr = source_list[*src]->get_field(*type, &f, &options);

  if(*ierr == FIO_SUCCESS) {
    *handle = field_list.size();
    field_list.push_back(f);
  }
}

void fio_get_options_(const int* src, int* ierr)
{
  *ierr = source_list[*src]->get_field_options(&options);
}

void fio_open_source_(const int* type, const char* filename, 
		  int* handle, int* ierr)
{
  fio_source* src;

  *ierr = fio_open_source(&src, *type, filename);

  if(*ierr != FIO_SUCCESS) return;
  
  *handle = source_list.size();
  source_list.push_back(src);
  *ierr = FIO_SUCCESS;
}

void fio_set_int_option_(const int* iopt, const int* v, int* ierr)
{
  *ierr = options.set_option(*iopt, *v);
}
void fio_set_str_option_(const int* iopt, const char* v, int* ierr)
{
  std::string str(v);
  *ierr = options.set_option(*iopt, str);
}
void fio_set_real_option_(const int* iopt, const double* v, int* ierr)
{
  *ierr = options.set_option(*iopt, *v);
}

*/
