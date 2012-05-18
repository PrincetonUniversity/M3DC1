#include "fusion_io.h"

#include <deque>

static std::deque<fio_source*> source_list;
static std::deque<fio_field*> field_list;
static fio_option_list options;

extern "C" {
  void add_field_(const int*, const int*, const int*, const double*, int*);
  void close_field_(const int*, int*);
  void close_source_(const int*, int*);
  void create_compound_field_(int*, int*);
  void eval_field_(const int*, const double*, double*, int*);
  void get_options_(const int*, int*);
  void get_field_(const int*, const int*, int*, int*);
  void open_source_(const int*, const char*, int*, int*);
  void set_int_option_(const int*, const int*, int*);
  void set_str_option_(const int*, const char*, int*);
  void set_real_option_(const int*, const double*, int*);
}

void add_field_(const int* icfield, const int* ifield, 
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

void close_field_(const int* ifield, int* ierr)
{
  if(field_list[*ifield])
    delete(field_list[*ifield]);

  field_list[*ifield] = (fio_field*)0;

  *ierr = FIO_SUCCESS;
}

void close_source_(const int* handle, int* ierr)
{
  if(source_list[*handle]) 
    delete(source_list[*handle]);

  source_list[*handle] = (fio_source*)0;
  *ierr = FIO_SUCCESS;
}

void create_compound_field_(int* ifield, int* ierr)
{
  fio_field* f = new fio_compound_field();

  *ifield = field_list.size();
  field_list.push_back(f);
  *ierr = FIO_SUCCESS;
}

void eval_field_(const int* ifield, const double* x, double* v, int* ierr)
{
  *ierr = field_list[*ifield]->eval(x, v);
}

void get_field_(const int* src, const int* type, 
		int* handle, int* ierr)
{
  fio_field* f;
  *ierr = source_list[*src]->get_field(*type, &f, &options);

  if(*ierr == FIO_SUCCESS) {
    *handle = field_list.size();
    field_list.push_back(f);
  }
}

void get_options_(const int* src, int* ierr)
{
  *ierr = source_list[*src]->get_field_options(&options);
}

void open_source_(const int* type, const char* filename, 
		  int* handle, int* ierr)
{
  fio_source* src;

  switch(*type) {
  case(FIO_M3DC1_SOURCE):
    src = new m3dc1_source();
    *ierr = src->open(filename);
    if(*ierr != FIO_SUCCESS) {
      delete(src);
      return;
    }
    break;

  default:
    std::cerr << "Source type " << type << " unsupported." << std::endl;
    *ierr = FIO_UNSUPPORTED;
    return;
  }
  
  *handle = source_list.size();
  source_list.push_back(src);
  *ierr = FIO_SUCCESS;
}

void set_int_option_(const int* iopt, const int* v, int* ierr)
{
  *ierr = options.set_option(*iopt, *v);
}
void set_str_option_(const int* iopt, const char* v, int* ierr)
{
  std::string str(v);
  *ierr = options.set_option(*iopt, str);
}
void set_real_option_(const int* iopt, const double* v, int* ierr)
{
  *ierr = options.set_option(*iopt, *v);
}

