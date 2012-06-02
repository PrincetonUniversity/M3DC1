#include "fusion_io.h"

#include <deque>

static std::deque<fio_source*> source_list;
static std::deque<fio_field*> field_list;
static fio_option_list options;

int fio_add_field(const int icfield, const int ifield, 
		const int op, const double fac)
{
  fio_field* f = field_list[icfield];
  if(typeid(*f) != typeid(fio_compound_field)) {
    std::cerr << "Error: can only add a field to a compound field"
	      << std::endl;
    return 1;
  }
  return ((fio_compound_field*)f)->add_field(field_list[ifield], op, fac);
}

int fio_close_field(const int ifield)
{
  if(field_list[ifield])
    delete(field_list[ifield]);

  field_list[ifield] = (fio_field*)0;

  return FIO_SUCCESS;
}

int fio_close_source(const int ifield)
{
  return fio_close_source(&(source_list[ifield]));
}

int fio_create_compound_field(int* ifield)
{
  fio_field* f = new fio_compound_field();

  *ifield = field_list.size();
  field_list.push_back(f);
  return FIO_SUCCESS;
}

int fio_eval_field(const int ifield, const double* x, double* v)
{
  return field_list[ifield]->eval(x, v);
}

int fio_get_field(const int isrc, const int type, int* handle)
{
  fio_field* f;
  int ierr;

  ierr = source_list[isrc]->get_field(type, &f, &options);

  if(ierr == FIO_SUCCESS) {
    *handle = field_list.size();
    field_list.push_back(f);
  }
  return ierr;
}

int fio_get_options(const int isrc)
{
  return source_list[isrc]->get_field_options(&options);
}

int fio_open_source(const int itype, const char* filename, int* handle)
{
  fio_source* src;
  int ierr;

  ierr = fio_open_source(&src, itype, filename);
  if(ierr != FIO_SUCCESS) return ierr;
  
  *handle = source_list.size();
  source_list.push_back(src);
  return FIO_SUCCESS;
}

int fio_set_int_option(const int iopt, const int v)
{
  return options.set_option(iopt, v);
}
int fio_set_str_option(const int iopt, const char* v)
{
  std::string str(v);
  return options.set_option(iopt, str);
}
int fio_set_real_option(const int iopt, const double v)
{
  return options.set_option(iopt, v);
}

