#ifndef GATO_SOURCE_H
#define GATO_SOURCE_H

#include "fusion_io_source.h"
#include <fstream>

class gato_source : public fio_source {
  int ntor;
  int jpsi, itht;
  double psimax, xma, zma;
  double paxe, faxe, ppace, ffpace, qace, fqpace, fqmaxe;
  double dnaxe, pfaxe, dnnorm, pfnorm;
  double psilim, fpli, qlim, psisep, xsep, zsep;
  double btor, totcur, omega0;
  double *r_bound, *t_bound;

 public:
  double *psival, *pressure, *ftor, *pprime, *ffprime;
  double *rcc, *zcc, *psimesh, *dpsidr, *dpsidz;

  static int scan_until(std::ifstream&, const char*);

  int set_element_bounds();
  int get_element(const double*, int*);

 public:
  int interpolate_psi(const double*, double*, int* e=0);
  int interpolate_flux_function(const double*, const double, double*) const;
  
 public:
  gato_source();
  virtual ~gato_source();
  virtual int open(const char*);
  virtual int close();

  virtual int get_field_options(fio_option_list*) const;
  virtual int get_field(const field_type, fio_field**, const fio_option_list*);
  virtual int get_available_fields(fio_field_list*) const;
};

#endif
