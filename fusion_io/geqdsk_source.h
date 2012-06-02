#ifndef GEQDSK_SOURCE_H
#define GEQDSK_SOURCE_H

#include "fusion_io_source.h"

#include <string>

class geqdsk_source : public fio_source {
 public:
  double rmaxis, zmaxis;

  int nw, nh;
  double dx, dz;
  double rleft, zmid;

  double* psi;
  double* psirz;
  double* fpol; 

 public:
  std::string filename;
  
  geqdsk_source();
  virtual ~geqdsk_source();

  virtual int get_field_options(fio_option_list*);
  virtual int get_field(const field_type, fio_field**, const fio_option_list*);

  int open(const char*);
  int close();

  int center(double* r0, double* z0) const;
  int extent(double* r0, double* r1, double* z0, double* z1) const;
};

#endif
