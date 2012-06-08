#ifndef GEQDSK_FIELD_H
#define GEQDKS_FIELD_H

#include "geqdsk_source.h"

class geqdsk_field : public fio_field {
 protected:
  geqdsk_source* source;

 public:
  geqdsk_field(geqdsk_source* s) 
    { source = s; }
  virtual ~geqdsk_field()
    { }

  virtual int dimension() const;
  virtual int eval(const double*, double*);
};


class geqdsk_magnetic_field : public geqdsk_field {
 public:
  geqdsk_magnetic_field(geqdsk_source* s)
    : geqdsk_field(s) { }
  int dimension() const { return 3; }
  int eval(const double*, double*);
};

#endif
