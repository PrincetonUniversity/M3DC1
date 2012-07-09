#ifndef GEQDSK_FIELD_H
#define GEQDSK_FIELD_H

#include "geqdsk_source.h"

class geqdsk_field : public fio_field {
 protected:
  geqdsk_source* source;

 public:
  geqdsk_field(geqdsk_source* s) 
    { source = s; }
  virtual ~geqdsk_field()
    { }

  virtual int dimension() const = 0;
  virtual int eval(const double*, double*) = 0;
};

class geqdsk_current_density : public geqdsk_field {
 public:
  geqdsk_current_density(geqdsk_source* s)
    : geqdsk_field(s) { }
  geqdsk_current_density* clone() const 
  { return new geqdsk_current_density(*this); }
  int dimension() const { return 3; }
  int eval(const double*, double*);
};

class geqdsk_magnetic_field : public geqdsk_field {
 public:
  geqdsk_magnetic_field(geqdsk_source* s)
    : geqdsk_field(s) { }
  geqdsk_magnetic_field* clone() const 
  { return new geqdsk_magnetic_field(*this); }
  int dimension() const { return 3; }
  int eval(const double*, double*);
};

class geqdsk_pressure_field : public geqdsk_field {
 public:
  geqdsk_pressure_field(geqdsk_source* s)
    : geqdsk_field(s) { }
  geqdsk_pressure_field* clone() const 
  { return new geqdsk_pressure_field(*this); }
  int dimension() const { return 1; }
  int eval(const double*, double*);
};


#endif
