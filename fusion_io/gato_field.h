#ifndef GATO_FIELD_H
#define GATO_FIELD_H

#include "gato_source.h"

class gato_field : public fio_field {
 protected:
  gato_source* source;
  int part;

 public:
  gato_field(gato_source* s, int p) 
    { source = s; part = p; }
  virtual ~gato_field()
    { }

  virtual int dimension() const = 0;
  virtual int eval(const double*, double*) = 0;
};

class gato_magnetic_field : public gato_field {
 public:
 gato_magnetic_field(gato_source* s, int p)
   : gato_field(s, p) { }
  gato_magnetic_field* clone() const 
  { return new gato_magnetic_field(*this); }
  int dimension() const { return 3; }
  int eval(const double*, double*);
};

class gato_pressure_field : public gato_field {
 public:
 gato_pressure_field(gato_source* s, int p)
   : gato_field(s, p) { }
  gato_pressure_field* clone() const 
  { return new gato_pressure_field(*this); }
  int dimension() const { return 1; }
  int eval(const double*, double*);
};


#endif
