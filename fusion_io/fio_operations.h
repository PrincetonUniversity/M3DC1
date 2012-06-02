#ifndef FIO_OPERATIONS_H
#define FIO_OPERATIONS_H

#include "fusion_io_field.h"

class fio_field_add : public fio_field {
  fio_field* f;
  double val;
 public:
  fio_field_add(const fio_field*, const double);
  fio_field_add* clone() const { return (fio_field_add*)this; }
  virtual int dimension() const;
  virtual int eval(const double*, double*);
};

class fio_field_multiply : public fio_field {
  fio_field* f;
  double val;
 public:
  fio_field_multiply(const fio_field*, const double);
  fio_field_multiply* clone() const { return (fio_field_multiply*)this; }
  virtual int dimension() const;
  virtual int eval(const double*, double*);
};

class fio_field_sum : public fio_field {
  fio_field* f[2];
  int dim;
 public:
  fio_field_sum(const fio_field*, const fio_field*);
  fio_field_sum* clone() const { return (fio_field_sum*)this; }
  ~fio_field_sum();
  int dimension() const;
  int eval(const double*, double*);
};


class fio_field_product : public fio_field {
  fio_field* f[2];
  int dim;
 public:
  fio_field_product(const fio_field*, const fio_field*);
  ~fio_field_product();
  fio_field_product* clone() const { return (fio_field_product*)this; }
  int dimension() const;
  int eval(const double*, double*);
};

#endif
