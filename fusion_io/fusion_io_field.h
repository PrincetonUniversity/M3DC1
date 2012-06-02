#ifndef FUSION_IO_FIELD_H
#define FUSION_IO_FIELD_H

class fio_field {
 public:
  virtual ~fio_field()
    { }

  virtual fio_field* clone() const = 0;
  virtual int dimension() const = 0;
  virtual int eval(const double*, double*) = 0;

  fio_field& operator+(const fio_field&);
  fio_field& operator*(const fio_field&);
};

#endif
