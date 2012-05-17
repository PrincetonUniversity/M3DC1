#ifndef FUSION_IO_FIELD
#define FUSION_IO_FIELD

class fio_field {
 public:
  virtual ~fio_field()
    { }

  virtual int get_dimensions() const = 0;
  virtual int eval(const double*, double*) = 0;
};


#endif
