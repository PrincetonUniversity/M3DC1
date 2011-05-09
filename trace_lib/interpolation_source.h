#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "trace_source.h"

class interpolation_source: public trace_field_source {

 private:
  double r_extent[2], z_extent[2];
  int n[3];
  double* Br;
  double* Bphi;
  double* Bz;

  void free_grid();
  bool trilinear_interpolation(const double[3], const int[3],
			       double*, double*, double*) const;
  bool triquadratic_interpolation(const double[3], const int[3],
				  double*, double*, double*) const;
  bool tricubic_interpolation(const double[3], const int[3],
			      double*, double*, double*) const;


 public:
  trace_source_list sources;

  interpolation_source()
    {
      Br = Bphi = Bz = 0;
    }
  virtual ~interpolation_source();

  virtual bool load();
  virtual bool center(double*, double*) const;
  virtual bool extent(double*, double*, double*, double*) const;
  virtual bool eval(const double, const double, const double,
		    double*, double*, double*);

  void set_extent(const double, const double,
		  const double, const double,
		  const int, const int, const int);
};

#endif
