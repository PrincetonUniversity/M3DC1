#ifndef TRACE_H
#define TRACE_H

#include <deque>
#include <math.h>

class trace_field_source {
 public:
  bool interpolate;
  bool toroidal;

  trace_field_source() { interpolate = false; toroidal=true; }
  virtual ~trace_field_source() { }

  virtual bool load() = 0;

  virtual bool eval(const double r, const double phi, const double z, 
		    double* b_r, double* b_phi, double* b_z) = 0;
  virtual bool eval_psi(const double r, const double z, double* psi)
  { return false; }
  virtual bool psibound(double* psi_axis, double* psi_lcfs) const
  { return false; }

  virtual bool center(double* r0, double* z0) const = 0;
  virtual bool extent(double* r0, double* r1, double* z0, double* z1) 
    const = 0;

  virtual bool get_surface(const double r0, const double phi0, const double z0,
			   double ds, double** r, double** z, int* n);
  virtual double get_period() const
  { return 2.*M_PI; }
};

typedef std::deque<trace_field_source*> trace_source_list;

#endif
