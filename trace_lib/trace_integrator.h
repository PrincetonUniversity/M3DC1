#ifndef TRACE_INTEGRATOR_H
#define TRACE_INTEGRATOR_H

#include <fstream>

#include "trace_source.h"
#include "interpolation_source.h"

class trace_integrator : private trace_field_source {
  double R, Phi, Z;

  double dr[4], dz[4];

  std::fstream file;

 public:

  struct integrator_data {
    int toroidal_transits;
    int poloidal_transits;
    double q;
    double distance;
  };

  double plane;
  trace_source_list sources;
  interpolation_source interp_source;

  trace_integrator();
  virtual ~trace_integrator();

  virtual bool load();
  virtual bool eval(const double r, const double phi, const double z,
		    double* b_r, double* b_phi, double* b_z);
  virtual bool center(double* r0, double* z0) const;
  virtual bool extent(double* r0, double* r1, double* z0, double* z1) const;
  virtual bool get_surface(const double r0, const double phi0, const double z0,
			   const double ds, double** r, double** z, int* n);


  double find_min_bn(const double, const double);

  bool set_pos(const double r, const double phi, const double z);
  void get_pos(double *r, double* phi, double* z) const
  { *r = R; *phi = Phi; *z = Z; }
  bool step_euler(double dphi);
  bool step_rk4(double dphi);
  bool step_predcorr(double dphi);
  bool integrate(int transits, int steps_per_transit, 
		 integrator_data* data=0);

  void plot_point();

  bool open_file(const char*);
  bool close_file();
};

#endif
