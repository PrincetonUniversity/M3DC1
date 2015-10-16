#ifndef TRACE_M3DC1_SOURCE_H
#define TRACE_M3DC1_SOURCE_H

#include "trace_source.h"
#include "m3dc1_file.h"

#include <string>

class m3dc1_source : public trace_field_source {
  m3dc1_file file;
  m3dc1_field *psi, *f, *g;
  m3dc1_field *psi_x, *f_x, *g_x;

  bool use_f, use_g;
  int extsubtract, eqsubtract, version, itor;
  double bzero, rzero;
  double R_axis, Z_axis;
  double psi_axis, psi_lcfs;
  double period;

 public:
  std::string filename;
  int time;
  double factor;
  double shift;

  m3dc1_source();
  m3dc1_source(std::string f, int t);
  ~m3dc1_source();

  bool load();
  bool eval(const double r, const double phi, const double z,
	    double* b_r, double* b_phi, double* b_z);
  bool eval_psi(const double r, const double z, double* p);
  bool psibound(double* psi0, double* psi1) const;
  virtual bool center(double* r0, double* z0) const;
  virtual bool extent(double* r0, double* r1, double* z0, double* z1) const;
  virtual double get_period() const
  { return period; }
};

#endif
