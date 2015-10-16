#ifndef GEQDSK_H
#define GEQDSK_h

#include "trace_source.h"

#include <string>

class geqdsk_source : public trace_field_source {
  double rmaxis, zmaxis;

  int nw, nh;
  double dx, dz;
  double rleft, zmid;

  double* psi;
  double** psirz;
  double* fpol;
  

 public:
  std::string filename;
  
  geqdsk_source();
  ~geqdsk_source();

  bool load();
  bool eval(const double r, const double phi, const double z,
	    double *b_r, double *b_phi, double *b_z);
  bool center(double* r0, double* z0) const;
  bool extent(double* r0, double* r1, double* z0, double* z1) const;
};

#endif
