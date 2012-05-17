#ifndef GEQDSK_SOURCE_H
#define GEQDSK_SOURCE_H

#include "source.h"

#include <string>

class geqdsk_source : public field_source {
  double rmaxis, zmaxis;

  int nw, nh;
  double dx, dz;
  double rleft, zmid;

  double* psi;
  double* psirz;
  double* fpol;
  

 public:
  std::string filename;
  
  geqdsk_source();
  ~geqdsk_source();

  int load();
  bool eval(const double r, const double phi, const double z,
	    double *b_r, double *b_phi, double *b_z);
  bool center(double* r0, double* z0) const;
  bool extent(double* r0, double* r1, double* z0, double* z1) const;

  virtual int eval_scalar_field(const scalar_field, const double*, double*);
  virtual int eval_vector_field(const vector_field, const double*, double*);
  virtual int eval_scalar_field_grad(const scalar_field,const double*,double*);
};

#endif
