#ifndef COIL_SOURCE_H
#define COIL_SOURCE_H

#include <list>

#include "source.h"

struct coil_segment {
private:
  static void biot_savart_integrand(const double, const double, const double,
				    const double, const double, const double,
				    const double, const double, const double,
				    double*, double*, double*);
public:
  double current;
  double R[2];
  double Phi[2];
  double Z[2];

  coil_segment()
  { }
  coil_segment(double i, 
	       double r0, double r1, 
	       double phi0, double phi1,
	       double z0, double z1)
  {
    current = i;
    R[0] = r0;     R[1] = r1;
    Phi[0] = phi0; Phi[1] = phi1;
    Z[0] = z0;     Z[1] = z1;
  }
  coil_segment(const coil_segment& c)
  {
    current = c.current;
    R[0] = c.R[0];     R[1] = c.R[1];
    Phi[0] = c.Phi[0]; Phi[1] = c.Phi[1];
    Z[0] = c.Z[0];     Z[1] = c.Z[1];
  }

  bool eval(const double, const double, const double,
	    double*, double*, double*);
};


class coil_source : public field_source, public std::list<coil_segment> {
  
 private:

 public:
  coil_source();
  virtual ~coil_source();

  int load();
  bool eval(const double, const double, const double,
	    double*, double*, double*);
  virtual bool center(double*, double*) const
  { return false;}
  virtual bool extent(double*, double*, double*, double*) const
  { return false;}
};

#endif
