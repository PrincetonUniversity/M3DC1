#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "fusion_io.h"
#include "options.h"

#include <iostream>

int main()
{
  int result;
  fio_source* src = new m3dc1_source();
  fio_field *pressure, *density, *magnetic_field;
  fio_option_list opt;

  result = src->open("C1.h5");
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    delete(src);
    return result;
  };

  src->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, 1);
  opt.set_option(FIO_PERTURBED_ONLY, 1);
  opt.set_option(FIO_LINEAR_SCALE, 10.);

  result = src->get_field(FIO_PRESSURE, &pressure, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening pressure field" << std::endl;
    delete(src);
    return result;
  };

  result = src->get_field(FIO_DENSITY, &density, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening density field" << std::endl;
    delete(src);
    return result;
  };

  result = src->get_field(FIO_MAGNETIC_FIELD, &magnetic_field, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field" << std::endl;
    delete(src);
    return result;
  };



  int npts = 10;
  double R0 = 1.6;
  double R1 = 2.1;
  double Z0 = 0.0;
  double Z1 = 0.0;
  double phi0 = 0.;
  double phi1 = 0.;
  double x[3];
  double p, n, b[3];
  
  for(int i=0; i<npts; i++) {
    x[0] = R0 + (R1-R0)*i/(npts-1);
    x[1] = phi0 + (phi1-phi0)*i/(npts-1);
    x[2] = Z0 + (Z1-Z0)*i/(npts-1);
    
    std::cout << "(" << x[0] << ", " << x[1] << ", " << x[2] << "):\n";

    result = pressure->eval(x, &p);
    std::cout << "\tpressure = " << p << "\n";

    result = density->eval(x, &n);
    std::cout << "\tdensity = " << n << "\n";

    result = magnetic_field->eval(x, b);
    std::cout << "\tB = " "(" << b[0] << ", " << b[1] << ", " << b[2] << "):\n";
  }

  src->close();
  delete(pressure);
  delete(density);
  delete(magnetic_field);
  delete(src);

  return 0;
}


