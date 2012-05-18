#include "m3dc1_source.h"
#include "m3dc1_field.h"
#include "fusion_io.h"
#include "options.h"
#include "compound_field.h"

#include <iostream>

int main()
{
  int result;
  fio_source* src = new m3dc1_source();
  fio_field *pressure, *density, *magnetic_field;
  fio_option_list opt;

  // Open an m3dc1 source
  result = src->open("C1.h5");
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    delete(src);
    return result;
  };

  // set options for fields obtained from this source
  src->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, 1);
  opt.set_option(FIO_PERTURBED_ONLY, 1);
  opt.set_option(FIO_LINEAR_SCALE, 10.);

  // open fields
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

  
  // Open a second source
  fio_source* src2 = new m3dc1_source();
  fio_field* magnetic_field2;
  result = src2->open("C2.h5");
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    delete(src2);
    return result;
  };
  // set options
  src2->get_field_options(&opt);
  opt.set_option(FIO_TIMESLICE, 0);
  opt.set_option(FIO_PERTURBED_ONLY, 1);
  opt.set_option(FIO_LINEAR_SCALE, 10.);

  // open field
  result = src2->get_field(FIO_MAGNETIC_FIELD, &magnetic_field2, &opt);
  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening magnetic field" << std::endl;
    delete(src2);
    return result;
  };

  // create compound field
  fio_compound_field total_field;
  total_field.add_field(magnetic_field,  FIO_ADD);
  total_field.add_field(magnetic_field2, FIO_ADD, -1.);

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
    std::cout << "\tB = (" << b[0] << ", " << b[1] << ", " << b[2] << "):\n";

    result = magnetic_field2->eval(x, b);
    std::cout << "\tB = (" << b[0] << ", " << b[1] << ", " << b[2] << "):\n";

    result = total_field.eval(x, b);
    std::cout << "\tB = (" << b[0] << ", " << b[1] << ", " << b[2] << "):\n";
  }

  src->close();
  delete(pressure);
  delete(density);
  delete(magnetic_field);
  delete(magnetic_field2);
  delete(src);
  delete(src2);

  return 0;
}


