#ifndef INTERPOLATE_H
#define INTERPOLATE_H

extern "C" {
  void bicubic_interpolation_(const int*, const int*, 
			      const double*, const double*, 
			      const double*, const double*, const double*, 
			      double*, double*, double*, int*);
  void bicubic_interpolation_coeffs_(const double*, const int*, const int*, 
				     const double*, const double*, 
				     double*, int*);
  void cubic_interpolation_(const int*, const double*, const double*, 
			    const double*, double*);
}

#endif
