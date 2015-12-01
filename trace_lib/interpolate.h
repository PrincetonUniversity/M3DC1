bool bicubic_interpolation_coeffs(const double** x, const int m, const int n,
				  const int i, const int j, double a[4][4]);
bool cubic_interpolation_coeffs(const double* x, const int m, const int i,
				double* a);
bool cubic_interpolation(const int m, const double* p, const double p0, 
			 const double* f, double* f0);
