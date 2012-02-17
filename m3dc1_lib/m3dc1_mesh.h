#ifndef M3DC1_MESH_H
#define M3DC1_MESH_H

#include <math.h>

#define TOL 1e-3

class m3dc1_mesh {
 private:
  int memory_depth;

  int hits, misses;
  
  int last_elm; 
  int** next_elm;

  void clear_memory();

 protected:
  //  static const double tol = 1e-1;

  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    double t = (a[i] + b[i] + c[i])*TOL;
    if(eta + t < 0.) return false;
    if(eta - t > c[i]) return false;
    double x = 1. - eta/c[i];
    if(xi + t < -b[i]*x) return false;
    if(xi - t > a[i]*x) return false;
    return true;
  }
  virtual void global_to_local(const int i, 
			       const double X,const double Phi,const double Z,
			       double* xi, double* zi, double* eta) const
  {
    *xi  =  (X - x[i])*co[i] + (Z - z[i])*sn[i] - b[i];
    *eta = -(X - x[i])*sn[i] + (Z - z[i])*co[i];
  }

 public:
  int nelms;
  double* a;
  double* b;
  double* c;
  double* co;
  double* sn;
  double* x;
  double* z;
  int* bound;

 public:
  m3dc1_mesh(int n);
  virtual ~m3dc1_mesh(); 

  bool set_memory_depth(int d);

  bool is_in_element(const int i, 
		     const double X, const double Phi, const double Z,
		     double* xi_out=0, double* zi_out=0, double* eta_out=0)
    const
  {
    double xi, zi, eta;
    global_to_local(i, X, Phi, Z, &xi, &zi, &eta);
    if(is_in_element_local(i, xi, zi, eta)) {
      if(xi_out) *xi_out = xi;
      if(zi_out) *zi_out = zi;
      if(eta_out) *eta_out = eta;
      return true;
    }
    return false;
  }
  virtual int in_element(double X, double Phi, double Z, 
			 double* xi=0, double* zi=0, double* eta=0, 
			 int guess=-1);

  virtual void extent(double *X0, double* X1,
		      double *Phi0, double* Phi1,
		      double *Z0, double* Z1) const;

};



class m3dc1_3d_mesh : public m3dc1_mesh {
 protected:
  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    if(!m3dc1_mesh::is_in_element_local(i, xi, zi, eta)) return false;
    if(zi + TOL*d[i] < 0.) return false;
    if(zi - TOL*d[i] > d[i]) return false;
    return true;
  }
  virtual void global_to_local(const int i, 
			       const double X,const double Phi,const double Z,
			       double* xi, double* zi, double* eta) const
  {
    m3dc1_mesh::global_to_local(i, X, Phi, Z, xi, zi, eta);
    *zi = Phi - phi[i];
  }

 public:
  double *phi;
  double *d;

  virtual int in_element(double X, double Phi, double Z, 
			 double* xi=0, double* zi=0, double* eta=0, 
			 int guess=-1)
  {
    while(Phi < 0       ) Phi += 2.*M_PI;
    while(Phi >= 2.*M_PI) Phi -= 2.*M_PI;

    return m3dc1_mesh::in_element(X, Phi, Z, xi, zi, eta, guess);
  }


  m3dc1_3d_mesh(int n);
  virtual ~m3dc1_3d_mesh();
};


#endif
