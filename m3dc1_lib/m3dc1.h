#ifndef M3DC1_H
#define M3DC1_H

#include <hdf5.h>
#include <map>
#include <string>
#include <math.h>


class m3dc1_mesh {
 protected:
  static const double tol = 1e-1;

  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    double t = (a[i] + b[i] + c[i])*tol;
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
		 int* guess=0, const int num_guesses=0) const
  {
    for(int i=0; i<num_guesses; i++)
      if(is_in_element(guess[i],X,Phi,Z,xi,zi,eta))
	return guess[i];

    for(int i=0; i<nelms; i++)
      if(is_in_element(i,X,Phi,Z,xi,zi,eta))
	return i;

    return -1;
  }
  virtual void extent(double *X0, double* X1,
		      double *Phi0, double* Phi1,
		      double *Z0, double* Z1) const;

 public:
  m3dc1_mesh(int n);
  virtual ~m3dc1_mesh(); 
};

class m3dc1_3d_mesh : public m3dc1_mesh {
 protected:
  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    if(!m3dc1_mesh::is_in_element_local(i, xi, zi, eta)) return false;
    if(zi + tol*d[i] < 0.) return false;
    if(zi - tol*d[i] > d[i]) return false;
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
			 int* guess=0, const int num_guesses=0) const
  {
    while(Phi < 0       ) Phi += 2.*M_PI;
    while(Phi >= 2.*M_PI) Phi -= 2.*M_PI;

    return m3dc1_mesh::in_element(X, Phi, Z, xi, zi, eta, guess, num_guesses);
  }


  m3dc1_3d_mesh(int n);
  virtual ~m3dc1_3d_mesh();
};


class m3dc1_field {
 public:
  double* data;

  static const int nbasis = 20;
  static const int mi[nbasis];
  static const int ni[nbasis];

  enum m3dc1_get_op {
    GET_VAL    = 1,
    GET_DVAL   = 2,
    GET_DDVAL  = 4,
    GET_PVAL   = 8,
    GET_PPVAL  = 16
  };

  enum m3dc1_op {
    OP_1     = 0,
    OP_DR    = 1,
    OP_DZ    = 2,
    OP_DRR   = 3,
    OP_DRZ   = 4,
    OP_DZZ   = 5,
    OP_DP    = 6,
    OP_DRP   = 7,
    OP_DZP   = 8,
    OP_DRRP  = 9,
    OP_DRZP  = 10,
    OP_DZZP  = 11,
    OP_DPP   = 12,
    OP_DRPP  = 13,
    OP_DZPP  = 14,
    OP_DRRPP = 15,
    OP_DRZPP = 16,
    OP_DZZPP = 17,
    OP_NUM   = 18
  };

 public:
  double time;
  m3dc1_mesh* mesh;

  m3dc1_field() { mesh=0; data=0; }
  m3dc1_field(m3dc1_mesh* m);
  virtual ~m3dc1_field();

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element, int* guess, const int num_guesses);
};

class m3dc1_complex_field : public m3dc1_field {
 public:
  double* data_i;
  int ntor;

 public:
  m3dc1_complex_field(m3dc1_mesh* m, int ntor);
  virtual ~m3dc1_complex_field();

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element, int* guess, const int num_guesses);
};

class m3dc1_3d_field : public m3dc1_field {
 public:
  static const int tbasis = 4;
  static const int li[tbasis];

 public:
  m3dc1_3d_field(m3dc1_mesh* m);

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element, int* guess, const int num_guesses);
};

class m3dc1_timeslice {
 public:
  typedef std::map<std::string, m3dc1_field*> m3dc1_field_map;
  m3dc1_field_map field_map;

  int ntor;
  int is_3d;

  double time;
  m3dc1_mesh* mesh;

  m3dc1_timeslice();
  ~m3dc1_timeslice();

  bool get_field(const char*, m3dc1_field**) const;
};

#endif
