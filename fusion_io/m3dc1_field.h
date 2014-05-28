#ifndef M3DC1_FIELD_H
#define M3DC1_FIELD_H

#include <m3dc1_field.h>
#include <m3dc1_file.h>

#include "fusion_io_field.h"
#include "m3dc1_source.h"
#include "options.h"

class m3dc1_fio_series : public fio_series {
  m3dc1_scalar_list* data;
  m3dc1_scalar_list* time;
  std::string name;
  m3dc1_source* source;
  double factor;

 public:
  m3dc1_fio_series(m3dc1_source *s, const char* n, const double f)
    { source = s; name = n; factor = f; }
  int load();
  int eval(const double, double*);
};



class m3dc1_fio_field : public fio_field {
 protected:
  int time;
  double linfac, phase;
  bool eqsub, extsub, use_f;
  m3dc1_source* source;
 public:
  m3dc1_fio_field(m3dc1_source* s) 
    { source = s; }
  virtual int load(const fio_option_list*) = 0;
};


// Scalar (explicitly named field)
class m3dc1_scalar_field : public m3dc1_fio_field {
  m3dc1_field *f0, *f1, *fx;
  double factor;
  std::string name;

 public:
 m3dc1_scalar_field(m3dc1_source* s, const char* n, const double f) 
   : m3dc1_fio_field(s) { name = n; factor = f; }
  m3dc1_scalar_field* clone() const { return new m3dc1_scalar_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 1; }
  int eval(const double*, double*); 
};

// ALPHA
class m3dc1_alpha_field : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *i0, *i1, *ix;
  m3dc1_field *f0, *f1, *fx;
 public:
  m3dc1_alpha_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_alpha_field* clone() const { return new m3dc1_alpha_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 1; }
  virtual int eval(const double*, double*);
};

// SCALAR POTENTIAL
class m3dc1_phi_field : public m3dc1_fio_field {
 public:
  m3dc1_phi_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_phi_field* clone() const { return new m3dc1_phi_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 1; }
  virtual int eval(const double*, double*);
};

// ION_PRESSURE
class m3dc1_pi_field : public m3dc1_fio_field {
  m3dc1_scalar_field *p, *pe;
 public:
  m3dc1_pi_field(m3dc1_source* s)
    : m3dc1_fio_field(s) { }
  virtual ~m3dc1_pi_field();
  m3dc1_pi_field* clone() const { return new m3dc1_pi_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 1; }
  virtual int eval(const double*, double*);
};


// ELECTRIC FIELD
class m3dc1_electric_field : public m3dc1_fio_field {
  m3dc1_field *E0[3];
  m3dc1_field *E1[3];
 public:
  m3dc1_electric_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_electric_field* clone() const 
  { return new m3dc1_electric_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  virtual int eval(const double*, double*);
};

// V
class m3dc1_fluid_velocity : public m3dc1_fio_field {
  m3dc1_field *phi0, *phi1;
  m3dc1_field *w0, *w1;
  m3dc1_field *chi0, *chi1;
 public:
  m3dc1_fluid_velocity(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_fluid_velocity* clone() const 
  { return new m3dc1_fluid_velocity(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*);
};


// A
class m3dc1_vector_potential : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *f0, *f1, *fx;
 public:
  m3dc1_vector_potential(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_vector_potential* clone() const 
  { return new m3dc1_vector_potential(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*);
};


// B
class m3dc1_magnetic_field : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *i0, *i1, *ix;
  m3dc1_field *f0, *f1, *fx;
 public:
  m3dc1_magnetic_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_magnetic_field* clone() const 
  { return new m3dc1_magnetic_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*);
  int eval_deriv(const double*, double*);
};

// J
class m3dc1_current_density : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *i0, *i1, *ix;
  m3dc1_field *f0, *f1, *fx;
 public:
  m3dc1_current_density(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_current_density* clone() const 
  { return new m3dc1_current_density(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*);
};


// V
class m3dc1_velocity_field : public m3dc1_fio_field {
  m3dc1_field *u0, *u1;
  m3dc1_field *v0, *v1;
  m3dc1_field *x0, *x1;
 public:
  m3dc1_velocity_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_velocity_field* clone() const
  { return new m3dc1_velocity_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 3; }
  int eval(const double*, double*);
};

// Grad(A)
class m3dc1_grad_vector_potential : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *i0, *i1, *ix;
  m3dc1_field *f0, *f1, *fx;
 public:
  m3dc1_grad_vector_potential(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_grad_vector_potential* clone() const 
  { return new m3dc1_grad_vector_potential(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 9; }
  int eval(const double*, double*);
};


#endif
