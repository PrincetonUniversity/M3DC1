#ifndef M3DC1_FIELD_H
#define M3DC1_FIELD_H

#include <m3dc1_field.h>
#include <m3dc1_file.h>

#include "fusion_io_field.h"
#include "m3dc1_source.h"
#include "options.h"

class m3dc1_fio_field : public fio_field {
 protected:
  int time;
  double linfac;
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

// ELECTRIC FIELD
class m3dc1_electric_field : public m3dc1_fio_field {
 public:
  m3dc1_electric_field(m3dc1_source* s) 
    : m3dc1_fio_field(s) { }
  m3dc1_electric_field* clone() const 
  { return new m3dc1_electric_field(*this); }
  int load(const fio_option_list*);
  int dimension() const { return 1; }
  virtual int eval(const double*, double*);
};


// A
class m3dc1_vector_potential : public m3dc1_fio_field {
  m3dc1_field *psi0, *psi1, *psix;
  m3dc1_field *i0, *i1, *ix;
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
