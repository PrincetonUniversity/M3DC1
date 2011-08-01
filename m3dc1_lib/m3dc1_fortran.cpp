#include "m3dc1_file.h"

#include <iostream>
#include <deque>

static m3dc1_file file;
static m3dc1_field *psi, *g, *f, *psi0, *g0;
static int eqsubtract;

struct field_data {
  std::string name;
  int time;
  m3dc1_field *field, *field0;
};

static typedef std::deque<field_data> handle_list;
static handle_list handles;


extern "C" void m3dc1_read_parameter_int_(const char* name, int* p, int* ierr)
{
  *ierr = 0;

  if(!file.read_parameter(name, p))
    *ierr = 1;
}

extern "C" void m3dc1_read_parameter_double_(const char* name, 
					     double* p, int* ierr)
{
  *ierr = 0;

  if(!file.read_parameter(name, p))
    *ierr = 1;
}

extern "C" void m3dc1_open_file_(const char* filename, int* ierr)
{
  *ierr = 0;

  if(!file.open(filename)) {
    *ierr = 1;
    return;
  }

  // determine if eqsubtract==1
  if(!file.read_parameter("eqsubtract", &eqsubtract)) {
    *ierr = 2;
    return;
  }

}

extern "C" void m3dc1_close_file_()
{
  file.close();
}

extern "C" void m3dc1_load_field_(const char* n, int* time, int* h, int* ierr)
{
  field_data fd;

  *ierr = 0;
  fd.field = file.load_field(n, *time);
  if(!fd.field) {
    *ierr = 1;
    return;
  }

  if(eqsubtract==1) {
    fd.field0 = file.load_field(n, -1);
    if(!fd.field0) {
      *ierr = 2;
      return;
    }
  }
 
  *h = handles.size();
  fd.name = n;
  fd.time = *time;
  handles.push_back(fd);
}

extern "C" void m3dc1_unload_field_(int* h, int* ierr)
{
  *ierr = 0;
  
  handle_list::const_reference i = handles.at(*h);

  if(!file.unload_field(i.name.c_str(), i.time)) {
    *ierr = 1;
  }
}

extern "C" void m3dc1_eval_field_(const int* h, 
				  const double* r, 
				  const double* phi,
				  const double* z,
				  double* v, 
				  int* ierr)
{
  const m3dc1_field::m3dc1_get_op op = 
    (m3dc1_field::m3dc1_get_op)(m3dc1_field::GET_VAL);

  *ierr = 0;

  handle_list::const_reference i = handles.at(*h);

  double val[m3dc1_field::OP_NUM];

  if(!i.field->eval(*r, *phi, *z, op, val)) {
    *ierr = 1;
    return;
  }

  *v = val[m3dc1_field::OP_1];

  if(eqsubtract==1) {
    if(!i.field0->eval(*r, *phi, *z, op, val)) {
      *ierr = 2;
      return;
    }

    *v += val[m3dc1_field::OP_1];
  }
}

extern "C" void m3dc1_load_magnetic_field_(int* time, int* ierr)
{
  *ierr = 0;

  // read fields at appropriate time
  psi = g = f = 0;
  psi = file.load_field("psi", *time);
  g = file.load_field("I", *time);
  f = file.load_field("f", *time);
  if(!psi || !g || !f) {
    *ierr = 1;
    return;
  }

  // if so, load equilibrium fields also
  if(eqsubtract==1) {
    psi0 = g0 = 0;
    psi0 = file.load_field("psi", -1);
    g0 = file.load_field("I", -1);
    if(!psi0 || !g0) {
      *ierr = 2;
      return;
    }
  }
}

extern "C" void m3dc1_unload_magnetic_field_(int* time, int* ierr)
{
  *ierr = 0;

  // unload fields
  if(!file.unload_field("psi", *time)) *ierr=1;
  if(!file.unload_field("I", *time)) *ierr=1;
  if(!file.unload_field("f", *time)) *ierr=1;

  // if so, unloadload equilibrium fields also
  if(eqsubtract==1) {
    if(!file.unload_field("psi", -1)) *ierr=1;
    if(!file.unload_field("I", -1)) *ierr=1;
  }
}

extern "C" void m3dc1_eval_magnetic_field_(const double* r,
					   const double* phi,
					   const double* z,
					   double* b_r, 
					   double* b_phi, 
					   double* b_z, 
					   int* ierr)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op gget = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op fget = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  double val[m3dc1_field::OP_NUM];

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  *b_r = 0;
  *b_z = 0;
  *b_phi = 0;

  if(!psi->eval(*r, *phi, *z, psiget, val)) {
    *ierr = 1;
    return;
  }
  *b_r -= val[m3dc1_field::OP_DZ] / *r;
  *b_z += val[m3dc1_field::OP_DR] / *r;

  if(!g->eval(*r, *phi, *z, gget, val)) {
    *ierr = 2;
    return;
  }
  *b_phi = val[m3dc1_field::OP_1] / *r;

  if(!f->eval(*r, *phi, *z, fget, val)) {
    *ierr = 3;
    return;
  }
  *b_r -= val[m3dc1_field::OP_DRP];
  *b_z -= val[m3dc1_field::OP_DZP];

  if(eqsubtract==1) {
    if(!psi0->eval(*r, *phi, *z, psiget, val)) {
      *ierr = 4;
      return;
    }
    *b_r -= val[m3dc1_field::OP_DZ] / *r;
    *b_z += val[m3dc1_field::OP_DR] / *r;

    if(!g0->eval(*r, *phi, *z, gget, val)) {
      *ierr = 5;
      return;
    }
    *b_phi = val[m3dc1_field::OP_1] / *r;
  }

  *ierr = 0;
}
