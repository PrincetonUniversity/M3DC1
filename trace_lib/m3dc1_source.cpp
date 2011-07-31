#include "m3dc1_source.h"
#include <iostream>

m3dc1_source::m3dc1_source()
  : psi(0), f(0), g(0)
{
  filename = "C1.h5";
  time = -1;
  factor = 1.;
}

m3dc1_source::m3dc1_source(std::string f, int t)
  : psi(0), f(0), g(0)
{
  filename = f;
  time = t;
  factor = 1.;
}

m3dc1_source::~m3dc1_source()
{
}

bool m3dc1_source::load()
{
  if(!file.open(filename.data()))
    return false;

  file.read_parameter("bzero", &bzero);
  file.read_parameter("rzero", &rzero);

  std::cerr << "bzero = " << bzero << std::endl;
  std::cerr << "rzero = " << rzero << std::endl;

  m3dc1_scalar_list* xmag = file.read_scalar("xmag");
  m3dc1_scalar_list* zmag = file.read_scalar("zmag");
  if(!xmag || !zmag)
    return false;

  R_axis = xmag->at(0);
  Z_axis = zmag->at(0);

  std::cerr << "reading fields" << std::endl;
  psi = file.read_field("psi", time);
  g = file.read_field("I", time);
  f = file.read_field("f", time);

  if(!psi || !g || !f)
    return false;

  if(!file.close())
    return false;

  return true;
}

bool m3dc1_source::eval(const double r, const double phi, const double z,
			double* b_r, double* b_phi, double* b_z)
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

  if(!psi->eval(r, phi, z, psiget, val))
    return false;

  *b_r -= factor*val[m3dc1_field::OP_DZ]/r;
  *b_z += factor*val[m3dc1_field::OP_DR]/r;

  if(!g->eval(r, phi, z, gget, val))
    return false;
  *b_phi += factor*val[m3dc1_field::OP_1]/r;

  if(!f->eval(r, phi, z, fget, val))
    return false;

  *b_r -= factor*val[m3dc1_field::OP_DRP];
  *b_z -= factor*val[m3dc1_field::OP_DZP];

  return true;
}

bool m3dc1_source::center(double* r0, double* z0) const {
  *r0 = R_axis;
  *z0 = Z_axis;
  return true;
}

bool m3dc1_source::extent(double* r0, double* r1, double* z0, double* z1) const
{
  double phi0, phi1;

  psi->mesh->extent(r0, r1, &phi0, &phi1, z0, z1);

  return false;
}
