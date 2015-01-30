#include "m3dc1_source.h"
#include <iostream>

m3dc1_source::m3dc1_source()
  : psi(0), f(0), g(0), psi_x(0), f_x(0), g_x(0)
{
  filename = "C1.h5";
  time = -1;
  factor = 1.;
  shift = 0.;
  use_g = true;
}

m3dc1_source::m3dc1_source(std::string f, int t)
  : psi(0), f(0), g(0), psi_x(0), f_x(0), g_x(0)
{
  filename = f;
  time = t;
  factor = 1.;
  shift = 0.;
  use_g = true;
}

m3dc1_source::~m3dc1_source()
{
}

bool m3dc1_source::load()
{
  int icomplex, i3d;

  if(!file.open(filename.data()))
    return false;

  file.read_parameter("bzero", &bzero);
  file.read_parameter("rzero", &rzero);
  file.read_parameter("extsubtract", &extsubtract);
  file.read_parameter("eqsubtract", &eqsubtract);
  file.read_parameter("icomplex", &icomplex);
  file.read_parameter("3d", &i3d);
  file.read_parameter("version", &version);
  file.read_parameter("itor", &itor);

  use_f = ((icomplex==1) || (i3d==1));
  toroidal = (itor == 1);
  period = toroidal ? 2.*M_PI : 2.*M_PI*rzero;

  if(time >= 0 && eqsubtract==1) bzero = 0.;

  std::cerr << "bzero = " << bzero << std::endl;
  std::cerr << "rzero = " << rzero << std::endl;
  std::cerr << "extsubtract = " << extsubtract << std::endl;
  std::cerr << "eqsubtract = " << eqsubtract << std::endl;
  std::cerr << "icomplex = " << icomplex << std::endl;
  std::cerr << "i3d = " << i3d << std::endl;

  m3dc1_scalar_list* xmag = file.read_scalar("xmag");
  m3dc1_scalar_list* zmag = file.read_scalar("zmag");
  m3dc1_scalar_list* psi0 = file.read_scalar("psimin");
  m3dc1_scalar_list* psi1 = file.read_scalar("psi_lcfs");
  if(!xmag || !zmag || !psi0 || !psi1)
    return false;

  std::cerr << "reading fields at time " << time << std::endl;
  psi = file.load_field("psi", time);
  if(!psi) return false;

  g = file.load_field("I", time);
  if(!g) {
    std::cerr << 
      "WARNING: Field 'I' not found.  Using bzero for toroidal field"
	      << std::endl;
    use_g = false;
    use_f = false;
  }

  if(use_f) {
    f = file.load_field("f", time);
    if(!f) return false;
  }

  if(time >= 0 && extsubtract==1  && version<8) {
    std::cerr << "reading external fields" << std::endl;
    psi_x = file.load_field("psi_ext", time);
    if(!psi_x) return false;

    if(use_g) {
      g_x = file.load_field("I_ext", time);
      if(!g_x) return false;
    }

    if(use_f) {
      f_x = file.load_field("f_ext", time);
      if(!f_x) return false;
    }
  } else extsubtract = 0;


  int index;
  if(eqsubtract) {
    index = 0;
  } else {
    index = -1;
    m3dc1_scalar_list* time = file.read_scalar("time");
    for(int i=0; i<time->size(); i++) 
      if(psi->time == time->at(i)) {
	index = i;
	break;
      }
    if(index==-1) {
      std::cerr << "ERROR: can't find scalar index for time = " << psi->time
		<< std::endl;
      return false;
    }
  }

  R_axis = xmag->at(index);
  Z_axis = zmag->at(index);
  psi_axis = psi0->at(index);
  psi_lcfs = psi1->at(index);

  std::cerr << "Magnetic axis = ( " << R_axis << ", " << Z_axis << " )"
	    << std::endl;
  std::cerr << "Psi at axis, lcfs = " << psi_axis << ", " << psi_lcfs
	    << std::endl;

  
  if(!file.close())
    return false;

  return true;
}

bool m3dc1_source::psibound(double* psi0, double* psi1) const
{
  *psi0 = psi_axis;
  *psi1 = psi_lcfs;
  return true;
}

bool m3dc1_source::eval_psi(const double r, const double z, double* p)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double val[m3dc1_field::OP_NUM];

  double phi = 0.;

  if(!psi->eval(r, phi, z, psiget, val))
    return false;

  *p = factor*val[m3dc1_field::OP_1];
  return true;
}

bool m3dc1_source::eval(const double r, const double phi0, const double z,
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

  double phi = phi0 - shift;
  double rr = (itor==1) ? r : 1.;
  double r0 = (itor==0) ? rzero : 1.;

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  if(!psi->eval(r, phi, z, psiget, val))
    return false;

  *b_r -= factor*val[m3dc1_field::OP_DZ]/rr;
  *b_z += factor*val[m3dc1_field::OP_DR]/rr;

  if(use_g) {
    if(!g->eval(r, phi, z, gget, val))
      return false;
    *b_phi += factor*val[m3dc1_field::OP_1]/rr;
  } else {
    *b_phi += factor*bzero;
  }

  if(use_f) {
    if(!f->eval(r, phi, z, fget, val))
      return false;

    *b_r -= factor*val[m3dc1_field::OP_DRP]/r0;
    *b_z -= factor*val[m3dc1_field::OP_DZP]/r0;
  }

  if(extsubtract==1) {
    if(!psi_x->eval(r, phi, z, psiget, val))
      return false;

    *b_r -= factor*val[m3dc1_field::OP_DZ]/rr;
    *b_z += factor*val[m3dc1_field::OP_DR]/rr;

    if(use_g) {
      if(!g_x->eval(r, phi, z, gget, val))
	return false;
      *b_phi += factor*val[m3dc1_field::OP_1]/rr;
    }

    if(use_f) {
      if(!f_x->eval(r, phi, z, fget, val))
	return false;

      *b_r -= factor*val[m3dc1_field::OP_DRP]/r0;
      *b_z -= factor*val[m3dc1_field::OP_DZP]/r0;
    }
  }

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
