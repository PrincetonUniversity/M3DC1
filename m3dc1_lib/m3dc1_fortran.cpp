#include "m3dc1_file.h"

#include <deque>
#include <iostream>

static m3dc1_file file;
static m3dc1_field *psi, *g, *f, *psi0, *g0, *psi_ext, *g_ext, *f_ext;
static int eqsubtract;
static int extsubtract;
static double scale_factor = 1;

struct field_data {
  std::string name;
  int time;
  m3dc1_field *field, *field0;
};

typedef std::deque<field_data> handle_list;
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

  // determine if extsubtract==1
  if(!file.read_parameter("extsubtract", &extsubtract)) {
    extsubtract = 0;
  }
}

extern "C" void m3dc1_close_file_()
{
  file.close();
}

extern "C" void m3dc1_extent_(int *time, 
			      double* r0, double* r1, 
			      double* phi0, double* phi1,
			      double* z0, double* z1, 
			      int* ierr)
{
  *ierr = 0;
  if(!file.extent(*time, r0, r1, phi0, phi1, z0, z1)) 
    *ierr = 1;
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

  int guess = -1;
  double val[m3dc1_field::OP_NUM];

  if(!i.field->eval(*r, *phi, *z, op, val, &guess)) {
    *ierr = 1;
    return;
  }

  *v = val[m3dc1_field::OP_1];

  if(eqsubtract==1) {
    if(!i.field0->eval(*r, *phi, *z, op, val, &guess)) {
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

  // read external fields
  if(extsubtract==1) {
    psi_ext = g_ext = f_ext = 0;
    psi_ext = file.load_field("psi_ext", *time);
    g_ext = file.load_field("I_ext", *time);
    f_ext = file.load_field("f_ext", *time);
    if(!psi_ext || !g_ext || !f_ext) {
      *ierr = 3;
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

  if(extsubtract==1) {
    if(!file.unload_field("psi_ext", *time)) *ierr=1;
    if(!file.unload_field("I_ext", *time)) *ierr=1;
    if(!file.unload_field("f_ext", *time)) *ierr=1;
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

  int guess = -1;

  *b_r = 0;
  *b_z = 0;
  *b_phi = 0;

  if(!psi->eval(*r, *phi, *z, psiget, val, &guess)) {
    *ierr = 1;
    return;
  }
  *b_r -= scale_factor*val[m3dc1_field::OP_DZ] / *r;
  *b_z += scale_factor*val[m3dc1_field::OP_DR] / *r;

  if(!g->eval(*r, *phi, *z, gget, val, &guess)) {
    *ierr = 2;
    return;
  }
  *b_phi += scale_factor*val[m3dc1_field::OP_1] / *r;

  if(!f->eval(*r, *phi, *z, fget, val, &guess)) {
    *ierr = 3;
    return;
  }
  *b_r -= scale_factor*val[m3dc1_field::OP_DRP];
  *b_z -= scale_factor*val[m3dc1_field::OP_DZP];

  if(eqsubtract==1) {
    if(!psi0->eval(*r, *phi, *z, psiget, val, &guess)) {
      *ierr = 4;
      return;
    }
    *b_r -= val[m3dc1_field::OP_DZ] / *r;
    *b_z += val[m3dc1_field::OP_DR] / *r;

    if(!g0->eval(*r, *phi, *z, gget, val, &guess)) {
      *ierr = 5;
      return;
    }
    *b_phi += val[m3dc1_field::OP_1] / *r;
  }

  if(extsubtract==1) {
    if(!psi_ext->eval(*r, *phi, *z, psiget, val, &guess)) {
      *ierr = 6;
      return;
    }
    *b_r -= scale_factor*val[m3dc1_field::OP_DZ] / *r;
    *b_z += scale_factor*val[m3dc1_field::OP_DR] / *r;
    
    if(!g_ext->eval(*r, *phi, *z, gget, val, &guess)) {
      *ierr = 7;
      return;
    }
    *b_phi += scale_factor*val[m3dc1_field::OP_1] / *r;
    
    if(!f_ext->eval(*r, *phi, *z, fget, val, &guess)) {
      *ierr = 8;
      return;
    }
    *b_r -= scale_factor*val[m3dc1_field::OP_DRP];
    *b_z -= scale_factor*val[m3dc1_field::OP_DZP];
  }

  *ierr = 0;
}

extern "C" void m3dc1_eval_alpha_(const double* r,
				  const double* phi,
				  const double* z,
				  double* alpha, 
				  int* ierr)
{
  const m3dc1_field::m3dc1_get_op getdval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op getval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double psi0_val[m3dc1_field::OP_NUM], g0_val[m3dc1_field::OP_NUM];
  double psi1_val[m3dc1_field::OP_NUM], f1_val[m3dc1_field::OP_NUM];
  double temp[m3dc1_field::OP_NUM];
  double b2, ab, r2;

  int i;
  int guess = -1;

  if(eqsubtract != 1) {
    std::cerr << "Calculation of alpha only supported for eqsubtract=1" 
	 << std::endl;
    *ierr = 1;
    return;
  }
    
  if(!psi0->eval(*r, *phi, *z, getdval, psi0_val, &guess)) {
    *ierr = 2;
    return;
  }
  if(!g0->eval(*r, *phi, *z, getval, g0_val, &guess)) {
    *ierr = 3;
    return;
  }
  if(!psi->eval(*r, *phi, *z, getval, psi1_val, &guess)) {
    *ierr = 4;
    return;
  }
  if(!f->eval(*r, *phi, *z, getdval, f1_val, &guess)) {
    *ierr = 5;
    return;
  }
  if(extsubtract==1) {
    if(!psi_ext->eval(*r, *phi, *z, getval, temp, &guess)) {
      *ierr = 6;
      return;
    }
    for(i=0; i<m3dc1_field::OP_NUM; i++) psi1_val[i] += temp[i];

    if(!f_ext->eval(*r, *phi, *z, getdval, temp, &guess)) {
      *ierr = 7;
      return;
    }
    for(i=0; i<m3dc1_field::OP_NUM; i++) f1_val[i] += temp[i];
  }

  // alpha = B0.A1 / B0.B0

  r2 = (*r)*(*r);
  ab = g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_1] / r2
    - psi0_val[m3dc1_field::OP_DR]*f1_val[m3dc1_field::OP_DR]
    - psi0_val[m3dc1_field::OP_DZ]*f1_val[m3dc1_field::OP_DZ];
  b2 = (psi0_val[m3dc1_field::OP_DR]*psi0_val[m3dc1_field::OP_DR] +
	psi0_val[m3dc1_field::OP_DZ]*psi0_val[m3dc1_field::OP_DZ] + 
	g0_val[m3dc1_field::OP_1]*g0_val[m3dc1_field::OP_1])/r2;
  
  *alpha = scale_factor*ab/b2;
    
  *ierr = 0;
}

extern "C" void m3dc1_eval_grad_alpha_(const double* r,
				       const double* phi,
				       const double* z,
				       double* dadR,
				       double* dadPhi,
				       double* dadZ,
				       int* ierr)
{
  const m3dc1_field::m3dc1_get_op getdpval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_DDVAL | m3dc1_field::GET_PVAL);

  const m3dc1_field::m3dc1_get_op getpval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  const m3dc1_field::m3dc1_get_op getdval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_DDVAL);

  const m3dc1_field::m3dc1_get_op getval = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL);

  double psi0_val[m3dc1_field::OP_NUM], g0_val[m3dc1_field::OP_NUM];
  double psi1_val[m3dc1_field::OP_NUM], f1_val[m3dc1_field::OP_NUM];
  double temp[m3dc1_field::OP_NUM];
  double b2, ab, r2, dabdR, dabdPhi, dabdZ, db2dR, db2dZ;

  int i;
  int guess = -1;

  if(eqsubtract != 1) {
    std::cerr << "Calculation of alpha only supported for eqsubtract=1" 
	 << std::endl;
    *ierr = 1;
    return;
  }
    
  if(!psi0->eval(*r, *phi, *z, getdval, psi0_val, &guess)) {
    *ierr = 2;
    return;
  }
  if(!g0->eval(*r, *phi, *z, getval, g0_val, &guess)) {
    *ierr = 3;
    return;
  }
  if(!psi->eval(*r, *phi, *z, getpval, psi1_val, &guess)) {
    *ierr = 4;
    return;
  }
  if(!f->eval(*r, *phi, *z, getdpval, f1_val, &guess)) {
    *ierr = 5;
    return;
  }
  if(extsubtract==1) {
    if(!psi_ext->eval(*r, *phi, *z, getpval, temp, &guess)) {
      *ierr = 6;
      return;
    }
    for(i=0; i<m3dc1_field::OP_NUM; i++) psi1_val[i] += temp[i];

    if(!f_ext->eval(*r, *phi, *z, getdpval, temp, &guess)) {
      *ierr = 7;
      return;
    }
    for(i=0; i<m3dc1_field::OP_NUM; i++) f1_val[i] += temp[i];
  }

  // alpha = B0.A1 / B0.B0

  r2 = (*r)*(*r);
  ab = g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_1] / r2
    - psi0_val[m3dc1_field::OP_DR]*f1_val[m3dc1_field::OP_DR]
    - psi0_val[m3dc1_field::OP_DZ]*f1_val[m3dc1_field::OP_DZ];
  b2 = (psi0_val[m3dc1_field::OP_DR]*psi0_val[m3dc1_field::OP_DR] +
	psi0_val[m3dc1_field::OP_DZ]*psi0_val[m3dc1_field::OP_DZ] + 
	g0_val[m3dc1_field::OP_1]*g0_val[m3dc1_field::OP_1])/r2;

  dabdR = g0_val[m3dc1_field::OP_DR]*psi1_val[m3dc1_field::OP_1] / r2
    - psi0_val[m3dc1_field::OP_DRR]*f1_val[m3dc1_field::OP_DR]
    - psi0_val[m3dc1_field::OP_DRZ]*f1_val[m3dc1_field::OP_DZ]
    + g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_DR] / r2
    - psi0_val[m3dc1_field::OP_DR]*f1_val[m3dc1_field::OP_DRR]
    - psi0_val[m3dc1_field::OP_DZ]*f1_val[m3dc1_field::OP_DRZ] -
    2.*g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_1] / ((*r)*r2);

  dabdPhi = g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_DP] / r2
    - psi0_val[m3dc1_field::OP_DR]*f1_val[m3dc1_field::OP_DRP]
    - psi0_val[m3dc1_field::OP_DZ]*f1_val[m3dc1_field::OP_DZP];

  dabdZ = g0_val[m3dc1_field::OP_DZ]*psi1_val[m3dc1_field::OP_1] / r2
    - psi0_val[m3dc1_field::OP_DRZ]*f1_val[m3dc1_field::OP_DR]
    - psi0_val[m3dc1_field::OP_DZZ]*f1_val[m3dc1_field::OP_DZ] 
    + g0_val[m3dc1_field::OP_1]*psi1_val[m3dc1_field::OP_DZ] / r2
    - psi0_val[m3dc1_field::OP_DR]*f1_val[m3dc1_field::OP_DRZ]
    - psi0_val[m3dc1_field::OP_DZ]*f1_val[m3dc1_field::OP_DZZ];

  db2dR = 2.*(psi0_val[m3dc1_field::OP_DRR]*psi0_val[m3dc1_field::OP_DR] +
	      psi0_val[m3dc1_field::OP_DRZ]*psi0_val[m3dc1_field::OP_DZ] + 
	      g0_val[m3dc1_field::OP_DR]*g0_val[m3dc1_field::OP_1])/r2  
    - 2.*b2/(*r);

  db2dZ = 2.*(psi0_val[m3dc1_field::OP_DRZ]*psi0_val[m3dc1_field::OP_DR] +
	      psi0_val[m3dc1_field::OP_DZZ]*psi0_val[m3dc1_field::OP_DZ] + 
	      g0_val[m3dc1_field::OP_DZ]*g0_val[m3dc1_field::OP_1])/r2;
  
  *dadR = scale_factor*(dabdR/b2 - ab*db2dR/(b2*b2));
  *dadPhi = scale_factor*(dabdPhi/b2);
  *dadZ = scale_factor*(dabdZ/b2 - ab*db2dZ/(b2*b2));
    
  *ierr = 0;
}


extern "C" void m3dc1_eval_equilibrium_magnetic_field_(const double* r,
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

  double val[m3dc1_field::OP_NUM];

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  int guess = -1;

  *b_r = 0;
  *b_z = 0;
  *b_phi = 0;

  if(eqsubtract==1) {
    if(!psi0->eval(*r, *phi, *z, psiget, val, &guess)) {
      *ierr = 4;
      return;
    }
    *b_r -= val[m3dc1_field::OP_DZ] / *r;
    *b_z += val[m3dc1_field::OP_DR] / *r;

    if(!g0->eval(*r, *phi, *z, gget, val, &guess)) {
      *ierr = 5;
      return;
    }
    *b_phi += val[m3dc1_field::OP_1] / *r;
  }

  *ierr = 0;
}


extern "C" void m3dc1_eval_perturbed_magnetic_field_(const double* r,
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

  int guess = -1;

  *b_r = 0;
  *b_z = 0;
  *b_phi = 0;

  if(!psi->eval(*r, *phi, *z, psiget, val, &guess)) {
    *ierr = 1;
    return;
  }
  *b_r -= scale_factor*val[m3dc1_field::OP_DZ] / *r;
  *b_z += scale_factor*val[m3dc1_field::OP_DR] / *r;

  if(!g->eval(*r, *phi, *z, gget, val, &guess)) {
    *ierr = 2;
    return;
  }
  *b_phi += scale_factor*val[m3dc1_field::OP_1] / *r;

  if(!f->eval(*r, *phi, *z, fget, val, &guess)) {
    *ierr = 3;
    return;
  }
  *b_r -= scale_factor*val[m3dc1_field::OP_DRP];
  *b_z -= scale_factor*val[m3dc1_field::OP_DZP];

  if(extsubtract==1) {
    if(!psi_ext->eval(*r, *phi, *z, psiget, val, &guess)) {
      *ierr = 4;
      return;
    }
    *b_r -= scale_factor*val[m3dc1_field::OP_DZ] / *r;
    *b_z += scale_factor*val[m3dc1_field::OP_DR] / *r;
    
    if(!g_ext->eval(*r, *phi, *z, gget, val, &guess)) {
    *ierr = 5;
    return;
    }
    *b_phi += scale_factor*val[m3dc1_field::OP_1] / *r;

    if(!f_ext->eval(*r, *phi, *z, fget, val, &guess)) {
    *ierr = 6;
    return;
    }
    *b_r -= scale_factor*val[m3dc1_field::OP_DRP];
    *b_z -= scale_factor*val[m3dc1_field::OP_DZP];
  }

  *ierr = 0;
}


// the following are provided for backward compatibility
extern "C" void m3dc1_load_file_(int* time,  int* ierr)
{
  m3dc1_open_file_("C1.h5", ierr);
  if(*ierr != 0) return;

  m3dc1_load_magnetic_field_(time, ierr);
}

extern "C" void m3dc1_unload_file_()
{
  m3dc1_close_file_();
}

extern "C" void m3dc1_get_field_(const double* R, const double* Phi, const double* Z, 
				 double* Br, double* Bphi, double* Bz)
{
  int ierr;
  m3dc1_eval_magnetic_field_(R, Phi, Z, Br, Bphi, Bz, &ierr);
}

extern "C" void m3dc1_get_field0_(const double* R, const double* Phi, const double* Z, 
				 double* Br, double* Bphi, double* Bz)
{
  int ierr;
  m3dc1_eval_equilibrium_magnetic_field_(R, Phi, Z, Br, Bphi, Bz, &ierr);
}

extern "C" void m3dc1_get_field1_(const double* R, const double* Phi, const double* Z, 
				 double* Br, double* Bphi, double* Bz)
{
  int ierr;
  m3dc1_eval_perturbed_magnetic_field_(R, Phi, Z, Br, Bphi, Bz, &ierr);
}


extern "C" void m3dc1_get_num_timesteps_(int* n, int* ierr)
{
  m3dc1_scalar_list* scalar_list = file.read_scalar("time");

  if(!scalar_list) {
    *ierr = 1;
    return;
  }

  *n = scalar_list->size();
  *ierr = 0;
}

extern "C" void m3dc1_read_scalar_(const char* name, double* scalar, 
				   const int* n, int* ierr)
{
  m3dc1_scalar_list* scalar_list = file.read_scalar(name);

  *ierr = 0;

  if(!scalar_list) {
    *ierr = 1;
    return;
  }

  int sz = scalar_list->size();
  if(sz > *n) sz = *n;

  for(int i=0; i<sz; i++)
    scalar[i] = scalar_list->at(i);
}

extern "C" void m3dc1_set_scale_factor_(const double* d)
{
  scale_factor = *d;
}


