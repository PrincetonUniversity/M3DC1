#include "fusion_io.h"
#include <iostream>

int m3dc1_fio_field::load(const fio_option_list* opt)
{
  int ilin;

  opt->get_option(FIO_TIMESLICE, &time);
  opt->get_option(FIO_LINEAR_SCALE, &linfac);
  opt->get_option(FIO_PERTURBED_ONLY, &ilin);

  extsub = (source->extsubtract==1);
  eqsub = (source->eqsubtract==1) && (ilin==0);
  use_f = (source->i3d==1 || source->icomplex==1);

  if(time==-1) {
    eqsub = false;   // equilibrium fields need not be added in
    use_f = false;   // equilibrium is assumed axisymmetric
  }

  if(linfac != 1.) {
    if(source->linear == 0) {
      std::cerr << "Linear scale linfac is ignored for nonlinear data."
		<< std::endl;
      linfac = 1.;
    }
    if(time==-1) {
      std::cerr << "Linear scale linfac is ignored for equilibrium data." 
		<< std::endl;
      linfac = 1.;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_scalar_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  f1 = source->file.load_field(name.c_str(), time);
  if(!f1) return 1;

  if(eqsub) {
    f0 = source->file.load_field(name.c_str(), -1);
    if(!f0) return 1;
  }

  return FIO_SUCCESS;
}


int m3dc1_scalar_field::eval(const double* x, double* v)
{
  const m3dc1_field::m3dc1_get_op get = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  double val[m3dc1_field::OP_NUM];

  if(!f1->eval(x[0], x[1], x[2], get, val)) {
    return FIO_OUT_OF_BOUNDS;
  }
  *v = linfac*val[m3dc1_field::OP_1];

  if(eqsub) {
    if(!f0->eval(x[0], x[1], x[2], get, val)) {
      return FIO_OUT_OF_BOUNDS;
    }
    *v += val[m3dc1_field::OP_1];
  }

  *v *= factor;

  return FIO_SUCCESS;
}


int m3dc1_phi_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  return FIO_SUCCESS;
}

int m3dc1_phi_field::eval(const double* x, double* v)
{
  *v = 0.;

  return FIO_SUCCESS;
}

int m3dc1_electric_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  return FIO_SUCCESS;
}

int m3dc1_electric_field::eval(const double* x, double* v)
{
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;

  return FIO_SUCCESS;
}

int m3dc1_vector_potential::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_vector_potential::eval(const double* x, double* v)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op fget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  double val[m3dc1_field::OP_NUM];

  // A_R   =  R (df/dZ)
  // A_Z   = -R (df/dR) - F0 ln(R)
  // A_Phi = psi/R

  if(!psi1->eval(x[0], x[1], x[2], psiget, val))
    return FIO_OUT_OF_BOUNDS;

  v[1] = linfac*val[m3dc1_field::OP_1]/x[0];

  if(use_f) {
    if(!f1->eval(x[0], x[1], x[2], fget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] =  linfac*x[0]*val[m3dc1_field::OP_DZ];
    v[2] = -linfac*x[0]*val[m3dc1_field::OP_DR];
  } else {
    v[0] = v[2] = 0.;
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[1] += val[m3dc1_field::OP_1]/x[0];

    v[2] -= source->bzero*source->rzero*log(x[0]);
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[1] += linfac*val[m3dc1_field::OP_1]/x[0];

    if(use_f) {
      if(!fx->eval(x[0], x[1], x[2], fget, val))
        return FIO_OUT_OF_BOUNDS;

      v[0] += linfac*x[0]*val[m3dc1_field::OP_DZ];
      v[2] -= linfac*x[0]*val[m3dc1_field::OP_DR];
    }
  }

  // convert to mks
  v[0] *= source->B0*source->L0;
  v[1] *= source->B0*source->L0;
  v[2] *= source->B0*source->L0;
  
  return FIO_SUCCESS;
}

int m3dc1_magnetic_field::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_magnetic_field::eval(const double* x, double* v)
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

  if(!psi1->eval(x[0], x[1], x[2], psiget, val))
    return FIO_OUT_OF_BOUNDS;

  v[0] = -linfac*val[m3dc1_field::OP_DZ]/x[0];
  v[2] =  linfac*val[m3dc1_field::OP_DR]/x[0];

  if(!i1->eval(x[0], x[1], x[2], gget, val))
    return FIO_OUT_OF_BOUNDS;
  v[1] =  linfac*val[m3dc1_field::OP_1]/x[0];

  if(use_f) {
    if(!f1->eval(x[0], x[1], x[2], fget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DRP];
    v[2] -= linfac*val[m3dc1_field::OP_DZP];
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= val[m3dc1_field::OP_DZ]/x[0];
    v[2] += val[m3dc1_field::OP_DR]/x[0];

    if(!i0->eval(x[0], x[1], x[2], gget, val))
      return FIO_OUT_OF_BOUNDS;
    v[1] += val[m3dc1_field::OP_1]/x[0];
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DZ]/x[0];
    v[2] += linfac*val[m3dc1_field::OP_DR]/x[0];

    if(!ix->eval(x[0], x[1], x[2], gget, val))
      return FIO_OUT_OF_BOUNDS;
    v[1] += linfac*val[m3dc1_field::OP_1]/x[0];

    if(use_f) {
      if(!fx->eval(x[0], x[1], x[2], fget, val))
        return FIO_OUT_OF_BOUNDS;

      v[0] -= linfac*val[m3dc1_field::OP_DRP];
      v[1] -= linfac*val[m3dc1_field::OP_DZP];
    }
  }

  v[0] *= source->B0;
  v[1] *= source->B0;
  v[2] *= source->B0;
  
  return FIO_SUCCESS;
}


int m3dc1_current_density::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_current_density::eval(const double* x, double* v)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_DDVAL);

  const m3dc1_field::m3dc1_get_op gget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_PPVAL);

  double val[m3dc1_field::OP_NUM];

  // J_R   = -(d(F+f'')/dZ)/R + (d(psi')/dR)/R^2
  // J_Z   =  (d(F+f'')/dR)/R + (d(psi')/dZ)/R^2
  // B_Phi = -Del*[psi]/R

  if(!psi1->eval(x[0], x[1], x[2], psiget, val))
    return FIO_OUT_OF_BOUNDS;

  v[0] =  linfac*val[m3dc1_field::OP_DRP]/(x[0]*x[0]);
  v[2] =  linfac*val[m3dc1_field::OP_DZP]/(x[0]*x[0]);
  v[1] = -linfac*(val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ] -
		  val[m3dc1_field::OP_DR]/x[0])/x[0]; 

  if(!i1->eval(x[0], x[1], x[2], gget, val))
    return FIO_OUT_OF_BOUNDS;
  v[0] -= linfac*val[m3dc1_field::OP_DZ]/x[0];
  v[2] += linfac*val[m3dc1_field::OP_DR]/x[0];

  if(use_f) {
    if(!f1->eval(x[0], x[1], x[2], fget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] -= linfac*val[m3dc1_field::OP_DZPP]/x[0];
    v[2] += linfac*val[m3dc1_field::OP_DRPP]/x[0];
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] +=  val[m3dc1_field::OP_DRP]/(x[0]*x[0]);
    v[2] +=  val[m3dc1_field::OP_DZP]/(x[0]*x[0]);
    v[1] -= (val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ] -
	     val[m3dc1_field::OP_DR]/x[0])/x[0]; 

    if(!i0->eval(x[0], x[1], x[2], gget, val))
      return FIO_OUT_OF_BOUNDS;
    v[0] -= val[m3dc1_field::OP_DZ]/x[0];
    v[2] += val[m3dc1_field::OP_DR]/x[0];
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] += linfac*val[m3dc1_field::OP_DRP]/(x[0]*x[0]);
    v[2] += linfac*val[m3dc1_field::OP_DZP]/(x[0]*x[0]);
    v[1] -= linfac*(val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DZZ] -
		    val[m3dc1_field::OP_DR]/x[0])/x[0]; 

    if(!ix->eval(x[0], x[1], x[2], gget, val))
      return FIO_OUT_OF_BOUNDS;
    v[0] -= linfac*val[m3dc1_field::OP_DZ]/x[0];
    v[2] += linfac*val[m3dc1_field::OP_DR]/x[0];

    if(use_f) {
      if(!fx->eval(x[0], x[1], x[2], fget, val))
        return FIO_OUT_OF_BOUNDS;

      v[0] -= linfac*val[m3dc1_field::OP_DZPP]/x[0];
      v[2] += linfac*val[m3dc1_field::OP_DRPP]/x[0];
    }
  }

  v[0] *= source->J0;
  v[1] *= source->J0;
  v[2] *= source->J0;
  
  return FIO_SUCCESS;
}



int m3dc1_grad_vector_potential::load(const fio_option_list* opt)
{
  m3dc1_fio_field::load(opt);

  psi1 = source->file.load_field("psi", time);
  if(!psi1) return 1;
  i1 = source->file.load_field("I", time);
  if(!i1) return 1;
  if(use_f) {
    f1 = source->file.load_field("f", time);
    if(!f1) return 1;
  }

  if(eqsub) {
    psi0 = source->file.load_field("psi", -1);
    if(!psi0) return 1;
    i0 = source->file.load_field("I", -1);
    if(!i0) return 1;
  }

  if(extsub) {
    psix = source->file.load_field("psi_ext", time);
    if(!psix) return 1;
    ix = source->file.load_field("I_ext", time);
    if(!ix) return 1;
    if(use_f) {
      fx = source->file.load_field("f_ext", time);
      if(!fx) return 1;
    }
  }

  return FIO_SUCCESS;
}


int m3dc1_grad_vector_potential::eval(const double* x, double* v)
{
  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);

  const m3dc1_field::m3dc1_get_op fget =
    (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL | m3dc1_field::GET_DDVAL);

  double val[m3dc1_field::OP_NUM];

  if(!psi1->eval(x[0], x[1], x[2], psiget, val))
    return FIO_OUT_OF_BOUNDS;

  v[1] = linfac*(val[m3dc1_field::OP_DR] - val[m3dc1_field::OP_1]/x[0])/x[0];
  v[4] = linfac*val[m3dc1_field::OP_DP]/x[0];
  v[7] = linfac*val[m3dc1_field::OP_DZ]/x[0];

  if(use_f) {
    if(!f1->eval(x[0], x[1], x[2], fget, val))
      return FIO_OUT_OF_BOUNDS;

    v[0] =  linfac*(x[0]*val[m3dc1_field::OP_DRZ] + val[m3dc1_field::OP_DZ]);
    v[2] = -linfac*(x[0]*val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DR]);
    v[3] =  linfac*(x[0]*val[m3dc1_field::OP_DZP]);
    v[5] = -linfac*(x[0]*val[m3dc1_field::OP_DRP]);
    v[6] =  linfac*(x[0]*val[m3dc1_field::OP_DZZ]);
    v[8] = -linfac*(x[0]*val[m3dc1_field::OP_DRZ]);
  } else {
    v[0] = v[2] = v[3] = v[5] = v[6] = v[8] = 0.;
  }

  if(eqsub) {
    if(!psi0->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[1] += (val[m3dc1_field::OP_DR]-val[m3dc1_field::OP_1]/x[0])/x[0];
    v[4] += val[m3dc1_field::OP_DP]/x[0];
    v[7] += val[m3dc1_field::OP_DZ]/x[0];

    v[2] -= source->bzero*source->rzero/x[0];
  }

  if(extsub) {
    if(!psix->eval(x[0], x[1], x[2], psiget, val))
      return FIO_OUT_OF_BOUNDS;

    v[1] += linfac*(val[m3dc1_field::OP_DR]-val[m3dc1_field::OP_1]/x[0])/x[0];
    v[4] += linfac*val[m3dc1_field::OP_DP]/x[0];
    v[7] += linfac*val[m3dc1_field::OP_DZ]/x[0];

    if(use_f) {
      if(!fx->eval(x[0], x[1], x[2], fget, val))
        return FIO_OUT_OF_BOUNDS;

      v[0] += linfac*(x[0]*val[m3dc1_field::OP_DRZ] + val[m3dc1_field::OP_DZ]);
      v[2] -= linfac*(x[0]*val[m3dc1_field::OP_DRR] + val[m3dc1_field::OP_DR]);
      v[3] += linfac*(x[0]*val[m3dc1_field::OP_DZP]);
      v[5] -= linfac*(x[0]*val[m3dc1_field::OP_DRP]);
      v[6] += linfac*(x[0]*val[m3dc1_field::OP_DZZ]);
      v[8] -= linfac*(x[0]*val[m3dc1_field::OP_DRZ]);
    }
  }

  // convert to mks
  v[0] *= source->B0;
  v[1] *= source->B0;
  v[2] *= source->B0;
  v[3] *= source->B0*source->L0;
  v[4] *= source->B0*source->L0;
  v[5] *= source->B0*source->L0;
  v[6] *= source->B0;
  v[7] *= source->B0;
  v[8] *= source->B0;
  
  return FIO_SUCCESS;
}
