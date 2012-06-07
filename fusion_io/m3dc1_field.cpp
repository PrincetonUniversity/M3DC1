#include "m3dc1_field.h"
#include "fusion_io_defs.h"

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

