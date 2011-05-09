#include "m3dc1_source.h"
#include <iostream>

m3dc1_source::m3dc1_source()
  : next_elm(0), psi(0), f(0), g(0)
{
  filename = "C1.h5";
  time = -1;
  last_elm = 0;
  factor = 1.;
  misses = 0;
  hits = 0;
  evals = 0;
  lost_shifts = 0;
}

m3dc1_source::m3dc1_source(std::string f, int t)
  : next_elm(0), psi(0), f(0), g(0)
{
  filename = f;
  time = t;
  last_elm = 0;
  factor = 1.;
  misses = 0;
  hits = 0;
  evals = 0;
  lost_shifts = 0;
}

m3dc1_source::~m3dc1_source()
{
  std::cerr << "hits = " << hits
	    << " (" << 100.*(double)hits/(double)evals << "%)\n"
	    << "misses = " << misses
	    << " (" << 100.*(double)misses/(double)evals << "%)\n"
	    << "lost shifts = " << lost_shifts 
	    << " (" << 100.*(double)lost_shifts/(double)evals << "%)\n"
	    << "evals = " << evals << '\n'
	    << "hits+misses = " << hits+misses
	    << std::endl;

  for(int d=0; d<memory_depth; d++)
    delete[] next_elm[d];

  if(next_elm) delete[] next_elm;
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

  next_elm = new int*[memory_depth];
  for(int d=0; d<memory_depth; d++) {
    next_elm[d] = new int[psi->mesh->nelms];
    for(int i=0; i<psi->mesh->nelms; i++) 
      next_elm[d][i] = -1;
  }

  return true;
}

bool m3dc1_source::eval(const double r, const double phi, const double z,
			double* b_r, double* b_phi, double* b_z)
{
  int elm;
  int guess[memory_depth+1];
  int nguess;

  const m3dc1_field::m3dc1_get_op psiget = (m3dc1_field::m3dc1_get_op)
    (m3dc1_field::GET_DVAL);

  const m3dc1_field::m3dc1_get_op gget = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_VAL);

  const m3dc1_field::m3dc1_get_op fget = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);


  double val[m3dc1_field::OP_NUM];

  if(last_elm < 0) std::cerr << "ARRRGH" << std::endl;

  // guess which element contains (r, phi, z)
  // first guess is the element that cointained the previous point
  guess[0] = last_elm;
  nguess = 1;

  // next guess is the "next_elm" of the element that contained 
  // the previous point
  for(int d=0; d<memory_depth; d++) {
    if(next_elm[d][last_elm] == -1) break;
    guess[d+1] = next_elm[d][last_elm];
    nguess++;
  }

  // calculate psi
  if(!(psi->eval(r, phi, z, psiget, val, &elm, guess, nguess)))
    return false;

  // if first guess was wrong, then set the "next_elm" of the element that
  // contained the previous point to be the correct answer
  if(elm != guess[0]) {
    // if there is room in memory, add new elm to memory

    // make sure elm is not already in memory
    bool was_guess = false;
    for(int d=0; d<memory_depth; d++) {
      if(elm==guess[d+1]) {
	was_guess = true;
	break;
      }
    }
    if(!was_guess) {
      misses++;
      if(next_elm[memory_depth-1][last_elm] != -1)
	lost_shifts++;
      for(int d=1; d<memory_depth; d++) 
	next_elm[d][last_elm] = next_elm[d-1][last_elm];
      next_elm[0][last_elm] = elm;
    } else hits++;
    last_elm = elm;
  } else hits++;

  evals++;

  // B_R   = -(dpsi/dZ)/R - (d2f/dRdphi)
  // B_Z   =  (dpsi/dR)/R - (d2f/dZdphi)
  // B_Phi =  F/R

  *b_r -= factor*val[m3dc1_field::OP_DZ]/r;
  *b_z += factor*val[m3dc1_field::OP_DR]/r;

  if(!g->eval(r, phi, z, gget, val, &elm, &last_elm, 1))
    return false;
  *b_phi += factor*val[m3dc1_field::OP_1]/r;

  if(!f->eval(r, phi, z, fget, val, &elm, &last_elm, 1))
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
