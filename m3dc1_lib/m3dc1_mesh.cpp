#include "m3dc1_mesh.h"

#include <iostream>

m3dc1_mesh::m3dc1_mesh(int n)
{
  a = new double[n];
  b = new double[n];
  c = new double[n];
  co = new double[n];
  sn = new double[n];
  x = new double[n];
  z = new double[n];
  bound = new int[n];
  nelms = n;
  memory_depth = 0;
  last_elm = -1;
  next_elm = 0;
  hits = 0;
  misses = 0;
}

m3dc1_mesh::~m3dc1_mesh()
{
  /*
  int evals = hits+misses;
  std::cerr << "hits = " << hits
	    << " (" << 100.*(double)hits/(double)evals << "%)\n"
	    << "misses = " << misses
	    << " (" << 100.*(double)misses/(double)evals << "%)\n"
	    << "hits+misses = " << evals
	    << std::endl;
  */
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] co;
  delete[] sn;
  delete[] x;
  delete[] z;
  delete[] bound;
  clear_memory();
}

void m3dc1_mesh::extent(double* X0, double* X1,
			double* Phi0, double* Phi1,
			double* Z0, double* Z1) const
{
  if(nelms==0) {
    *X0 = *X1 = *Phi0 = *Phi1 = *Z0 = *Z1 = 0.;
    return;
  }

  double R[3], Z[3];

  *X0 = x[0];
  *X1 = x[0];
  *Z0 = z[0];
  *Z1 = z[0];

  for(int i=0; i<nelms; i++) {
    R[0] = x[i];
    R[1] = x[i] + (a[i]+b[i])*co[i];
    R[2] = x[i] + b[i]*co[i] - c[i]*sn[i];
    Z[0] = z[i];
    Z[1] = z[i] + (a[i]+b[i])*sn[i];
    Z[2] = z[i] + c[i]*co[i] + b[i]*sn[i];

    for(int j=0; j<3; j++) {
      if(R[j] < *X0) *X0 = R[j];
      if(R[j] > *X1) *X1 = R[j];
      if(Z[j] < *Z0) *Z0 = Z[j];
      if(Z[j] > *Z1) *Z1 = Z[j];
    }
  }

  *Phi0 = 0.;
  *Phi1 = 0.;
}

void m3dc1_mesh::clear_memory()
{
  if(next_elm) {
    for(int d=0; d<memory_depth; d++)
      delete[] next_elm[d];
    
    delete[] next_elm;

    next_elm = 0;
  }
  last_elm = 0;
  memory_depth = 0;
}

bool m3dc1_mesh::set_memory_depth(const int d)
{
  if(memory_depth) clear_memory();

  memory_depth = d;
  if(memory_depth==0) return true;

  next_elm = new int*[memory_depth];
  for(int d=0; d<memory_depth; d++) {
    next_elm[d] = new int[nelms];
    for(int i=0; i<nelms; i++) 
      next_elm[d][i] = -1;
  }

  return true;
}

int m3dc1_mesh::in_element(double X, double Phi, double Z, 
			   double* xi, double* zi, double* eta,
			   int guess)
{
  // if a guess is provided, test it
  if(guess >= 0) {
    if(is_in_element(guess,X,Phi,Z,xi,zi,eta)) {
      hits++;
      last_elm = guess;
      return guess;
    }
  }

  if(last_elm >= 0) {
    // first, check last elm
    if(is_in_element(last_elm,X,Phi,Z,xi,zi,eta)) {
      hits++;
      return last_elm;
    }

    // now, check each element in memory
    for(int m=0; m<memory_depth; m++) {
      int e = next_elm[m][last_elm];
      if(e==-1) break;
      if(is_in_element(e,X,Phi,Z,xi,zi,eta)) {
	hits++;
	last_elm = e;
	return e;
      }
    }
  }

  misses++;

  // Now, search randomly.
  for(int e=0; e<nelms; e++) {
    if(is_in_element(e,X,Phi,Z,xi,zi,eta)) {     
      // append this to memory
      if(memory_depth>0 && last_elm >= 0) {
	for(int m=memory_depth-1; m>0; m--)
	  next_elm[m][last_elm] = next_elm[m-1][last_elm];
	next_elm[0][last_elm] = e;
      }
      last_elm = e;
      return e;
    }
  }

  // failed to find an elm containing coordinates.
  return -1;
}



m3dc1_3d_mesh::m3dc1_3d_mesh(int n)
  : m3dc1_mesh(n)
{
  d = new double[n];
  phi = new double[n];
}

m3dc1_3d_mesh::~m3dc1_3d_mesh()
{
  delete[] d;
  delete[] phi;
}
