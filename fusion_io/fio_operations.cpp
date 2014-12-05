#include "fusion_io.h"
#include <iostream>


fio_field_sum::fio_field_sum(const fio_field* f1, const fio_field* f2)
{
  f[0] = f1->clone();
  f[1] = f2->clone();

  if(f1->dimension() != f2->dimension()) {
    std::cerr << "Error: incompatible dimensions" << std::endl;
  } 
}

fio_field_sum::~fio_field_sum()
{
  delete(f[0]);
  delete(f[1]);
}

int fio_field_sum::dimension() const
{
  return f[0]->dimension();
}


int fio_field_sum::eval(const double* d, double *v)
{
  double *v1, *v2;
  int result;

  v1 = new double[dimension()];
  v2 = new double[dimension()];

  result = f[0]->eval(d, v1);
  if(result != FIO_SUCCESS) return result;
  result = f[1]->eval(d, v2);
  if(result != FIO_SUCCESS) return result;

  for(int i=0; i<dimension(); i++)
    v[i] = v1[i] + v2[i];
 
  delete[] v1;
  delete[] v2;

  return FIO_SUCCESS;
}

fio_field_product::fio_field_product(const fio_field* f1, const fio_field* f2)
{
  f[0] = f1->clone();
  f[1] = f2->clone();

  if(f1->dimension() != f2->dimension()) {
    std::cerr << "Error: incompatible dimensions" << std::endl;
  } 
}

fio_field_product::~fio_field_product()
{
  delete(f[0]);
  delete(f[1]);
}

int fio_field_product::dimension() const
{
  return 1;
}

int fio_field_product::eval(const double* d, double *v)
{
  double *v1, *v2;
  int result;

  v1 = new double[f[0]->dimension()];
  v2 = new double[f[1]->dimension()];


  result = f[0]->eval(d, v1);
  if(result != FIO_SUCCESS) return result;
  result = f[1]->eval(d, v2);
  if(result != FIO_SUCCESS) return result;

  *v = 0;
  for(int i=0; i<dimension(); i++)
    *v += v1[i]*v2[i];

  delete[] v1;
  delete[] v2;
 
  return FIO_SUCCESS;
}


