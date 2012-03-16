#include "m3dc1_field.h"
#include <math.h>

double pow(double x, int p)
{
  double r = 1.;
  if(p < 0) return 0.;
  for(int i=0; i<p; i++) r = r*x;
  return r;
}

const int m3dc1_field::mi[m3dc1_field::nbasis] = 
  { 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0 };
const int m3dc1_field::ni[m3dc1_field::nbasis] = 
  { 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5 };
const int m3dc1_3d_field::li[m3dc1_3d_field::tbasis] =
  { 0,1,2,3 };

m3dc1_field::m3dc1_field(m3dc1_mesh* m)
{
  mesh = m;
  data = new double[mesh->nelms*nbasis];
}

m3dc1_field::~m3dc1_field()
{
  delete[] data;
}

bool m3dc1_field::eval(const double r, const double phi, const double z,
		       const m3dc1_get_op op, double* val, int* element)
{
  double xi, eta, temp, co2, sn2, cosn;
  double xipow[10], etapow[10];
  double *xip, *etap;
  int e;
  
  int guess = (element ? *element : -1);

  e = mesh->in_element(r, phi, z, &xi, 0, &eta, guess);
  
  if(element) *element = e;

  if(e < 0) return false;
  
  // create array of powers (xip[i] = xi^i, etc..)
  xip = &(xipow[2]);
  etap = &(etapow[2]);
  for(int i=-2; i<=0; i++) {
    xip[i] = 1.;
    etap[i] = 1.;
  }
  for(int i=1; i<=5; i++) {
    xip[i] = xip[i-1]*xi;
    etap[i] = etap[i-1]*eta;
  }

  for(int i=0; i<OP_NUM; i++) val[i] = 0.;
  
  co2  = mesh->co[e]*mesh->co[e];
  sn2  = mesh->sn[e]*mesh->sn[e];
  cosn = mesh->co[e]*mesh->sn[e];

  for(int p=0, j=e*nbasis; p<nbasis; p++, j++) {
      
    if((op & GET_VAL) == GET_VAL) {
      val[OP_1] = val[OP_1] + data[j]*xip[mi[p]]*etap[ni[p]];
    }
      
    if((op & GET_DVAL) == GET_DVAL) {
      temp = data[j]*mi[p]*xip[mi[p]-1]*etap[ni[p]];
      val[OP_DR] += mesh->co[e]*temp;
      val[OP_DZ] += mesh->sn[e]*temp;
      
      temp = data[j]*ni[p]*xip[mi[p]]*etap[ni[p]-1];
      val[OP_DR] -= mesh->sn[e]*temp;
      val[OP_DZ] += mesh->co[e]*temp;	
    }

    if((op & GET_DDVAL) == GET_DDVAL) {
      temp = data[j]*(mi[p]-1)*mi[p]*xip[mi[p]-2]*etap[ni[p]];
      val[OP_DRR] +=  co2*temp;
      val[OP_DRZ] += cosn*temp;
      val[OP_DZZ] +=  sn2*temp;

      temp = data[j]*mi[p]*ni[p]*xip[mi[p]-1]*etap[ni[p]-1];
      val[OP_DRR] -=     2.*cosn*temp;
      val[OP_DRZ] += (co2 - sn2)*temp;
      val[OP_DZZ] +=     2.*cosn*temp;

      temp = data[j]*(ni[p]-1)*ni[p]*xip[mi[p]]*etap[ni[p]-2];
      val[OP_DRR] +=  sn2*temp;
      val[OP_DRZ] -= cosn*temp;
      val[OP_DZZ] +=  co2*temp;
    }
  }
  
  return true;
}

m3dc1_complex_field::m3dc1_complex_field(m3dc1_mesh* m, int n)
  : m3dc1_field(m)
{
  ntor = n;
  data_i = new double[mesh->nelms*nbasis];
}

m3dc1_complex_field::~m3dc1_complex_field()
{
  delete[] data_i;
}

bool m3dc1_complex_field::eval(const double r, const double phi, 
			       const double z,
			       const m3dc1_get_op op, double* val,
			       int* element)
{
  m3dc1_get_op new_op = (m3dc1_get_op)(op - (op & GET_PVAL) 
				       - (op & GET_PPVAL));

  double val_r[OP_NUM];
  double val_i[OP_NUM];
  double* temp;
  double co, sn;

  // get real values
  if(!m3dc1_field::eval(r,phi,z,new_op,val_r,element))
    return false;

  // get imaginary values
  temp = data;
  data = data_i;
  bool result = 
    m3dc1_field::eval(r,phi,z,new_op,val_i,element);
  data = temp;

  if(!result) 
    return false;

  co = cos(phi*ntor);
  sn = sin(phi*ntor);

  for(int i=0; i<OP_DP; i++)
    val[i] = val_r[i]*co - val_i[i]*sn;

  if((op & GET_PVAL) == GET_PVAL)
    for(int i=0; i<OP_DP; i++)
      val[i+OP_DP] = -ntor*(val_i[i]*co + val_r[i]*sn);

  if((op & GET_PPVAL) == GET_PPVAL)
    for(int i=0; i<OP_DP; i++) 
      val[i+OP_DPP] = -ntor*ntor*val[i];

  return true;
}

m3dc1_3d_field::m3dc1_3d_field(m3dc1_mesh* m)
{
  mesh = m;
  data = new double[mesh->nelms*nbasis*tbasis];
}

bool m3dc1_3d_field::eval(const double r, const double phi, const double z,
			  const m3dc1_get_op op, double* val, int* element)
{
  double xi, zi, eta, temp, co2, sn2, cosn;
  double v[6];
  int e;

  e = mesh->in_element(r, phi, z, &xi, &zi, &eta);
  
  if(element) *element = e;

  if(e < 0) return false;
  
  for(int i=0; i<OP_NUM; i++) val[i] = 0.;
  
  co2  = mesh->co[e]*mesh->co[e];
  sn2  = mesh->sn[e]*mesh->sn[e];
  cosn = mesh->co[e]*mesh->sn[e];

  int j = e*nbasis*tbasis;
  for(int q=0; q<tbasis; q++) {
    for(int p=0; p<nbasis; p++) {      
      if((op & GET_VAL) == GET_VAL) {
	temp = data[j]*pow(xi,mi[p])*pow(eta,ni[p]);
	v[OP_1] = temp;
      }
      
      if((op & GET_DVAL) == GET_DVAL) {
	temp = data[j]*mi[p]*pow(xi,mi[p]-1)*pow(eta,ni[p]);
	v[OP_DR] = mesh->co[e]*temp;
	v[OP_DZ] = mesh->sn[e]*temp;
	
	temp = data[j]*ni[p]*pow(xi,mi[p])*pow(eta,ni[p]-1);
	v[OP_DR] -= mesh->sn[e]*temp;
	v[OP_DZ] += mesh->co[e]*temp;	
      }

      if((op & GET_DDVAL) == GET_DDVAL) {
	temp = data[j]*(mi[p]-1)*mi[p]*pow(xi,mi[p]-2)*pow(eta,ni[p]);
	v[OP_DRR] =  co2*temp;
	v[OP_DRZ] = cosn*temp;
	v[OP_DZZ] =  sn2*temp;
	
	temp = data[j]*mi[p]*ni[p]*pow(xi,mi[p]-1)*pow(eta,ni[p]-1);
	v[OP_DRR] -=     2.*cosn*temp;
	v[OP_DRZ] += (co2 - sn2)*temp;
	v[OP_DZZ] +=     2.*cosn*temp;
	
	temp = data[j]*(ni[p]-1)*ni[p]*pow(xi,mi[p])*pow(eta,ni[p]-2);
	v[OP_DRR] +=  sn2*temp;
	v[OP_DRZ] -= cosn*temp;
	v[OP_DZZ] +=  co2*temp;
      }

      temp = pow(zi, li[q]);
      if((op & GET_VAL) == GET_VAL) {
	val[OP_1  ] += v[OP_1  ]*temp;
      }
      if((op & GET_DVAL) == GET_DVAL) {
	val[OP_DR ] += v[OP_DR ]*temp;
	val[OP_DZ ] += v[OP_DZ ]*temp;
      }
      if((op & GET_DDVAL) == GET_DDVAL) {
	val[OP_DRR] += v[OP_DRR]*temp;
	val[OP_DRZ] += v[OP_DRZ]*temp;
	val[OP_DZZ] += v[OP_DZZ]*temp;
      }

      if((op & GET_PVAL) == GET_PVAL) {
	temp = pow(zi, li[q]-1)*li[q];
	if((op & GET_VAL) == GET_VAL) {
	  val[OP_DP  ] += v[OP_1  ]*temp;
	}
	if((op & GET_DVAL) == GET_DVAL) {
	  val[OP_DRP ] += v[OP_DR ]*temp;
	  val[OP_DZP ] += v[OP_DZ ]*temp;
	}
	if((op & GET_DDVAL) == GET_DDVAL) {
	  val[OP_DRRP] += v[OP_DRR]*temp;
	  val[OP_DRZP] += v[OP_DRZ]*temp;
	  val[OP_DZZP] += v[OP_DZZ]*temp;
	}
      }

      if((op & GET_PPVAL) == GET_PPVAL) {
	temp = pow(zi, li[q]-2)*li[q]*(li[q]-1);
	if((op & GET_VAL) == GET_VAL) {
	  val[OP_DPP  ] += v[OP_1  ]*temp;
	}
	if((op & GET_DVAL) == GET_DVAL) {
	  val[OP_DRPP ] += v[OP_DR ]*temp;
	  val[OP_DZPP ] += v[OP_DZ ]*temp;
	}
	if((op & GET_DDVAL) == GET_DDVAL) {
	  val[OP_DRRPP] += v[OP_DRR]*temp;
	  val[OP_DRZPP] += v[OP_DRZ]*temp;
	  val[OP_DZZPP] += v[OP_DZZ]*temp;
	}
      }

      j++;
    }
  }

  return true;
}

bool m3dc1_compound_field::eval(const double r,const double phi,const double z,
				const m3dc1_get_op op,double* val,int* element)
{
  double temp[OP_NUM];

  for(int i=0; i<OP_NUM; i++) val[i] = 0.;

  subfield_list::iterator j = subfield.begin();

  while(j != subfield.end()) {
    if(!(*j)->eval(r,phi,z,op,temp,element)) return false;
    for(int i=0; i<OP_NUM; i++) val[i] += temp[i];
    j++;
  }

  return true;
}



m3dc1_timeslice::m3dc1_timeslice()
  : mesh(0)
{
  is_3d = 0;
}

m3dc1_timeslice::~m3dc1_timeslice()
{
  std::map<std::string, m3dc1_field*>::iterator i = field_map.begin();
  while(i != field_map.end()) {
    delete(i->second);
    i++;
  }

  if(mesh) delete(mesh);
}
