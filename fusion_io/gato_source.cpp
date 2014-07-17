#include "fusion_io.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

gato_source::gato_source()
  : r_bound(0), t_bound(0), 
    psival(0), pressure(0), ftor(0), pprime(0), ffprime(0), rcc(0), zcc(0),
    psimesh(0), dpsidr(0), dpsidz(0)
{
}

gato_source::~gato_source()
{
  close();
}

int gato_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  opt->add_option(FIO_LINEAR_SCALE, 1.);
  opt->add_option(FIO_PART, FIO_TOTAL);

  return FIO_SUCCESS;
}

int gato_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_TOTAL_PRESSURE);

  return FIO_SUCCESS;
}

int gato_source::get_field(const field_type t,fio_field** f,
			     const fio_option_list* opt)
{
  int part;

  opt->get_option(FIO_PART, &part);
  if(part != FIO_TOTAL && 
     part != FIO_EQUILIBRIUM_ONLY &&
     part != FIO_PERTURBED_ONLY) {
    std::cerr << "Unrecognized option for FIO_PART: " << part << std::endl;
    return FIO_UNSUPPORTED;
  }

  switch(t) {
  case(FIO_MAGNETIC_FIELD):
    *f = new gato_magnetic_field(this, part);
    break;

  case(FIO_TOTAL_PRESSURE):
    *f = new gato_pressure_field(this, part);
    break;

  default:
    return FIO_UNSUPPORTED;
  };

  return FIO_SUCCESS;
}

int gato_source::get_element(const double* x, int* e)
{
  double r = sqrt((x[0]-xma)*(x[0]-xma) + (x[2]-zma)*(x[2]-zma));
  double t = atan2(x[2]-zma,x[0]-xma);

  for(int i=0; i<(jpsi-1)*itht; i++) {
    if(r >= r_bound[i] && r <= r_bound[i+1]) {
      int l = i+jpsi;
      if(l > jpsi*itht) l -= jpsi*itht;
      if(t >= t_bound[i] && t <= t_bound[l]) {
	*e = i;
	return FIO_SUCCESS;
      }
    }
  };

  return FIO_OUT_OF_BOUNDS;
}

int gato_source::interpolate_psi(const double* x, double* psi, int* eout)
{
  int e, ierr;
  

  ierr = get_element(x, &e);
  if(ierr != FIO_SUCCESS) {
    std::cerr << "Couldn't find element" << std::endl; 
    return ierr;
  }

  double dr = x[0] - rcc[e];
  double dz = x[2] - zcc[e];
  *psi = psimesh[e]*(psival[jpsi-1]-psival[0]) + psival[0]
    + dr*dpsidr[e] + dz*dpsidz[e];

  if(eout) 
    *eout = e;

  return FIO_SUCCESS;
}

int gato_source::interpolate_flux_function(const double* f, const double psi0,
					   double* f0) const
{
  for(int i=0; i<jpsi-1; i++) {
    if(psi0 >= psival[i] && psi0 <= psival[i+1]) {
      *f0 = (f[i]*(psival[i+1]-psi0) + f[i+1]*(psi0-psival[i]))
	/(psival[i+1]-psival[i]);
      return FIO_SUCCESS;
    }
    if(psi0 <= psival[i] && psi0 >= psival[i+1]) {
      *f0 = (f[i]*(psival[i+1]-psi0) + f[i+1]*(psi0-psival[i]))
	/(psival[i+1]-psival[i]);
      return FIO_SUCCESS;
    }
  }
  std::cerr << "Error interpolating from psi" << std::endl;
  return FIO_OUT_OF_BOUNDS;
}

int gato_source::scan_until(std::ifstream& ifile, const char* line)
{
  const int nsize = 256;
  char buff[nsize];
  int size = strlen(line);
  
  if(size > nsize) {
    std::cerr << "Warning: truncating line " << line;
    size = nsize;
  }

  do {
    ifile.getline(buff, nsize);
    if(ifile.eof()) {
      std::cerr << "Error: " << line << " not found." << std::endl;
      return FIO_FILE_ERROR;
    }

    if(ifile.fail()) 
      ifile.clear();

  } while(strncmp(line, buff, size) != 0);

  return FIO_SUCCESS;
}

int gato_source::open(const char* filename)
{
  std::ifstream ifile;
  
  ifile.open(filename, std::ifstream::in);
  if(!ifile)
    return FIO_FILE_ERROR;

  if(scan_until(ifile, "  ntor nev  ncase norm") != FIO_SUCCESS)
    return FIO_FILE_ERROR;
  ifile >> ntor;
  std::cout << "ntor = " << ntor << std::endl;

  if(scan_until(ifile, "   jpsi     itht") != FIO_SUCCESS)
    return FIO_FILE_ERROR;
  ifile >> jpsi >> itht;
  std::cout << "jpsi = " << jpsi << ", itht = " << itht << std::endl;

  if(scan_until(ifile, " psimax xma zma") != FIO_SUCCESS)
    return FIO_FILE_ERROR;
  ifile >> psimax >> xma >> zma;
  std::cout << "psimax = " << psimax
	    << ", xma = " << xma
	    << ", zma = " << zma << std::endl;


  psival = new double[jpsi+1];
  pressure = new double[jpsi+1];
  ftor = new double[jpsi+1];
  pprime = new double[jpsi+1];
  ffprime = new double[jpsi+1];
  rcc = new double[jpsi*itht];
  zcc = new double[jpsi*itht];
  dpsidr = new double[jpsi*itht];  
  dpsidz = new double[jpsi*itht];
  psimesh = new double[jpsi*itht];

  if(scan_until(ifile, " psival(j)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<=jpsi; i++) ifile >> psival[i];

  if(scan_until(ifile, " pressure(j)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<=jpsi; i++) ifile >> pressure[i];

  if(scan_until(ifile, " ftor(j)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<=jpsi; i++) ifile >> ftor[i];

  if(scan_until(ifile, " pprime(j)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<=jpsi; i++) ifile >> pprime[i];

  if(scan_until(ifile, " ffprime(j)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<=jpsi; i++) ifile >> ffprime[i];

  if(scan_until(ifile, " rcc(j,i)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<jpsi*itht; i++) ifile >> rcc[i];

  if(scan_until(ifile, " zcc(j,i)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<jpsi*itht; i++) ifile >> zcc[i];

  if(scan_until(ifile, " dpsidr(j,i)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<jpsi*itht; i++) ifile >> dpsidr[i];

  if(scan_until(ifile, " dpsidz(j,i)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<jpsi*itht; i++) ifile >> dpsidz[i];

  if(scan_until(ifile, " psimesh(j,i)") != FIO_SUCCESS) return FIO_FILE_ERROR;
  for(int i=0; i<jpsi*itht; i++) ifile >> psimesh[i];



  ifile.close();

  std::ofstream ofile;
  ofile.open("testout", std::ofstream::out | std::ofstream::trunc);
  if(!ofile){
    std::cerr << "ERROR!" << std::endl;
  }
  for(int i=0; i<jpsi*itht; i++)  {
    //    std::cerr << i;
    //    std::cout << rcc[i] << " " << zcc[i] << '\n';
    ofile << rcc[i]  << '\t' << zcc[i] << '\n';
    //    std::cerr<< jpsi*itht << " " << zcc << " " << ftor;
    //    ofile << "hi";
  }
  ofile.close();

  set_element_bounds();

  int e;
  double x[3] = { 1.6, 0.0, 0.3 };
  
  if(get_element(x, &e)==FIO_SUCCESS) {
    std::cerr << e << " " << rcc[e] << " " << zcc[e] << std::endl;
  }

  return FIO_SUCCESS;
}

int gato_source::close()
{
  if(psival) delete[] psival;
  if(pressure) delete[] pressure;
  if(ftor) delete[] ftor;
  if(pprime) delete[] pprime;
  if(ffprime) delete[] ffprime;
  if(rcc) delete[] rcc;
  if(zcc) delete[] zcc;
  if(r_bound) delete[] r_bound;
  if(t_bound) delete[] t_bound;
  if(psimesh) delete[] psimesh;
  if(dpsidr) delete[] dpsidr;
  if(dpsidz) delete[] dpsidz;

  psival = 0;
  pressure = 0;
  ftor = 0;
  pprime = 0;
  ffprime = 0;
  rcc = 0;
  zcc = 0;
  r_bound = 0;
  t_bound = 0;
  psimesh = 0;
  dpsidr= 0;
  dpsidz = 0;

  return FIO_SUCCESS;
}

int gato_source::set_element_bounds()
{
  r_bound = new double[jpsi*itht];
  t_bound = new double[jpsi*itht];

  int k=0;
  double r[2];
  double t[2];
  for(int i=0; i<itht; i++) {
    for(int j=0; j<jpsi; j++) {

      if(i==0) {
	t_bound[k] = 0.;
      }	else {
	int l = k-jpsi;
	if(l < 0) l+=jpsi*itht; 
	t[0] = atan2(zcc[l]-zma, rcc[l]-xma);
	t[1] = atan2(zcc[k]-zma, rcc[k]-xma);
	if(t[0] < 0) t[0] += 2.*M_PI;
	if(t[1] < 0) t[1] += 2.*M_PI;
	if(t[1] < t[0]) t[1] += 2.*M_PI;
	if(abs(t[1]-t[0]) > .1) 
	  std::cerr << "Large discrepancy in theta " 
		    << i << " " 
		    << t[0] << " " << t[1] << std::endl;
	t_bound[k] = (t[0]+t[1])/2.;
      }

      if(j==0) {
	r_bound[k] = 0.;
      } else {
	int l = k-1;
	r[0] = sqrt((zcc[l]-zma)*(zcc[l]-zma)+(rcc[l]-xma)*(rcc[l]-xma));
	r[1] = sqrt((zcc[k]-zma)*(zcc[k]-zma)+(rcc[k]-xma)*(rcc[l]-xma));
	/*
	if(r[1] < r[0]) std::cerr << "Error: r[0] > r[1] " << j << " "
				  << r[0] << " " << r[1] << " " 
				  << rcc[k] << " " << zcc[k] << std::endl;
	*/
	r_bound[k] = (r[0]+r[1])/2.;
      }

      k++;
    }
  }

  return FIO_SUCCESS;
}
