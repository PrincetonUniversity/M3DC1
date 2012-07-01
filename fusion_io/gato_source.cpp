#include "fusion_io.h"
#include <iostream>
#include <iomanip>

gato_source::gato_source()
  : psival(0), pressure(0), ftor(0), pprime(0), ffprime(0)
{
}

gato_source::~gato_source()
{
  close();
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

  psival = new double[jpsi+1];
  pressure = new double[jpsi+1];
  ftor = new double[jpsi+1];
  pprime = new double[jpsi+1];
  ffprime = new double[jpsi+1];
  rcc = new double[jpsi*itht];
  zcc = new double[jpsi*itht];

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

  ifile.close();

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

  psival = 0;
  pressure = 0;
  ftor = 0;
  pprime = 0;
  ffprime = 0;
  rcc = 0;
  zcc = 0;

  return FIO_SUCCESS;
}

int gato_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  return FIO_SUCCESS;
}

int gato_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();

  return FIO_SUCCESS;
}

int gato_source::get_field(const field_type t,fio_field** f,
			    const fio_option_list* opt)
{
  *f = 0;

  switch(t) {

  default:
    return FIO_UNSUPPORTED;
  };

  return FIO_SUCCESS;
}
