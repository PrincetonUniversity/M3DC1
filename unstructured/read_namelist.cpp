#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"

typedef std::map<std::string, double*> double_variable_list;
typedef std::map<std::string, int*> int_variable_list;

static double_variable_list double_variables;
static int_variable_list int_variables;

static void communicate_ints()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int n = int_variables.size();
  int* vals = new int[n];

  if(rank==0) {
    int j=0;
    int_variable_list::iterator i = int_variables.begin();
    while(i != int_variables.end()) {
      vals[j] = *(i->second);
      i++;
      j++; 
    }
  }
  
  MPI_Bcast(vals, n, MPI_INTEGER, 0, MPI_COMM_WORLD);
  
  if(rank != 0) {
    int j=0;
    int_variable_list::iterator i = int_variables.begin();
    while(i != int_variables.end()) {
      *(i->second) = vals[j];
      i++;
      j++; 
    }
  }
}
  
static void communicate_doubles()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int n = double_variables.size();
  double* vals = new double[n];
  
  if(rank==0) {
    int j=0;
    double_variable_list::iterator i = double_variables.begin();
    while(i != double_variables.end()) {
      vals[j] = *(i->second);
      i++;
      j++; 
    }
  }
  
  MPI_Bcast(vals, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(rank != 0) {
    int j=0;
    double_variable_list::iterator i = double_variables.begin();
    while(i != double_variables.end()) {
      *(i->second) = vals[j];
      i++;
      j++; 
    }
  }
  
  delete[] vals;
}

extern "C" {
  void add_variable_int_(const char* name, int* v, int* def)
  {
    std::string s(name);
    *v = *def;
    int_variables[s] = v;
  }
 
  void add_variable_double_(const char* name, double* v, double* def)
  {
    std::string s(name);
    *v = *def;
    double_variables[s] = v;
  }

  // ierr = 0: success
  // ierr = 1: unknown variable
  // ierr = 2: can't open file
  void read_namelist_(const char* filename, int* ierr_out)
  {
    int ierr = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0) {
      std::string op1, op2, op3;
      std::fstream file;
      char buffer[256];
      file.open(filename, std::ios_base::in);
      if(file.fail()) {
	std::cout << "Error opening file " << filename << " for input" 
		  << std::endl;
	ierr = 2;
	
      } else {
	std::cout << "Reading input file " << filename << "..." << std::endl;

	while(!file.eof()) {
	  file.getline(buffer, 255);
	  std::string str(buffer);
	  
	  size_t peq = str.find_first_of("=");
	  if(peq==std::string::npos) continue;
	
	  size_t p1 = str.find_first_not_of(" \n\t=");
	  if(p1 == std::string::npos) continue;
	  size_t p2 = str.find_first_of(" \n\t=",p1);
	  if(p2 == std::string::npos) continue;

	  size_t p3 = str.find_first_not_of(" \n\t=",peq);
	  if(p3 == std::string::npos) continue;
	  size_t p4 = str.find_first_of(" \n\t=!",p3);
	  if(p4 == std::string::npos) p4 = str.size();

	  // remove comments
	  size_t pc = str.find_first_of("!");
	  if(pc < p4) continue;
	  
	  op1 = str.substr(p1,p2-p1);
	  op2 = str.substr(p3,p4-p3);
	  
	  double_variable_list::iterator di = double_variables.find(op1);
	  if(di != double_variables.end()) {
	    *(di->second) = atof(op2.c_str());
	    continue;
	  }

	  int_variable_list::iterator ii = int_variables.find(op1);
	  if(ii != int_variables.end()) {
	    *(ii->second) = atoi(op2.c_str());
	    continue;
	  }

	  std::cout << "Warning: variable not recognized: " 
		    << "|" << op1 << "|" << std::endl;
	  ierr = 1;
	}
      }
	
      file.close();
    }

    MPI_Allreduce(&ierr, ierr_out, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
      
    if(*ierr_out==2) return;

    communicate_doubles();
    communicate_ints();
  }

  void print_namelist_() 
  {
    int_variable_list::iterator ii = int_variables.begin();
    while(ii != int_variables.end()) {
      std::cout << ii->first << " = " << *(ii->second) << '\n';
      ii++;
    }

    double_variable_list::iterator di = double_variables.begin();
    while(di != double_variables.end()) {
      std::cout << di->first << " = " << *(di->second) << '\n';
      di++;
    }
  }
}

/*
int main(int argc, char** argv) {
  double kappat, amu, nowhere, etar;
  int numvar, gyro, eqsubtract;
  double def = 3.;
  int idef = 5;

  MPI_Init(&argc, &argv);

  add_variable_double_("etar", &etar, &def);
  add_variable_double_("kappat", &kappat, &def);
  add_variable_double_("amu", &amu, &def);
  add_variable_double_("nowhere", &nowhere, &def);
  add_variable_int_("numvar", &numvar, &idef);
  add_variable_int_("gyro", &gyro, &idef);
  add_variable_int_("eqsubtract", &eqsubtract, &idef);

  read_namelist_("C1input");

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 1) print_variables_();

  MPI_Finalize();

  return 0;
}
*/


