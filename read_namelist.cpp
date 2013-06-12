#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <deque>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define INT_TYPE    0
#define DOUBLE_TYPE 1
#define STRING_TYPE 2

#define MAX_STRING_LENGTH 256;

struct variable {
  std::string name;
  std::string description;
  int group;
  void* data;

  variable(std::string n, void* var, std::string desc, int g)
  { 
    name = n;
    description = desc;
    data = var;
    group = g;
  }
  int get_group() const { return group; }
  virtual int get_type() const = 0;
  virtual bool read(const std::string) = 0;
  virtual void printval() const = 0;
  virtual void print(bool print_desc) const
  {
    std::cout << std::setw(20) << std::left << name << " = ";
    printval();
    if(print_desc && description.length()) {
      size_t pos = 0;
      size_t oldpos = 0;
      do {
	pos = description.find_first_of("|", oldpos);
	if(pos==std::string::npos) pos=description.length();
	if(oldpos != 0) std::cout << std::setw(25) << "";
	std::cout << "\t! " << description.substr(oldpos, pos) << '\n';
	oldpos = pos+1;
      } while(pos != description.length());
    } else {
      std::cout << '\n';
    }
  }
};

struct string_variable : public variable {
  int size;
  char pad;

  string_variable(std::string n, char* var, const int s, const char* val, 
		  std::string desc, int g, const char p='\0')
    : variable(n, (void*)var, desc, g)
  { size = s; pad=p; strncpy(var, val, size); do_padding(); }
  virtual int get_type() const 
  { return STRING_TYPE; }
  virtual bool read(const std::string s)
  {
    strncpy((char*)data, s.c_str(), size);
    do_padding();
    return true;
  }
  virtual void printval() const
  { 
    int l = 0;
    if(pad != '\0') {
      for(int i=0; i<size; i++) {
	if(((char*)data)[i] != pad) l = i;
      }
      ((char*)data)[l+1] = '\0';
    }
    std::cout << (char*)data; 
    if(pad != '\0') {
      ((char*)data)[l+1] = pad;
    }
  }

private:
  void do_padding() {
    if(pad=='\0') return;
    for(int i=0; i<size; i++) {
      if(((char*)data)[i] == '\0') 
	((char*)data)[i] = pad;
    }
  }
};


struct int_variable : public variable {
  int_variable(std::string n, int* var, int val, std::string desc, int g)
    : variable(n, (void*)var, desc, g)
  { *((int*)var) = val;  }
  virtual int get_type() const 
  { return INT_TYPE; }
  virtual bool read(const std::string s)
  {
    *((int*)data) = atoi(s.c_str());
    return true;
  }
  virtual void printval() const
  { std::cout << *((int*)data); }
};

struct double_variable : public variable {
  double_variable(std::string n, double* var, double val, std::string desc, 
		  int g)
    : variable(n, (void*)var, desc, g)
  { *((double*)var) = val;  }
  virtual int get_type() const 
  { return DOUBLE_TYPE; }
  virtual bool read(const std::string s)
  { 
    *((double*)data) = atof(s.c_str());
    return true;
  }
  virtual void printval() const
  { std::cout << *((double*)data); }
};

class variable_list {
  typedef std::map<std::string, variable*> variable_map;
  typedef std::deque<std::string> group_deque;
  typedef std::map<int, int> type_map;
  variable_map variables;
  group_deque groups;
  type_map type_num;

public:
  ~variable_list() 
  {
    variable_map::iterator i = variables.begin();
    while(i != variables.end()) {
      delete(i->second);
      i++;
    }
  }

  int create_group(std::string n) {
    int g = groups.size();
    groups.push_back(n);
    return g;
  }

  bool create_variable(variable* v) {
    if(variables.find(v->name) != variables.end()) {
      std::cerr << "Warning: variable " << v->name << " already created."
		<< std::endl;
      delete v;
      return false;
    }
    if(v->get_group() < 0 || v->get_group() >= groups.size()) {
      std::cerr << "Warning: invalid group " << v->get_group()
		<< std::endl;
      delete v;
      return false;
    }

    variables[v->name] = v;

    if(type_num.find(v->get_type()) == type_num.end())
      type_num[v->get_type()] = 1;
    else
      type_num[v->get_type()]++;
    return true;
  }

  void print_group(int g, bool print_desc) {
    std::cout << "! " << groups[g] << '\n'
	      << "! ------------------------------\n";
    variable_map::iterator i = variables.begin();
    while(i != variables.end()) {
      if(i->second->get_group() == g) 
	i->second->print(print_desc);
      i++;
    }
  }

  void print_all_groups(bool print_desc) {
    for(int g=0; g<groups.size(); g++) {
      std::cout << '\n';
      print_group(g, print_desc);
    }
  }

  void print_all(bool print_desc) {
    variable_map::iterator i = variables.begin();
    while(i != variables.end()) {
      i->second->print(print_desc);
      i++;
    }   
  }

  int read_namelist(std::string filename) {
    int ierr = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0) {
      std::string op1, op2, op3;
      std::fstream file;
      char buffer[256];
      file.open(filename.c_str(), std::ios_base::in);
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

	  size_t pc = str.find_first_of("!");
	  if(pc < peq) continue;
	
	  size_t p1 = str.find_first_not_of(" \n\t=");
	  if(p1 == std::string::npos) continue;
	  size_t p2 = str.find_last_not_of(" \n\t=",peq);
	  if(p2 == std::string::npos) continue;

	  size_t p3 = str.find_first_not_of(" \n\t='\"",peq);
	  if(p3 == std::string::npos) continue;
	  size_t p4 = str.find_last_not_of(" \n\t='\"",pc);
	  if(p4 == std::string::npos) p4 = str.size();

	  if(pc < p4) continue;
	  
	  op1 = str.substr(p1,p2-p1+1);
	  op2 = str.substr(p3,p4-p3+1);
	  
	  variable_map::iterator i = variables.find(op1);

	  if(i == variables.end()) {
	    std::cout << "Warning: variable not recognized: " 
		      << "|" << op1 << "|" << std::endl;
	    ierr = 1;
	  } else {
	    i->second->read(op2);
	  }
	}
      }
	
      file.close();
    }

    int ierr_out;
    MPI_Allreduce(&ierr, &ierr_out, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
      
    if(ierr_out==2) return ierr_out;

    communicate();
    return ierr_out;
  }

  void communicate() {
    variable_map::iterator i;
    int rank, k;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Communicate doubles
    if(type_num[DOUBLE_TYPE] > 0) {
      double* doubles = new double[type_num[DOUBLE_TYPE]];
      if(rank==0) {
	i = variables.begin();
	k = 0;
	while(i != variables.end()) {
	  if(i->second->get_type() == DOUBLE_TYPE)
	    doubles[k++] = *((double*)(i->second->data));
	  i++;
	}
	if(k != type_num[DOUBLE_TYPE])
	  std::cerr << "Error: incorrect number of doubles (rank == 0)"
		    << std::endl;
      }
      
      MPI_Bcast(doubles, type_num[DOUBLE_TYPE], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      if(rank!=0) {
	i = variables.begin();
	k = 0;
	while(i != variables.end()) {
	  if(i->second->get_type() == DOUBLE_TYPE)
	    *((double*)(i->second->data)) = doubles[k++];
	  i++;
	}
	if(k != type_num[DOUBLE_TYPE])
	  std::cerr << "Error: incorrect number of doubles (rank != 0)"
		    << std::endl;
      }     
      delete[] doubles;
    }

    // Communicate integers
    if(type_num[INT_TYPE] > 0) {
      int* ints = new int[type_num[INT_TYPE]];
      if(rank==0) {
	i = variables.begin();
	k = 0;
	while(i != variables.end()) {
	  if(i->second->get_type() == INT_TYPE)
	    ints[k++] = *((int*)(i->second->data));
	  i++;
	}
	if(k != type_num[INT_TYPE])
	  std::cerr << "Error: incorrect number of ints (rank == 0)"
		    << std::endl;
      }
      
      MPI_Bcast(ints, type_num[INT_TYPE], MPI_INTEGER, 0, MPI_COMM_WORLD);
      
      if(rank!=0) {
	i = variables.begin();
	k = 0;
	while(i != variables.end()) {
	  if(i->second->get_type() == INT_TYPE)
	    *((int*)(i->second->data)) = ints[k++];
	  i++;
	}
	if(k != type_num[INT_TYPE])
	  std::cerr << "Error: incorrect number of ints (rank != 0)"
		    << std::endl;
      }
      delete[] ints;
    }

    // Communicate strings
    if(type_num[STRING_TYPE] > 0) {
      int len;

      i = variables.begin();

      k = 0;
      while(i != variables.end()) {
	if(i->second->get_type() == STRING_TYPE) {
	  MPI_Bcast(i->second->data, ((string_variable*)i->second)->size, 
	    MPI_CHAR, 0, MPI_COMM_WORLD);
	  k++;
	}
	i++;
      }
      if(k != type_num[STRING_TYPE])
	std::cerr << "Error: incorrect number of strings"
		  << " " << k << " vs " << type_num[STRING_TYPE] 
		  << std::endl;
    }
  }
};

static variable_list variables;

extern "C" void create_group_(const char* name, int* id)
{
  *id = variables.create_group(name);
}

extern "C" void add_variable_int_(const char* name, int* v, int* def, 
				  const char* description, const int* group)
{
  variables.create_variable(new int_variable(name, v, *def, 
					     description, *group));
}
 
extern "C" void add_variable_double_(const char* name, double* v, double* def, 
				     const char* description, const int* group)
{
  variables.create_variable(new double_variable(name, v, *def, 
						description, *group));
}

extern "C" void add_variable_string_(const char* name, char* v, const int* sz,
				     const char* def, const char* description,
				     const int* group, const char* pad)
{
  variables.create_variable(new string_variable(name, v, *sz, def,
						description, *group, *pad));
}

// ierr = 0: success
// ierr = 1: unknown variable
// ierr = 2: can't open file
extern "C" void read_namelist_(const char* filename, int* ierr)
{
  *ierr = variables.read_namelist(filename);
}

// style = 0: ungrouped, without descriptions
// style = 1: ungrouped, with desctiptions
// style = 2: grouped, without descriptions
// style = 3: grouped, with descriptions
extern "C" void print_variables_(int* style)
{
  switch(*style) {
  case(1): variables.print_all(true); break;
  case(2): variables.print_all_groups(false); break;
  case(3): variables.print_all_groups(true); break;
  default: variables.print_all(false); break;
  }
}

/*
int main(int argc, char** argv) {

  int ierr;

  int transport_group = 2;
  int model_group = 1;
  int misc_group = 0;

  double kappat, amu, nowhere, etar;
  int numvar, gyro, eqsubtract;
  double def = 3.;
  int idef = 5;
  
  MPI_Init(&argc, &argv);

  create_group_(&transport_group, "Transport Parameters");
  create_group_(&model_group, "Model Options");
  create_group_(&misc_group, "Miscellaneous");

  add_variable_double_("etar", &etar, &def, 
		       "Resistivity", &transport_group);
  add_variable_double_("kappat", &kappat, &def,
		       "Thermal conductivity", &transport_group);
  add_variable_double_("amu", &amu, &def, 
		       "Viscosity", &transport_group);
  add_variable_double_("nowhere", &nowhere, &def, 
		       "Nothing", &misc_group);
  add_variable_int_("numvar", &numvar, &idef, 
		    "Model", &model_group);
  add_variable_int_("gyro", &gyro, &idef, 
		    "Gyroviscosity", &model_group);
  add_variable_int_("eqsubtract", &eqsubtract, &idef, 
		    "1 = subtract equilibrium quantites", &transport_group);

  read_namelist_("C1input", &ierr);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 1) print_variables_(3);

  MPI_Finalize();

  return 0;
}
*/
