#include "m3dc1_file.h"
#include <math.h>

#include <iostream>


m3dc1_file::m3dc1_file()
{
  file = -1;
}

m3dc1_file::~m3dc1_file()
{
  close();
}

bool m3dc1_file::open(const char* filename)
{
  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file < 0) {
    std::cerr << "Error: couldn't open file " << filename << std::endl;
    return false;
  }

  return true; 
}

bool m3dc1_file::close()
{
  if(file < 0) return 0;
  if(H5Fclose(file) < 0) return false;
  file = -1;
  return true;
}

bool m3dc1_file::time_name(int time, char* name)
{
  if(time < 0)
    sprintf(name, "equilibrium");
  else
    sprintf(name, "time_%03d", time);
  return true;
}

hid_t m3dc1_file::open_timeslice(int t)
{
  char name[256];

  time_name(t, name);
  hid_t time_group = H5Gopen(file, name, H5P_DEFAULT);

  return time_group;
}

m3dc1_timeslice* m3dc1_file::read_timeslice(const int t)
{
  char name[256];
  m3dc1_timeslice_map::iterator i;

  i = timeslice_map.find(t);
  if(i != timeslice_map.end()) {
    std::cerr << "Timeslice " << t << " has already been read."  << std::endl;
    return &i->second;
  }

  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening time slice " << name << std::endl;
    return 0;
  }

  std::pair<m3dc1_timeslice_map::iterator,bool> tsi;
  tsi = timeslice_map.insert(
    m3dc1_timeslice_map::value_type(t, m3dc1_timeslice()));

  if(!tsi.second) {
    std::cerr << "Error inserting timeslice" << std::endl;
    return 0;
  }
  m3dc1_timeslice* ts = &tsi.first->second;

  if(t >= 0) read_parameter("ntor", &ts->ntor);
  else ts->ntor = 0;
  std::cerr << "ntor = " << ts->ntor << std::endl;

  read_parameter("3d", &ts->is_3d);
  std::cerr << "3d = " << ts->is_3d << std::endl;

  hid_t attr_id = H5Aopen(time_group, "time", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_DOUBLE, (void*)(&ts->time));
  H5Aclose(attr_id);
  
  ts->mesh = read_mesh(t);
  if(!ts->mesh) { 
    std::cerr << "Error reading mesh" << std::endl;
    return 0;
  }

  H5Gclose(time_group);

  return ts;
}

m3dc1_mesh* m3dc1_file::read_mesh(const int t)
{
  int nelms;

  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening time group" << std::endl;
    return 0;
  }

  hid_t mesh_group = H5Gopen(time_group, "mesh", H5P_DEFAULT);
  if(mesh_group < 0) {
    std::cerr << "Error opening mesh" << std::endl;
    return 0;
  }

  hid_t nelms_attr = H5Aopen(mesh_group, "nelms", H5P_DEFAULT);
  H5Aread(nelms_attr, H5T_NATIVE_INT, &nelms);
  H5Aclose(nelms_attr);
  if(nelms < 1) {
    std::cerr << "Error: too few elements" << std::endl;
    return 0;
  }

  int is_3d = 0;
  read_parameter("3d", &is_3d);
  std::cerr << "is_3d = " << is_3d << std::endl;

  int nfields;

  if(is_3d) nfields = 9;
  else nfields = 7;

  double* data = new double[nfields*nelms];

  hid_t mesh_dataset = H5Dopen(mesh_group, "elements", H5P_DEFAULT);
  H5Dread(mesh_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	  (void*)data);
  H5Dclose(mesh_dataset);

  m3dc1_mesh* mesh;
  if(is_3d) mesh = new m3dc1_3d_mesh(nelms);
  else mesh = new m3dc1_mesh(nelms);

  // set default mesh memory depth for searching
  mesh->set_memory_depth(3);

  for(int i=0; i<nelms; i++) {
    mesh->a[i]     =      data[i*nfields  ];
    mesh->b[i]     =      data[i*nfields+1];
    mesh->c[i]     =      data[i*nfields+2];
    mesh->co[i]    =  cos(data[i*nfields+3]);
    mesh->sn[i]    =  sin(data[i*nfields+3]);
    mesh->x[i]     =      data[i*nfields+4];
    mesh->z[i]     =      data[i*nfields+5];
    mesh->bound[i] = (int)data[i*nfields+6];
    if(is_3d) {
      ((m3dc1_3d_mesh*)mesh)->d[i]   =      data[i*nfields+7];
      ((m3dc1_3d_mesh*)mesh)->phi[i] =      data[i*nfields+8];
    }
  }
  delete[] data;
  
  H5Gclose(mesh_group);
  H5Gclose(time_group);

  return mesh;
}


bool m3dc1_file::read_parameter(const char* n, double* d)
{
  hid_t attr_id = H5Aopen(file, n, H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_DOUBLE, d);
  H5Aclose(attr_id);

  return true;  
}

bool m3dc1_file::read_parameter(const char* n, int* d)
{
  hid_t attr_id = H5Aopen(file, n, H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, d);
  H5Aclose(attr_id);

  return true;  
}

m3dc1_scalar_list* m3dc1_file::read_scalar(const char* name)
{
  m3dc1_scalar_map::iterator i;
  i = scalar_map.find(name);
  if(i != scalar_map.end()) {
    std::cerr << "Scalar " << name << " has already been read" << std::endl;
    return &i->second;
  }

  // Open fields group
  hid_t scalar_group = H5Gopen(file, "scalars", H5P_DEFAULT);
  if(scalar_group < 0) {
    std::cerr << "Error opening scalars group" << std::endl;
    return 0;
  }

  // Open field dataset
  hid_t scalar_dataset = H5Dopen(scalar_group, name, H5P_DEFAULT);
  if(scalar_dataset < 0) {
    std::cerr << "Error opening scalar " << name << std::endl;
    H5Gclose(scalar_group);
    return 0;
  }

  hid_t scalar_space = H5Dget_space(scalar_dataset);
  hssize_t n = H5Sget_simple_extent_npoints(scalar_space);
  H5Sclose(scalar_space);

  if(n==0) {
    std::cerr << "Error: no elements in dataset " << name << std::endl;
    return 0;
  }

  m3dc1_scalar_list* s = &(scalar_map[name] = m3dc1_scalar_list(n));

  // Read Field Data
  H5Dread(scalar_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  	  (void*)(&((*s)[0])));

  H5Dclose(scalar_dataset);
  H5Gclose(scalar_group);

  return s;
}

m3dc1_field* m3dc1_file::read_field(const char* n, const int t)
{
  hid_t fields_group, field_dataset, fieldi_dataset;
  std::map<std::string, m3dc1_field*>::iterator i;
  std::string name = n;
  std::string name_i = name + "_i";

  m3dc1_timeslice* ts = read_timeslice(t);
  if(!ts) {
    std::cerr << "Error reading timeslice" << std::endl;
    return 0;
  }

  // Check if field has already been read
  i = ts->field_map.find(name);
  if(i != ts->field_map.end()) {
    std::cerr << "Field " << name << " has already been read" << std::endl;
    return i->second;
  }

  bool is_3d = (ts->is_3d==1);
  bool is_complex = (ts->ntor!=0 && !is_3d);

  // open time group
  hid_t time_group = open_timeslice(t);
  if(time_group < 0) {
    std::cerr << "Error opening timeslice" << std::endl;
    return 0;
  }

  // Open fields group
  fields_group = H5Gopen(time_group, "fields", H5P_DEFAULT);
  if(fields_group < 0) {
    std::cerr << "Error opening fields group" << std::endl;
    return 0;
  }

  // Open field dataset
  field_dataset = H5Dopen(fields_group, name.data(), H5P_DEFAULT);
  if(field_dataset < 0) {
    std::cerr << "Error opening field " << name << std::endl;
    H5Gclose(fields_group);
    return 0;
  }

  // For complex fields, open imaginary field dataset
  if(is_complex) {
    fieldi_dataset = H5Dopen(fields_group, name_i.data(), H5P_DEFAULT);
    if(fieldi_dataset < 0) {
      std::cerr << "Error opening field " << name_i << std::endl;
      is_complex = false;
    }
  }

  // Create new field
  m3dc1_field* field;
  if(is_3d)
    field = new m3dc1_3d_field(ts->mesh);
  else if(is_complex)
    field = new m3dc1_complex_field(ts->mesh,ts->ntor);
  else
    field = new m3dc1_field(ts->mesh);

  // Read Field Data
  H5Dread(field_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  	  (void*)field->data);
  H5Dclose(field_dataset);

  if(is_complex) {
    H5Dread(fieldi_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	    (void*)((m3dc1_complex_field*)field)->data_i);
    H5Dclose(fieldi_dataset);
  }
  
  H5Gclose(fields_group);
  H5Gclose(time_group);

  field->time = ts->time;

  // Add field to list of fields that have been read
  ts->field_map.insert(
    m3dc1_timeslice::m3dc1_field_map::value_type(name, field));

  return field;
}

