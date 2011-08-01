#ifndef M3DC1_FILE_H
#define M3DC1_FILE_H

#include "m3dc1_field.h"
#include <map>
#include <string>
#include <vector>

typedef std::vector<double> m3dc1_scalar_list;

class m3dc1_file {
  hid_t file;

  typedef std::map<int, m3dc1_timeslice> m3dc1_timeslice_map;
  m3dc1_timeslice_map timeslice_map;

  typedef std::map<std::string, m3dc1_scalar_list> m3dc1_scalar_map;
  m3dc1_scalar_map scalar_map;

  bool time_name(int time, char* name);
  hid_t open_timeslice(int time);
  m3dc1_mesh* read_mesh(int time);
  m3dc1_timeslice* read_timeslice(const int time);

 public:
  m3dc1_file();
  ~m3dc1_file();

  bool open(const char* filename);
  bool close();

  m3dc1_scalar_list* read_scalar(const char*);
  m3dc1_field* read_field(const char* name, const int time);
  bool read_parameter(const char*, int*); 
  bool read_parameter(const char*, double*); 

};

#endif
