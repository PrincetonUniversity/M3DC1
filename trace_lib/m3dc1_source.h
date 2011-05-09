#ifndef TRACE_M3DC1_SOURCE_H
#define TRACE_M3DC1_SOURCE_H

#include "trace_source.h"
#include "m3dc1_file.h"

#include <string>

class m3dc1_source : public trace_field_source {
  static const int memory_depth = 2;
  int last_elm;
  int** next_elm;

  int hits, misses, evals, lost_shifts;

  m3dc1_file file;
  m3dc1_field *psi, *f, *g;

  double bzero, rzero;
  double R_axis, Z_axis;

 public:
  std::string filename;
  int time;
  double factor;

  m3dc1_source();
  m3dc1_source(std::string f, int t);
  ~m3dc1_source();

  bool load();
  bool eval(const double r, const double phi, const double z,
	    double* b_r, double* b_phi, double* b_z);
  virtual bool center(double* r0, double* z0) const;
  virtual bool extent(double* r0, double* r1, double* z0, double* z1) const;
};

#endif
