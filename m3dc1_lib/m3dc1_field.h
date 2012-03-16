#ifndef M3DC1_H
#define M3DC1_H

#include <hdf5.h>
#include <map>
#include <vector>
#include <string>

#include "m3dc1_mesh.h"

class m3dc1_field {
 public:
  double* data;

  static const int nbasis = 20;
  static const int mi[nbasis];
  static const int ni[nbasis];

  enum m3dc1_get_op {
    GET_VAL    = 1,
    GET_DVAL   = 2,
    GET_DDVAL  = 4,
    GET_PVAL   = 8,
    GET_PPVAL  = 16
  };

  enum m3dc1_op {
    OP_1     = 0,
    OP_DR    = 1,
    OP_DZ    = 2,
    OP_DRR   = 3,
    OP_DRZ   = 4,
    OP_DZZ   = 5,
    OP_DP    = 6,
    OP_DRP   = 7,
    OP_DZP   = 8,
    OP_DRRP  = 9,
    OP_DRZP  = 10,
    OP_DZZP  = 11,
    OP_DPP   = 12,
    OP_DRPP  = 13,
    OP_DZPP  = 14,
    OP_DRRPP = 15,
    OP_DRZPP = 16,
    OP_DZZPP = 17,
    OP_NUM   = 18
  };

 public:
  double time;
  m3dc1_mesh* mesh;

  m3dc1_field() { mesh=0; data=0; }
  m3dc1_field(m3dc1_mesh* m);
  virtual ~m3dc1_field();

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element=0);
};

class m3dc1_complex_field : public m3dc1_field {
 public:
  double* data_i;
  int ntor;

 public:
  m3dc1_complex_field(m3dc1_mesh* m, int ntor);
  virtual ~m3dc1_complex_field();

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element=0);
};

class m3dc1_3d_field : public m3dc1_field {
 public:
  static const int tbasis = 4;
  static const int li[tbasis];

 public:
  m3dc1_3d_field(m3dc1_mesh* m);

  virtual bool eval(const double r, const double phi, const double z, 
		    const m3dc1_field::m3dc1_get_op op, double* val, 
		    int* element=0);
};

class m3dc1_compound_field : public m3dc1_field {
 public:
  typedef std::vector<m3dc1_field*> subfield_list;
  subfield_list subfield;

 public:
  virtual bool eval(const double r, const double phi, const double z,
		    const m3dc1_field::m3dc1_get_op op, double* val,
		    int* element=0);
};

class m3dc1_timeslice {
 public:
  typedef std::map<std::string, m3dc1_field*> m3dc1_field_map;
  m3dc1_field_map field_map;

  int ntor;
  int is_3d;

  double time;
  m3dc1_mesh* mesh;

  m3dc1_timeslice();
  ~m3dc1_timeslice();

  bool get_field(const char*, m3dc1_field**) const;
};

#endif
