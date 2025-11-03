/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifdef M3DC1_PETSC
#ifndef M3DC1_SOLVER_H
#define M3DC1_SOLVER_H
// #include "superlu_ddefs.h" // gridinfo_t
#include "apf.h"
#include "apfMesh2.h"
#include "apfNumbering.h"
#include "m3dc1_scorec.h"
#include "petscksp.h"
#include <vector>

int copyField2PetscVec(FieldID field, Vec petscVec, int scalar_type);
int copyPetscVec2Field(Vec &petscVec, FieldID field, int scalar_type);
void printMemStat();
PetscErrorCode MyKSPMonitor(KSP, PetscInt, PetscReal, void *);
// NOTE: all field realted interaction is done through m3dc1 api rather than apf
class m3dc1_matrix {
public:
  m3dc1_matrix(int i, int s, FieldID field);
  virtual ~m3dc1_matrix();
  virtual int initialize() = 0; // create a matrix and solver object
  int destroy();                // delete a matrix and solver object
  int set_value(int row, int col, int operation, double real_val,
                double imag_val); // insertion/addition with global numbering
  // values use column-wise, size * size block
  int add_values(int rsize, PetscInt *rows, int csize, PetscInt *columns,
                 double *values);
  int get_values(std::vector<PetscInt> &rows, std::vector<int> &n_columns,
                 std::vector<PetscInt> &columns, std::vector<double> &values);

  int get_status() { return mat_status; }
  int get_scalar_type() { return scalar_type; }
  int get_fieldOrdering() { return fieldOrdering; }
  int write(const char *file_name);
  virtual int reset_values() = 0;
  virtual int update_values() = 0;
  virtual int get_type() const = 0;
  virtual int assemble() = 0;
  virtual int setupMat() = 0;
  virtual int preAllocate() = 0;
  void printInfo();
  PetscInt mymatrix_id;
  // PETSc data structures
  Mat _A;
  int fieldOrdering; // the field that provide numbering
  apf::Mesh2 *mesh;

protected:
  int setupSeqMat();
  int setupParaMat();
  int preAllocateSeqMat();
  int preAllocateParaMat();
  int id;
  // 0 is for real, 1 is for complex
  int scalar_type;
  // 0: A is not allocated/used
  // 1: A is allocated/used but not solved
  // 2: A is solved
  int mat_status;
};

class matrix_mult : public m3dc1_matrix {
public:
  matrix_mult(int i, int s, FieldID field);
  virtual int initialize();
  void set_mat_local(bool flag) { localMat = flag; }
  int is_mat_local() { return localMat; }
  int multiply(FieldID in_field, FieldID out_field);
  int reset_values() {
    MatZeroEntries(_A);
    mat_status = M3DC1_NOT_FIXED;
    return M3DC1_SUCCESS;
  };
  int update_values() { MatZeroEntries(_A); mat_status=M3DC1_NOT_FIXED; return M3DC1_SUCCESS; };
  virtual int get_type() const { return 0; } // M3DC1_MULTIPLY; }
  virtual int assemble();
  virtual int setupMat();
  virtual int preAllocate();

private:
  bool localMat;
};

class matrix_solve : public m3dc1_matrix {
public:
  matrix_solve(int i, int s, FieldID fieldOrdering);
  virtual int initialize();
  virtual ~matrix_solve();
  int solve(FieldID field_id);
  int solve_with_guess(FieldID field_id, FieldID xVec_guess);
  int set_bc(int row);
  int set_row(int row, int numVals, int *colums, double *vals);
  int add_blockvalues(int rbsize, PetscInt *rows, int cbsize, PetscInt *columns,
                      double *values);
  int reset_values();
  int update_values();
  virtual int get_type() const { return 1; }
  virtual int assemble();
  virtual int setupMat();
  virtual int preAllocate();
  PetscInt its;

private:
  int setUpRemoteAStruct();
  int setKspType();
  int _kspSet;
  KSP _ksp;
  Mat remoteA;

  // block mg in toroidal direction
  int _BgmgSet;         // only for mymatrix_id=5 or 17, the hard ones
  PetscInt mg_nlevels; // default = 2
  // PC *pc;
  Mat *mg_interp_mat;
  KSP *mg_level_ksp;
  PC *mg_level_pc;
  int setBgmgType();
  int mapping(int, int, int, int, int, int, int, int, int *, int *, int *);

  // remoteA related data
  std::set<int> *remotePidOwned;
  std::map<int, std::map<int, int>> *remoteNodeRow; // <pid, <locnode>, numAdj >
  std::map<int, int> *remoteNodeRowSize;
};

class m3dc1_solver {
public:
  // functions
  m3dc1_solver() : assembleOption(0) {
    matrix_container = new std::map<int, m3dc1_matrix *>;
    PetscMemorySetGetMaximumUsage();
  }
  ~m3dc1_solver();
  static m3dc1_solver *instance();
  m3dc1_matrix *get_matrix(int matrix_id);
  void add_matrix(int matrix_id, m3dc1_matrix *);
  // data
  std::map<int, m3dc1_matrix *> *matrix_container;
  int assembleOption; // 0 use scorec; 1 use petsc
private:
  static m3dc1_solver *_instance;
};

#endif
#endif // #ifndef M3DC1_MESHGEN
