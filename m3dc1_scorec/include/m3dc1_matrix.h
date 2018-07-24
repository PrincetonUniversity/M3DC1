/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifdef M3DC1_PETSC
#ifndef M3DC1_SOLVER_H
#define M3DC1_SOLVER_H
//#include "superlu_ddefs.h" // gridinfo_t
#include "apf.h"
#include "petscksp.h"
#include "apfNumbering.h"
#include "m3dc1_scorec.h"
#include <vector>

int copyField2PetscVec(FieldID field, Vec& petscVec, int scalar_type);
int copyPetscVec2Field(Vec& petscVec, FieldID field, int scalar_type);
void printMemStat();
// NOTE: all field realted interaction is done through m3dc1 api rather than apf
class m3dc1_matrix
{
public:
  m3dc1_matrix(int i, int s, FieldID field);
  virtual ~m3dc1_matrix();
  virtual int initialize()=0; // create a matrix and solver object
  int destroy(); // delete a matrix and solver object
  int set_value(int row, int col, int operation, double real_val, double imag_val); //insertion/addition with global numbering
  // values use column-wise, size * size block
  int add_values(int rsize, int * rows, int csize, int * columns, double* values);
  int get_values(std::vector<int>& rows, std::vector<int>& n_columns, std::vector<int>& columns, std::vector<double>& values);
  void set_status(int s) {mat_status=s;}
  int get_status() {return mat_status;}
  int get_scalar_type() { return scalar_type; }
  int get_fieldOrdering() { return fieldOrdering;}
  int write( const char* file_name);
  virtual void reset_values() = 0;
  virtual int get_type() const = 0;
  virtual int assemble() = 0;
  virtual int setupMat() =0;
  virtual int preAllocate() =0;
  int printInfo();
  // PETSc data structures
  Mat* A;
protected:
  int setupSeqMat();
  int setupParaMat();
  int preAllocateSeqMat();
  int preAllocateParaMat();
  int id;
  int scalar_type;
  int mat_status; 
  int fieldOrdering; // the field that provide numbering
};

class matrix_mult: public m3dc1_matrix
{
public:
  matrix_mult(int i, int s, FieldID field): m3dc1_matrix(i,s,field), localMat(1) { initialize();}
  virtual int initialize();
  void set_mat_local(bool flag) {localMat=flag;}
  int is_mat_local() {return localMat;}
  int multiply(FieldID in_field, FieldID out_field);
  void reset_values() { MatZeroEntries(*A);   set_status(M3DC1_NOT_FIXED); };
  virtual int get_type() const { return 0; } //M3DC1_MULTIPLY; }
  virtual int assemble();
  virtual int setupMat();
  virtual int preAllocate();
private:
  bool localMat;
};

class matrix_solve: public m3dc1_matrix
{
public:
  matrix_solve(int i, int s,  FieldID fieldOrdering);
  virtual int initialize();
  virtual ~matrix_solve();
  int solve(FieldID field_id);
  int set_bc( int row);
  int set_row( int row, int numVals, int* colums, double * vals);
  int add_blockvalues( int rbsize, int * rows, int cbsize, int * columns, double* values);
  void reset_values();
  virtual int get_type() const {return 1; }
  virtual int assemble(); 
  virtual int setupMat();
  virtual int preAllocate();
  int iterNum;
private:  
  int setUpRemoteAStruct();
  int setKspType(Vec *b, Vec *x);
  int kspSet;
  KSP* ksp; 
  Mat remoteA;
  std::set<int> remotePidOwned;
  std::map<int, std::map<int, int> > remoteNodeRow; // <pid, <locnode>, numAdj >
  std::map<int, int> remoteNodeRowSize;
};

class m3dc1_solver
{
public:
// functions
  m3dc1_solver():assembleOption(0) {matrix_container = new std::map<int, m3dc1_matrix*>;PetscMemorySetGetMaximumUsage();}
  ~m3dc1_solver();
  static m3dc1_solver* instance();
  m3dc1_matrix* get_matrix(int matrix_id);
  void add_matrix(int matrix_id, m3dc1_matrix*);
// data
  std::map<int, m3dc1_matrix*>* matrix_container;
  int assembleOption; // 0 use scorec; 1 use petsc
private:
  static m3dc1_solver* _instance;
};

#endif
#endif //#ifndef M3DC1_MESHGEN
