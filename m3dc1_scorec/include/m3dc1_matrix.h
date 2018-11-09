/******************************************************************************
  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifdef M3DC1_PETSC
#ifndef M3DC1_SOLVER_H
#define M3DC1_SOLVER_H
#include "m3dc1_scorec.h"
#include "m3dc1_field.h"
#include <apf.h>
#include <apfNumbering.h>
#include <petscksp.h>
#include <cassert>
#include <vector>
#include <ostream>
void las_init(int * argc, char ** argv[], MPI_Comm cm);
class m3dc1_matrix;
class m3dc1_mesh;

void double2petsc(const int n, double * in, PetscScalar *& out);
void petsc2double(const int n, PetscScalar * in, double *& out);

void field2Vec(m3dc1_field * fld, Vec v);
void vec2Field(Vec v, m3dc1_field * fld);

// does the allocated matrix locally own row rw
bool ownRow(Mat A, int rw);
// will the unallocated matrix locally own row rw
//bool willOwnRow(Mat A, int lcl_rws, int rw);
// get the rank of the owner of row rw for allocated matrices
int getRowOwner(Mat A, int rw);
// get the rank of the owner of row rw for unallocated matrices
//int calcRowOwner(Mat A, int lcl_rws, int rw);

void get_block_ids(m3dc1_field * fld, bool lcl, int ent_dim, int eid, int * blk_ids);


// the mat has the field implicitly, should make it explicit
void insert_element_blocks(m3dc1_matrix * mat,
                           int ent_dim,
                           int eid,
                           double * vals);
void insert_node_blocks(m3dc1_matrix * mat,
                        int ent_dim,
                        int eid,
                        int nd1,
                        int nd2,
                        double * vals);
void printMemStat();
class m3dc1_matrix
{
public:
  m3dc1_matrix(int i, int s, m3dc1_mesh * msh,  m3dc1_field * f, MPI_Comm cm);
  ~m3dc1_matrix();
  void multiply(m3dc1_field * in, m3dc1_field * out);
  void solve(m3dc1_field * lhs);
  bool willOwn(int rw);
  int  whoOwns(int rw);
  int  calcFirstRow();
  int  calcLastRowP1();
  bool is_fixed() { return fixed; }
  int get_scalar_type() { return scalar_type; }
  m3dc1_field * get_field() { return fld; }
  int get_block_size() { return blk_sz; }
  int is_parallel() { return is_par; }
  MPI_Comm getComm() { return cm; }
  void zero();
  void zero_rows(int rsize, int * rows);
  void add_values(int rsize, int * rows, int csize, int * cols, double * vals);
  void set_values(int rsize, int * rows, int csize, int * cols, double * vals);
  void get_values(std::vector<int>& rows, std::vector<int>& n_columns, std::vector<int>& columns, std::vector<double>& values);
  void add_blocks(int blk_rw_cnt, int * blk_rws, int blk_col_cnt, int * blk_cls, double * vals);
  int add_blockvalues( int rbsize, int * rows, int cbsize, int * columns, double* values);
  void fix();
  int solver_iteration_count();
  Vec * getRHS() { return &b; }
  Vec * getLHS() { return &x; }
  void printAssemblyInfo();
  void printNNZStats();
  void write(const char * fn);
protected:
  MPI_Comm cm;
  Mat A;
  Vec x;
  Vec b;
  KSP ksp;
  int id;
  int scalar_type;
  bool fixed;
  m3dc1_field * fld;
  int is_par;
  int blk_sz;
  int * ownership;
  int low_row;
  int hgh_row;
  int non_lcl_nnz_cnt;
  int lcl_nnz_cnt;
};

class m3dc1_solver
{
public:
  static m3dc1_solver * instance();
  m3dc1_matrix * get_matrix(int mid)
  {
    auto mit = matrix_container.find(mid);
    assert(mit != matrix_container.end() && "[M3D-C1 Error] Matrix with specified ID does not exist.");
    return mit->second;
  }
  void add_matrix(int mid, m3dc1_matrix * mat)
  {
    auto loc_new = matrix_container.insert(std::make_pair(mid,mat));
    assert(loc_new.second && "[M3D-C1 Error] Matrix with specified ID already exists!");
  }
  void destroy_matrix(int mid)
  {
    auto mit = matrix_container.find(mid);
    assert(mit != matrix_container.end() && "[M3D-C1 Error] Matrix with specified ID does not exist.");
    matrix_container.erase(mit);
  }
  void destroy_matrices()
  {
    matrix_container.clear();
  }
private:
  std::map<int, m3dc1_matrix*> matrix_container;
  static m3dc1_solver * _instance;
};
#endif
#endif //#ifndef M3DC1_MESHGEN
