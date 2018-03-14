/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef M3DC1_SCOREC_HEADER_H
#define M3DC1_SCOREC_HEADER_H

#define M3DC1_PI 3.141592653589793
#define FIXSIZEBUFF 1024
#define M3DC1_SUCCESS 0
#define M3DC1_FAILURE 1
#define C1TRIDOFNODE 6

#include "name_convert.h"
#include <mpi.h>
extern MPI_Comm M3DC1_COMM_WORLD;

#ifdef __cplusplus
extern "C" {
#endif
typedef int FieldID;
enum m3dc1_coord_system { /*0*/ M3DC1_RZPHI,  // default
                          /*1*/ M3DC1_XYZ};

// M3DC1_SOLVER_ORDER should be default -Fan
enum m3dc1_ordering { /*0*/ M3DC1_NO_ORDER=0,  // no reordering applied - default
                      /*2*/ M3DC1_ADJ_ORDER, // use adjaceny-based reordering on local mesh;
                      /*2*/ M3DC1_SOLVER_ORDER, // use solver's reordering;
                      /*3*/ M3DC1_ADJ_SOLVER_ORDER}; // use both adjaceny-based and solver's reordering

enum m3dc1_field_type { /*0*/ M3DC1_REAL=0,  // real number for field value
                        /*1*/ M3DC1_COMPLEX}; // complex number for field value

enum m3dc1_matrix_type { /*0*/ M3DC1_MULTIPLY=0,
                         /*1*/ M3DC1_SOLVE};

enum m3dc1_matrix_status { /*0*/ M3DC1_NOT_FIXED=0,
                           /*2*/ M3DC1_FIXED};

bool m3dc1_double_isequal(double A, double B);

void m3dc1_scorec_setcomm_c(MPI_Comm cm);
void m3dc1_scorec_setcomm_f(MPI_Fint * cm);
void m3dc1_scorec_init(int * argc, char ** argv[]);
void m3dc1_scorec_finalize();

/** plane functions */
void m3dc1_plane_setnum(int*);
void m3dc1_plane_getnum(int*);
void m3dc1_plane_getid(int * /* out */ plane_id);
void m3dc1_plane_setphirange(double* minValue, double* maxValue);
void m3dc1_plane_setphi(int* planeid, double* phi);
void m3dc1_plane_getphi(int* planeid, double* phi);

/** model functions */
void m3dc1_model_load(char* /* in */ model_file);
void m3dc1_model_print();
void m3dc1_model_setnumplane(int*);
void m3dc1_model_getnumplane(int*);
void m3dc1_model_getplaneid(int * /* out */ plane_id);

void m3dc1_model_getmincoord(double* /* out */ x_min, double* /* out */ y_min); //getmincoord2_
void m3dc1_model_getmaxcoord(double* /* out */ x_max, double* /* out */ y_max); //getmaxcoord2_

/** mesh functions */

void m3dc1_mesh_load(char* mesh_file);
void m3dc1_mesh_write(char* filename, int *option); // 0: vtk file with field; 1:smb file
void m3dc1_mesh_build3d(int* num_field, int* fid, int* num_dofs_per_value);

void m3dc1_ghost_create (int* num_layer );
void m3dc1_ghost_delete ();

void m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent);
void m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent); //numownedents_
void m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ global_num_ent); //numglobalents_
void m3dc1_mesh_getnumghostent (int* /* in*/ ent_dim, int* /* out */ num_ent);

void m3dc1_mesh_search(int* initial_simplex, double* final_position, int* final_simplex);

/* mesh entity functions */
void m3dc1_ent_getglobalid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ global_ent_id);
void m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id,
		            int* /* out */ geom_class_dim, int* /* out */ geom_class_id);
void m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* in */ adj_dim,
                      int* /* out */ adj_ent, int* /* in */ adj_ent_allocated_size, int* /* out */ num_adj_ent);
void m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* in */ adj_dim, int* /* out */ num_adj_ent);
void m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ owning_partid); //entprocowner_
void m3dc1_ent_isowner (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ ismine);
void m3dc1_ent_isghost(int* /* in */ ent_dim, int* /* in */ ent_id, int* isghost);

// node-specific functions
void m3dc1_node_getglobalid (int* ent_dim, int* /* in */ ent_id, int* /* out */ global_ent_id);
void m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord );
void m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyz);
void m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv);
void m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry);
void m3dc1_node_write (const char* filename, int* start_index);

// region-specific function
// only used in 3D
void m3dc1_region_getoriginalface( int * /* in */ elm, int * /* out */ fac);

/** field manangement */
void m3dc1_field_getnewid (FieldID* /*out*/ fid);
// ordering should be reused for field and matrix??? -Fan
// is num_dofs input or output?
// *value_type is either M3DC1_REAL or M3DC1_COMPLEX
void m3dc1_field_create (FieldID* /*in*/ fid, const char* /* in */ field_name, int* num_values, int* value_type, int* num_dofs_per_value);
void m3dc1_field_delete (FieldID* /*in*/ fid);

void m3dc1_field_getinfo(FieldID* /*in*/ fid, char* /* out*/ field_name, int* num_values, int* value_type, int* total_num_dof);

void m3dc1_field_exist(FieldID* fid, int * exist);//checkppplveccreated_
void m3dc1_field_sync (FieldID* /* in */ fid); // updatesharedppplvecvals_;
void m3dc1_field_sum (FieldID* /* in */ fid); // sumsharedppplvecvals_
void m3dc1_field_sumsq (FieldID* /* in */ fid, double* /* out */ sum);

/** field dof functions */
void m3dc1_field_getlocaldofid (FieldID* fid, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);
void m3dc1_field_getowndofid (FieldID* fid, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);
void m3dc1_field_getglobaldofid ( FieldID* fid, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);
void m3dc1_field_getghostdofid (FieldID* fid, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);

void m3dc1_field_getnumlocaldof (FieldID* fid, int* /* out */ num_local_dof);
void m3dc1_field_getnumowndof (FieldID* fid, int* /* out */ num_own_dof);
void m3dc1_field_getnumglobaldof (FieldID* fid, int* /* out */ num_global_dof);
void m3dc1_field_getnumghostdof (FieldID* fid, int* /* out */ num_ghost_dof);

void m3dc1_field_getdataptr (FieldID* fid, double** pts);

void m3dc1_field_add(FieldID* /*inout*/ fid1, FieldID* /*in*/ fid2);
void m3dc1_field_mult(FieldID* /*inout*/ fid, double* fac, int * scalar_type);
void m3dc1_field_assign(FieldID* /*inout*/ fid, double* fac, int * scalar_type);
void m3dc1_field_copy(FieldID* /* out */ fid1, FieldID* /* in */ fid2);

void m3dc1_field_retrieve (FieldID* /* in */ fid, double * /*out*/ data, int * /* in */size);
void m3dc1_field_set (FieldID* /* in */ fid, double * /*in*/ data, int * /* in */size);

void m3dc1_field_insert(FieldID* /* in */ fid, int* /* in */ s_dofid,
                        int* /* in */ size, double* /* in */ values,
                        int* scalar_type, int* /* in */ op);
void m3dc1_field_isnan(FieldID* /* in */ fid, int * isnan);
void m3dc1_field_compare(FieldID* fid1, FieldID* fid2);
// load fields from a file and return a field ID
void m3dc1_field_load (FieldID* /*out*/ fid, const char* /*in*/ filename);
void m3dc1_field_write(FieldID* fid, const char* filename, int* start_index);
void m3dc1_field_print(FieldID* fid);
void m3dc1_field_sum_plane (FieldID* /* in */ fid);
void m3dc1_field_printcompnorm(FieldID* /* in */ fid, const char * info);
void m3dc1_field_max (FieldID* fid, double * max_val, double * min_val);

void m3dc1_field_verify();

void m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* fid,
                       int* dof_id, int* /* out */ dof_cnt);
void m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* fid,
                       int* dof_id, int* /* out */ dof_cnt);
void m3dc1_ent_getnumdof (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* fid, int* /* out */ num_dof);
void m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* fid, int* /* ou
t */ num_dof, double* dof_data);
void m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* fid, int* /* out */ num_dof, double* dof_data);

#ifdef M3DC1_PETSC
/** matrix and solver functions with PETSc */
void m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* fid); //zerosuperlumatrix_
void m3dc1_matrix_assemble(int* matrix_id); //finalizematrix_
void m3dc1_matrix_delete(int* matrix_id); //deletematrix_
void m3dc1_matrix_reset(int* matrix_id);

void m3dc1_matrix_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
void m3dc1_matrix_add(int* matrix_id, int* row, int* column, int* scalar_type, double* val); //globalinsertval_

// insert all blocks related to the nodes effecting a given mesh entity of dimension ent_dim
//  for region elements this will insert all values into the matrix for all nodes effecting the element
//  this assembles an entire elemental matrix into the global matrix at once
void m3dc1_matrix_insertentblocks(int * mid, int ent_dim, int * eid, double * vals);
// insert all blocks related to a specific intersection of nodes from the specified element of
//  dimension ent_dim. for region elements this will insert the matrix blocks from rows associated
//  with node nd1 and columns from rows associated with node nd2.
//  this assembles specific sections of an elemental matrix into the global matrix
void m3dc1_matrix_insertnodeblocks(int * mid, int * ent_dim, int * eid, int * nd1, int * nd2, double * vals);
void m3dc1_matrix_setbc(int* matrix_id, int* row);
void m3dc1_matrix_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);

void m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol); //solveSysEqu_
void m3dc1_matrix_getnumiter(int* matrix_id, int * iter_num);
void m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, FieldID* outputvecid); //matrixvectormult_

// for performance test
void m3dc1_matrix_setassembleoption(int * op);
  void m3dc1_matrix_write(int * matrix_id, const char * file_name);
void m3dc1_matrix_print(int* matrix_id);
#endif // #ifdef M3DC1_PETSC

// adaptation
int adapt_by_field (int * fieldId, double* psi0, double * psil);
int set_adapt_p (double * pp);
int adapt_by_error_field (double * errorField, double * errorAimed, int* max_node, int* option); // option 0: local error control; 1 global

// for adaptation
int set_mesh_size_bound (double* abs_size, double * rel_size);
int set_adapt_smooth_factor (double* fac);
int output_face_data (int * size, double * data, const char * vtkfile);
int sum_edge_data (double * data, int * size);
int get_node_error_from_elm (double * elm_data, int * size, double* nod_data);

#ifdef M3DC1_TRILINOS
//=========================================================================
/** matrix and solver functions with TRILINOS */
//=========================================================================

void m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* fid);
void m3dc1_epetra_delete(int* matrix_id);
void m3dc1_epetra_reset(int* matrix_id);

void m3dc1_epetra_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
void m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values);

void m3dc1_epetra_setbc(int* matrix_id, int* row);
void m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);
void m3dc1_epetra_assemble(int* matrix_id);
void m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid);
void m3dc1_epetra_write(int* matrix_id, const char*, int* skip_zero, int* start_index);
void m3dc1_epetra_print(int* matrix_id);

void m3dc1_solver_aztec(int* matrix_id, FieldID* x_fieldid, FieldID*
		       b_fieldid, int* num_iter, double* tolerance,
		       const char* krylov_solver, const char*
		       preconditioner, const char* sub_dom_solver,
		       int* overlap, int* graph_fill, double*
		       ilu_drop_tol,  double* ilu_fill,
		       double* ilu_omega, int* poly_ord);

void m3dc1_solver_amesos(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid, const char* solver_name);
void m3dc1_solver_getnumiter(int* matrix_id, int * iter_num);
#endif //#ifdef M3DC1_TRILINOS

#ifdef __cplusplus
}
#endif
#endif
