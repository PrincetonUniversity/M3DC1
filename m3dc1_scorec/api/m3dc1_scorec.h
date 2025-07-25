/****************************************************************************** 


  (c) 2005-2023 Scientific Computation Research Center, 
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
#include "mpi.h"
#include <apf.h>
#include "m3dc1_sizeField.h"
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
                           /*1*/ M3DC1_FIXED,
                           /*2*/ M3DC1_SOLVED};

bool m3dc1_double_isequal(double A, double B);

int m3dc1_scorec_init();
int m3dc1_scorec_finalize();
void m3dc1_scorec_verbosity(int*);

/** plane functions */
int m3dc1_plane_setnum(int*);
int m3dc1_plane_getnum(int*);
int m3dc1_plane_getid(int * /* out */ plane_id);
int m3dc1_plane_setphirange(double* minValue, double* maxValue);
int m3dc1_plane_setphi(int* planeid, double* phi);
int m3dc1_plane_getphi(int* planeid, double* phi);

/** model functions */
int m3dc1_model_load(char* /* in */ model_file);
int m3dc1_modelinfo_load(char* /* in */ modelInfo_file);
void m3dc1_model_getgeometricloop(int* modelEdges, int* numberOfEdges, int* loopType);
void m3dc1_model_getmodeltype(int* modeltype);
int m3dc1_model_print();
int m3dc1_model_setnumplane(int*);
int m3dc1_model_getnumplane(int*);
void m3dc1_model_settopo();
int m3dc1_model_getmincoord(double* /* out */ x_min, double* /* out */ y_min); //getmincoord2_
int m3dc1_model_getmaxcoord(double* /* out */ x_max, double* /* out */ y_max); //getmaxcoord2_

/** mesh functions */

int m3dc1_mesh_load(char* mesh_file);
void m3dc1_mesh_load_3d(char* mesh_file, int* num_plane);

// option 0: vtk file with field; 1:smb file
int m3dc1_mesh_write(char* filename, int *option, int* /*time step*/);
void m3dc1_mesh_verify();

int m3dc1_mesh_build3d(int* num_field, int* field_id, int* num_dofs_per_value);

void m3dc1_mesh_adapt(int* field_id_h1, int* field_id_h2, double* dir);

int m3dc1_ghost_create (int* num_layer ); 
int m3dc1_ghost_delete ();

void m3dc1_mesh_getentid (int* /* in */ ent_dim, int* /* out */ ent_ids, 
            int* /* in */ allocated_size, int* /* out */ num_ent);
int m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent);
int m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent); //numownedents_
int m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ global_num_ent); //numglobalents_
int m3dc1_mesh_getnumghostent (int* /* in*/ ent_dim, int* /* out */ num_ent);

int m3dc1_mesh_search(int* initial_simplex, double* final_position, int* final_simplex);

// mesh-level operators with communications 
// at the moment, this works only for 2 kinds 2nd order adjacency
// which are face-edge-face or region-face-region (ent_dim and adj_dim are 2 for 2D, 3 for 3D)
void m3dc1_ent_getglobaladj (int* /* in */ ent_dim, 
                      int* /* in */ ent_ids, int* /* in */ num_ent,
                      int* /* in */ adj_dim,
                      int* /* out */ num_adj_ent, int* /* out */ adj_ent_gids, int* /* out */ adj_ent_pids,
                      int* /* in */ adj_ent_allocated_size, int* /* out */ adj_ent_size);
void m3dc1_ent_getnumglobaladj (int* /* in */ ent_dim, 
                      int* /* in */ ent_ids, int* /* in */ num_ent,
                      int* /* in */ adj_dim, int* /* out */ num_adj_ent);

void m3dc1_ent_getlocalid (int* /* in */ ent_dim, int* /* out */ ent_ids,
            int* /* in */ allocated_size, int* /* out */ num_ent);

// individual entity level operators
int m3dc1_ent_getglobalid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ global_ent_id);
int m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id,
                            int* /* out */ geom_class_dim, int* /* out */ geom_class_id);
int m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ entids, int* /* in */ adj_dim,
                      int* /* out */ adj_ent, int* /* in */ adj_ent_allocated_size, int* /* out */ adj_ent_size);
int m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* in */ adj_dim, int* /* out */ num_adj_ent);

int m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ owning_partid); //entprocowner_
int m3dc1_ent_isowner (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ ismine); 
int m3dc1_ent_isghost(int* /* in */ ent_dim, int* /* in */ ent_id, int* isghost);
void m3dc1_ent_measure(int* /* in */ ent_dim, int* /* in */ ent_id, double* /* out */ value);

// node-specific functions
void m3dc1_node_setfield (int* /* in */ node_id, int* /* in */ field_id, double* /* in */ values, int* /* in */ num_values);
void m3dc1_node_getfield (int* /* in */ node_id, int* /* in */ field_id, double* /* out */ values, int* /* out */ num_values);
int m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord ); 
int m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyz);
int m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv);
int m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry); 
int m3dc1_node_write (const char* filename, int* start_index);
void m3dc1_node_getNormVecOnNewVert(apf::MeshEntity* v, double* normalVec);

// region-specific function
int m3dc1_region_getoriginalface( int * /* in */ elm_id, int * /* out */ fac_id);

/** field manangement */
int m3dc1_field_getnewid (FieldID* /*out*/field_id);
// ordering should be reused for field and matrix??? -Fan
// is num_dofs input or output?
// *value_type is either M3DC1_REAL or M3DC1_COMPLEX
int m3dc1_field_create (FieldID* /*in*/ field_id, const char* /* in */ field_name, int* num_values, int* value_type, int* num_dofs_per_value);
int m3dc1_field_delete (FieldID* /*in*/ field_id); 

int m3dc1_field_getinfo(FieldID* /*in*/ field_id, char* /* out*/ field_name, int* num_values, int* value_type, int* total_num_dof);

int m3dc1_field_exist(FieldID* field_id, int * exist);//checkppplveccreated_
int m3dc1_field_sync (FieldID* /* in */ field_id); // updatesharedppplvecvals_;
int m3dc1_field_sum (FieldID* /* in */ field_id); // sumsharedppplvecvals_
int m3dc1_field_sumsq (FieldID* /* in */ field_id, double* /* out */ sum);

/** field dof functions */
int m3dc1_field_getlocaldofid (FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one); 
int m3dc1_field_getowndofid (FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);  
int m3dc1_field_getglobaldofid ( FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);  
int m3dc1_field_getghostdofid (FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);

int m3dc1_field_getnumlocaldof (FieldID* field_id, int* /* out */ num_local_dof);
int m3dc1_field_getnumowndof (FieldID* field_id, int* /* out */ num_own_dof);
int m3dc1_field_getnumglobaldof (FieldID* field_id, int* /* out */ num_global_dof);
int m3dc1_field_getnumghostdof (FieldID* field_id, int* /* out */ num_ghost_dof);

int m3dc1_field_getdataptr (FieldID* field_id, double** pts);

int m3dc1_field_add(FieldID* /*inout*/ field1, FieldID* /*in*/ field2);
int m3dc1_field_mult(FieldID* /*inout*/ field, double* fac, int * scalar_type);
int m3dc1_field_assign(FieldID* /*inout*/ field, double* fac, int * scalar_type);
int m3dc1_field_copy(FieldID* /* out */ filed1, FieldID* /* in */ field2);
int m3dc1_field_retrieve (FieldID* /* in */ filed1, double * /*out*/ data, int * /* in */size);
int m3dc1_field_set (FieldID* /* in */ filed1, double * /*in*/ data, int * /* in */size);
int m3dc1_field_insert(FieldID* /* in */ field, int /* in */ * local_dof, int * /* in */ size, double* /* in */ values, int * type, int * /* in */ op);
int m3dc1_field_isnan(FieldID* /* in */ field, int * isnan);
int m3dc1_field_compare(FieldID* field_id_1, FieldID* field_id_2);
int m3dc1_field_write(FieldID* field, const char* filename, int* start_index);
void m3dc1_dir_export(double*, int, int);
void m3dc1_dir_import(double*, int);
void m3dc1_field_import(int, int);
void m3dc1_field_export(int*, int*);
void m3dc1_field_importall();
void m3dc1_field_exportall();
int m3dc1_field_print(FieldID* field);
int m3dc1_field_sum_plane(FieldID* /* in */ field_id);
int m3dc1_field_printcompnorm(FieldID* /* in */ field_id, char* info);
int m3dc1_field_max(FieldID* field_id, double * max_val, double * min_val);

void m3dc1_field_verify();
void m3dc1_field_mark4tx(FieldID* /*in*/ field_id); //mark for solution transfer

int m3dc1_model_getplaneid(int * /* out */ plane_id);


int m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                       int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);  //entdofs_
int m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                       int* /* out */ start_global_dof_id, int* /* out */ end_global_dof_id_plus_one); //globalentdofs_
int m3dc1_ent_getnumdof (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, int* /* out */ num_dof);
int m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, int* /* ou
t */ num_dof, double* dof_data);
int m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, int* /* out */ num_dof, double* dof_data);

#ifdef M3DC1_PETSC
/** matrix and solver functions with PETSc */
int m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id); //zerosuperlumatrix_
int m3dc1_matrix_assemble(int* matrix_id); //finalizematrix_
int m3dc1_matrix_delete(int* matrix_id); //deletematrix_
void m3dc1_matrix_reset(int* matrix_id); // cleanMatrixValues_
void m3dc1_matrix_update(int* matrix_id); // updateAPMatrixValues_

int m3dc1_matrix_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
int m3dc1_matrix_add(int* matrix_id, int* row, int* column, int* scalar_type, double* val); //globalinsertval_

int m3dc1_matrix_insertblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values);
int m3dc1_matrix_setbc(int* matrix_id, int* row);
int m3dc1_matrix_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);

int m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol); //solveSysEqu_

int m3dc1_matrix_getnumiter(int* matrix_id, int * iter_num);
void m3dc1_matrix_solve_with_guess(int* matrix_id, FieldID* rhs_sol, FieldID* xVec_guess);
int m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, FieldID* outputvecid); //matrixvectormult_

// for performance test
int m3dc1_matrix_setassembleoption(int * op);
int m3dc1_matrix_write(int* matrix_id, const char* file_name, int* start_index);
int m3dc1_matrix_print(int* matrix_id);
#endif // #ifdef M3DC1_PETSC

// adaptation
int adapt_by_field (int * fieldId, double* psi0, double * psil);
int set_adapt_p (double * pp);
int adapt_by_error_field (double * errorField, double * errorAimed, int* max_node, 
		int* option); // option 0: local error control; 1 global
// adapt the mesh on specific model faces based on psi field
int adapt_model_face(int * fieldId, double* psi0, double * psil, int* iadaptFaceNumber);

// 3D Adaptation
void m3dc1_spr_adapt (int * fieldId, int * index, int * ts,
    double * ar, double * max_size, int * refine_level, int * coarsen_level, 
    bool* update);
int node_error_3d_mesh (double* elm_data, int* size, double* nod_data);
int find_sizefield(double* node_error, double * errorAimed, int * max_adapt_node, int * option);
// for adaptation
int set_mesh_size_bound (double* abs_size, double * rel_size);
int set_adapt_smooth_factor (double* fac);
void output_face_data (int * size, double * data, char * vtkfile);
int sum_edge_data (double * data, int * size);
int get_node_error_from_elm (double * elm_data, int * size, double* nod_data);

#ifdef M3DC1_TRILINOS
//=========================================================================
/** matrix and solver functions with TRILINOS */
//=========================================================================

int m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id);
int m3dc1_epetra_delete(int* matrix_id);
void m3dc1_epetra_reset(int* matrix_id);

int m3dc1_epetra_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
int m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values);

int m3dc1_epetra_setbc(int* matrix_id, int* row);
int m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);
int m3dc1_epetra_assemble(int* matrix_id); 
int m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid);
int m3dc1_epetra_write(int* matrix_id, const char*, int* skip_zero, int* start_index);
int m3dc1_epetra_print(int* matrix_id);

int m3dc1_solver_aztec(int* matrix_id, FieldID* x_fieldid, FieldID*
		       b_fieldid, int* num_iter, double* tolerance,
		       const char* krylov_solver, const char*
		       preconditioner, const char* sub_dom_solver,
		       int* overlap, int* graph_fill, double*
		       ilu_drop_tol,  double* ilu_fill,
		       double* ilu_omega, int* poly_ord);
  
int m3dc1_solver_amesos(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid, const char* solver_name);
int m3dc1_solver_getnumiter(int* matrix_id, int * iter_num);
#endif //#ifdef M3DC1_TRILINOS

// for internal debugging purpose
// checkMatrixStatus_
int m3dc1_matrix_getstatus (int* matrix_id, int* status);
// getMatrixLocalDofNum_
void m3dc1_matrix_getlocalnumdof(int* matrix_id, int *num_own_dof);
// getMatrixGlobalDofs_
void m3dc1_matrix_getglobalnumdof(int *matrix_id, int *num_global_dof);
// getMatrixPetscDnnzOnnz_
void m3dc1_matrix_getpetscdnnzonnz(int *matrix_id, int *valType, int *d_nnz, int *o_nnz);
// getMatrixNNZRowSize_
void gm3dc1_matrix_getnnzrowsize(int *matrix_id, int *valType, int *rowSize);
// getMatrixNNZRowId_
void m3dc1_matrix_getnnzrowid(int *matrix_id, int *valType, int* ith, int *rowId);
// getMatrixNNZColSize_
void  m3dc1_matrix_getnnzcolsize(int *matrix_id, int *valType, int *rowId, int *colSize);
// getMatrixNNZValues_
void m3dc1_matrix_getnnzvalue(int *matrix_id, int *valType, int *rowId, int *colId, double *dvalues);
// getMatrixFirstDof_
void m3dc1_matrix_getmatrixfirstdof(int *matrix_id, int *firstdof);
// setMatrixSoln_
void m3dc1_matrix_setsoln(int *matrix_id, int *valType, double *soln);
#ifdef __cplusplus
}
#endif
#ifdef NOM3DC1
#include "petscksp.h"
extern "C" void resizevec_(double *, int*) {}
extern "C" void space_() {}
extern "C" int setPETScMat(int matrixid, Mat * A)
{
  PetscErrorCode ierr;
  ierr = MatSetType(*A, MATMPIAIJ);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScMat %d to MATMPIAIJ\n", matrixid);
  ierr = MatSetFromOptions(*A);CHKERRQ(ierr);
  return 0;
}
extern "C" int setPETScKSP(int matrixid, KSP * ksp, Mat * A)
{
  PetscErrorCode ierr;
  ierr = KSPCreate(MPI_COMM_WORLD, ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(*ksp, *A, *A, SAME_PRECONDITIONER /*DIFFERENT_NONZERO_PATTERN*/);CHKERRQ(ierr);
  ierr = KSPSetTolerances(*ksp, .000001, .000000001,
                          PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScKSP for %d\n", matrixid);
  return 0;
}
#endif

#endif
