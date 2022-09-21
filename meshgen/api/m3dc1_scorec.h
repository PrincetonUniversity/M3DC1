/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
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
#define NODE_GLB_ORDER -1
#include "name_convert.h"

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

// Global variable to keep track of ghosted mesh
enum ghosting {uninit, init};
static enum ghosting global_ghost_state = uninit;

 int m3dc1_scorec_init();
int m3dc1_scorec_finalize();

//=========================================================================
/** plane functions */
//=========================================================================

int m3dc1_plane_setnum(int*);
int m3dc1_plane_getnum(int*);
int m3dc1_plane_getid(int * /* out */ plane_id);
int m3dc1_plane_setphirange(double* minValue, double* maxValue);
int m3dc1_plane_setphi(int* planeid, double* phi);
int m3dc1_plane_getphi(int* planeid, double* phi);

//=========================================================================
/** model functions */
//=========================================================================

int m3dc1_model_load(const char* /* in */ model_file);
int m3dc1_model_print();
int m3dc1_model_setnumplane(int*);
int m3dc1_model_getnumplane(int*);
int m3dc1_model_getplaneid(int * /* out */ plane_id);

int m3dc1_model_setpbc (int* /* in */ x_pbc, int* /* in */ y_pbc);
int m3dc1_model_getpbc (int* /* out */ x_pbc, int* /* out */ y_pbc); 
int m3dc1_model_getedge (int*  /* out */  left_edge, int*  /* out */  right_edge,
                         int*  /* out */  bottom_edge, int*  /* out */  top_edge); 
int m3dc1_model_getmincoord(double* /* out */ x_min, double* /* out */ y_min); 
int m3dc1_model_getmaxcoord(double* /* out */ x_max, double* /* out */ y_max); 

//=========================================================================
/** mesh functions */
//=========================================================================

int m3dc1_mesh_load(const char* mesh_file);
int m3dc1_mesh_build3d(int* num_field, int* field_id, int* num_dofs_per_value);

int m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent);
int m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent); 
int m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ global_num_ent); 

int m3dc1_mesh_setordering (int* option ); 
int m3dc1_mesh_getordering (int* option); 

// adjacency search
int m3dc1_mesh_search(int* initial_simplex, double* final_position, int* final_simplex);

void m3dc1_mesh_print();
int m3dc1_mesh_write(const char* filename, int *option); // 0: vtk file with field; 1:smb file

/** ghosting functions */
int m3dc1_ghost_load(int* nlayers); // For loading ghosted mesh with nlayers
int m3dc1_ghost_delete(); // For deleting a ghosted mesh

//=========================================================================
/* mesh entity functions */
//=========================================================================

int m3dc1_ent_getglobalid(int* ent_id, int* ent_dim, int* global_ent_id);
int m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id, 
		int* /* out */ geom_class_dim, int* /* out */ geom_class_id); 
int m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* in */ adj_dim,
		int* /* out */ adj_ent, int* /* in */ adj_ent_allocated_size, int* /* out */ num_adj_ent);
int m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* in */ adj_dim,
		int* /* out */ num_adj_ent);
int m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ owning_partid); 
int m3dc1_ent_isowned (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ ismine);
int m3dc1_ent_isghost (int* /* in */ ent_dim, int* /* in */ ent_id, int* /* out */ isghost);


int m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
		int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one); 
int m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
		int* /* out */ start_global_dof_id, int* /* out */ end_global_dof_id_plus_one);
int m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
		/* out */ int* num_dof, double* dof_data);
int m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
		int* /* out */ num_dof, double* dof_data);

// node functions
int m3dc1_node_getglobalid (int* /* in */ node_id, int* /* out */ global_node_id);   
int m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord ); 
int m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyz);
int m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv);
int m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry); 
int m3dc1_node_write(const char* filename, int* start_index);

// region functions
int m3dc1_region_getoriginalface( int * /* in */ elm, int * /* out */ fac);

//=========================================================================
/** field manangement */
//=========================================================================

int m3dc1_field_getnewid (FieldID* /*out*/field_id);

// ordering should be reused for field and matrix??? -Fan
// is num_dofs input or output?
// *value_type is either M3DC1_REAL or M3DC1_COMPLEX
// *num_dofs_per_value refers to the spatial derivatives
//                     stored for each "value", namely,
//                     a value "u" has these derivatives:
//                     u, u_x, u_y, u_xx, u_xy, u_yy
//                     not necessarily in that order.
// so, the hierarchy is:
//   values, for which we store
//   dofs which are either
//   real or complex scalars
int m3dc1_field_create (FieldID* /*in*/ field_id, const char* /* in */ field_name,
    		int* num_values, int* value_type, int* num_dofs_per_value);
int m3dc1_field_delete (FieldID* /*in*/ field_id); 

int m3dc1_field_getinfo(FieldID* /*in*/ field_id, char* /* out*/ field_name, int* num_values, int* value_type, int* total_num_dof);

int m3dc1_field_exist(FieldID* field_id, int * exist);
int m3dc1_field_sync (FieldID* /* in */ field_id); 
int m3dc1_field_sum (FieldID* /* in */ field_id); 
int m3dc1_field_sumsq (FieldID* /* in */ field_id, double* /* out */ sum);

/** field dof functions */
int m3dc1_field_getlocaldofid (FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one); 
int m3dc1_field_getowndofid (FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);  
int m3dc1_field_getglobaldofid ( FieldID* field_id, int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one);  

int m3dc1_field_getnumlocaldof (FieldID* field_id, int* /* out */ num_local_dof);
int m3dc1_field_getnumowndof (FieldID* field_id, int* /* out */ num_own_dof);
int m3dc1_field_getnumglobaldof (FieldID* field_id, int* /* out */ num_global_dof);
int m3dc1_field_getdataptr (FieldID* field_id, double** pts);

int m3dc1_field_add(FieldID* /*inout*/ field1, FieldID* /*in*/ field2);
int m3dc1_field_mult(FieldID* /*inout*/ field, double* fac, int * scalar_type);
int m3dc1_field_assign(FieldID* /*inout*/ field, double* fac, int * scalar_type);
int m3dc1_field_copy(FieldID* /* out */ filed1, FieldID* /* in */ field2);
int m3dc1_field_retrieve (FieldID* /* in */ filed1, double * /*out*/ data, int * /* in */size);
int m3dc1_field_set (FieldID* /* in */ filed1, double * /*in*/ data, int * /* in */size);
int m3dc1_field_insert(FieldID* /* in */ field, int /* in */ * local_dof, int * /* in */ size, 
		double* /* in */ values, int * type, int * /* in */ op);
int m3dc1_field_isnan(FieldID* /* in */ field, int * isnan);
int m3dc1_field_compare(FieldID* field_id_1, FieldID* field_id_2);
//int m3dc1_field_write(FieldID* field, const char* file_name);
int m3dc1_field_write(FieldID* field, const char* filename, int* start_index);
int m3dc1_field_print(FieldID* field);
int m3dc1_field_sum_plane (FieldID* /* in */ field_id);
int m3dc1_field_printcompnorm(FieldID* /* in */ field_id, const char* info);
int m3dc1_field_max (FieldID* field_id, double * max_val, double * min_val);

#ifdef M3DC1_PETSC
//=========================================================================
/** matrix and solver functions with PETSc */
//=========================================================================

int m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id); 
int m3dc1_matrix_freeze(int* matrix_id);
int m3dc1_matrix_delete(int* matrix_id);

int m3dc1_matrix_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
int m3dc1_matrix_add(int* matrix_id, int* row, int* column, int* scalar_type, double* val);

int m3dc1_matrix_insertblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values);
int m3dc1_matrix_setbc(int* matrix_id, int* row);
int m3dc1_matrix_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);

int m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol, int* solver_type, double* solver_tol); 
int m3dc1_matrix_getiternum(int* matrix_id, int * iter_num);
int m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, FieldID* outputvecid); 

int m3dc1_matrix_flush(int* matrix_id);
int m3dc1_matrix_setassembleoption(int * op);
int m3dc1_matrix_write(int* matrix_id, const char* file_name, int* start_index);
int m3dc1_matrix_print(int* matrix_id);
#endif // #ifdef M3DC1_PETSC

//=========================================================================
/** adaptation functions */
//=========================================================================

int adapt_by_field (int * fieldId, double* psi0, double * psil);
int set_adapt_p (double * pp);
int adapt_by_error_field (double * errorField, double * errorAimed, int* max_node, 
                 int* option); // option 0: local error control; 1 global
// helper functions for adaptation
int set_mesh_size_bound (double* abs_size, double * rel_size);
int set_adapt_smooth_factor (double* fac);
int sum_edge_data (double * data, int * size);
int get_node_error_from_elm (double * elm_data, int * size, double* nod_data);
int output_face_data (int * size, double * data, const char* vtkfile);

#ifdef M3DC1_TRILINOS
//=========================================================================
/** matrix and solver functions with TRILINOS */
//=========================================================================

int m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id);
int m3dc1_epetra_delete(int* matrix_id);

int m3dc1_epetra_insert(int* matrix_id, int* row, int* column, int* scalar_type, double* val);
int m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values);

int m3dc1_epetra_setbc(int* matrix_id, int* row);
int m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values);
int m3dc1_epetra_freeze(int* matrix_id); 
int m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid);
int m3dc1_epetra_write(int* matrix_id, const char*, int* skip_zero, int* start_index);
int m3dc1_epetra_print(int* matrix_id);

int m3dc1_solver_aztec(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid, int*, double*);
int m3dc1_solver_amesos(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid, const char* solver_name);
int m3dc1_solver_getnumiter(int* matrix_id, int * iter_num);
#endif //#ifdef M3DC1_TRILINOS
#ifdef __cplusplus
}
#endif
#endif


