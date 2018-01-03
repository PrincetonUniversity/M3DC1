/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef NOM3DC1
#define NOM3DC1
#endif
#include "m3dc1_scorec.h"
#include <iostream>
#include <assert.h>
#include "apf.h"
#include "petscksp.h"
#include "PCU.h"
#include "m3dc1_mesh.h" // debugging purpose
#include "apfMesh.h" // debugging purpose
using namespace std;
static char help[] = "testing solver functions; \n first do mat-vec product A*b=c; solve Ax=c; compare x and b\n\n";

int main( int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  int op=0;
  m3dc1_matrix_setassembleoption(&op);
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  PetscLogDouble mem;
  if (argc<4 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes solver_name real(0)/complex(1) "<<endl;
    return M3DC1_FAILURE;
  }
  double t1, t2, t3, t4, t5;
  double tf1,tf2;

#if 0
  int i, processid = getpid();
  if (!PCU_Comm_Self())
  {
    std::cout<<"Proc "<<PCU_Comm_Self()<<">> pid "<<processid<<" Enter any digit...\n";
    std::cin>>i;
  }
  else
    std::cout<<"Proc "<<PCU_Comm_Self()<<">> pid "<<processid<<" Waiting...\n";
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int num_plane=1;
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && PCU_Comm_Peers()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);
  }

  m3dc1_model_load(argv[1]);
  m3dc1_model_print();

  int scalar_type=0;
  if(argc>5) scalar_type=atoi(argv[5]);
#ifndef PETSC_USE_COMPLEX
  if(scalar_type!=0)
  {
    if(!PCU_Comm_Self()) cout<<"PETSc is not configured with --with-scalar-type=complex"<<endl;
    return 1;
  } 
;
#else
  if(!PCU_Comm_Self()) cout<<"PETSC_USE_COMPLEX defined"<<endl;
#endif

  m3dc1_mesh_load(argv[2]);
  //int three=3;
  //m3dc1_mesh_write("geoId", &three);

  printStats(m3dc1_mesh::instance()->mesh);
  if (num_plane>1)
  {
    int zero=0;
    m3dc1_mesh_build3d(&zero, &zero, &zero);
  }
  printStats(m3dc1_mesh::instance()->mesh);

  // set/get field dof values
  int num_vertex, num_own_vertex, vertex_dim=0;
  int num_elem, elem_dim=2;
  if(num_plane>1) elem_dim=3; // use wedge in 3D
  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);
  m3dc1_mesh_getnumownent (&vertex_dim, &num_own_vertex);
  m3dc1_mesh_getnument(&elem_dim, &num_elem);
#if 0
  for(int i=0; i<num_vertex; i++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&i, xyz);
    int ent_dim=0;
    int geom_class_dim;
    int geom_class_id;
    double normal[3];
    double curv;
    m3dc1_ent_getgeomclass ( &ent_dim, &i, &geom_class_dim, &geom_class_id);
    int bdy;
    m3dc1_node_isongeombdry(&i, &bdy);
    if(bdy) assert(geom_class_dim<2);
    if(geom_class_dim<2)
    {
      assert(bdy);
      m3dc1_node_getnormvec (&i, normal);
      m3dc1_node_getcurv (&i, &curv);
      cout<<" rank "<<PCU_Comm_Self()<<" node "<<i<<"("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<") geo class "<<geom_class_dim<<" "<<geom_class_id<<" normal curv "<<normal[0]<<" "<<normal[1]<<" "<<curv<<endl;
    }
    else
      cout<<" rank "<<PCU_Comm_Self()<<" node "<<i<<"("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<") geo class "<<geom_class_dim<<" "<<geom_class_id<<endl;
    // 2D mesh, z component =0
    if(num_plane==1) assert(m3dc1_double_isequal(xyz[2], 0));
  }
#endif

  int value_type[] = {scalar_type,scalar_type};

  int b_field=1, c_field=2, x_field=3;
#ifndef ASSEMBLEONLY
  int num_values=1;
#else
  int num_values=3;
#endif
  int num_dofs=6;
  if (num_plane>1) num_dofs=12;
  int num_dofs_node = num_values * num_dofs;
  m3dc1_field_create (&b_field, "b_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&c_field, "c_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&x_field, "x_field", &num_values, value_type, &num_dofs);
#ifdef M3DC1_TRILINOS
  int epetra_b_field=4, epetra_c_field=5, aztec_x_field=6;
  m3dc1_field_create (&epetra_b_field, "epetra_b_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&epetra_c_field, "epetra_c_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&aztec_x_field, "aztec_x_field", &num_values, value_type, &num_dofs);
#endif

  if (num_dofs%6==0) m3dc1_field_printcompnorm(&b_field, "b_field init info");
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<" set b field ..."<<endl;
  // fill b field
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    // 2D mesh, z component =0
    if(num_plane==1) assert(m3dc1_double_isequal(xyz[2], 0));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &b_field, &num_dofs_node, &dofs.at(0));
#ifdef M3DC1_TRILINOS
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &epetra_b_field, &num_dofs_node, &dofs.at(0));
#endif
  }

  if (num_dofs%6==0) m3dc1_field_printcompnorm(&b_field, "b_field after set info");
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<" set matrix ..."<<endl; 
  t1 = MPI_Wtime();
  // fill matrix 
  // the matrix is diagnal dominant; thus should be positive definite
  int matrix_mult=1, matrix_solve=2;
  int matrix_mult_type = M3DC1_MULTIPLY;
  int matrix_solve_type = M3DC1_SOLVE;
#ifndef ASSEMBLEONLY
  m3dc1_matrix_create(&matrix_mult, &matrix_mult_type, value_type, &b_field);
#endif
  m3dc1_matrix_create(&matrix_solve, &matrix_solve_type, value_type, &b_field);

#ifdef M3DC1_TRILINOS
  m3dc1_epetra_create(&matrix_mult, &matrix_mult_type, value_type, &epetra_b_field);
  if (!PCU_Comm_Self()) cout<<"> epetra mult matrix with b_field created ..."<<endl; 
  m3dc1_epetra_create(&matrix_solve, &matrix_solve_type, value_type, &epetra_b_field);
  if (!PCU_Comm_Self()) cout<<"> epetra solve matrix with b_field created ..."<<endl; 
#endif

  double diag_value=2.0, off_diag=1.0;
  int node_elm = 3;
  if(num_plane>1) node_elm=6;
  int num_dofs_element = num_dofs_node*node_elm;
  vector<double> block(num_dofs_element*num_dofs_element*(1+scalar_type),0);
  for(int i=0; i<num_dofs_element; i++)
  {
    for(int j=0; j<num_dofs_element; j++)
    {
      double val= (i==j? diag_value: off_diag);
      if(!scalar_type) block.at(i*num_dofs_element+j)=val;
      else
      {
        block.at(2*i*num_dofs_element+2*j)=val;
        block.at(2*i*num_dofs_element+2*j+1)=off_diag;
      }
    }
  }
  //check geo
#ifndef ASSEMBLEONLY
  if(num_plane >1)
  {
    for( int ielm = 0; ielm < num_elem; ielm++)
    {
      int ent_dim=3, adj_dim=0, num_adj_ent, adj_ent_allocated_size=6;
      int nodes[6];
      m3dc1_ent_getadj (&ent_dim, &ielm, &adj_dim, nodes, &adj_ent_allocated_size, &num_adj_ent);
      assert(num_adj_ent==adj_ent_allocated_size);
      double normal1[3], normal2[3], curv1, curv2;
      for(int i=0; i<3; i++)
      {
        int is_bdy=0;
        m3dc1_node_isongeombdry(nodes+i, &is_bdy);
        if(is_bdy)
        {
          int geom_class_dim,geom_class_id;
          m3dc1_ent_getgeomclass (&adj_dim, nodes+i, &geom_class_dim, &geom_class_id);
          m3dc1_node_getnormvec(nodes+i, normal1);
          m3dc1_node_getcurv(nodes+i,&curv1);
          if(geom_class_dim==0)
          {
             cout<<" Node classified on geometric vertex with normal "<<normal1[0]<<" "<<normal1[1]<<" curv "<<curv1<<endl;
          }
          m3dc1_node_isongeombdry(nodes+i+3, &is_bdy);
          assert(is_bdy);
          m3dc1_node_getnormvec(nodes+i+3, normal2);
          m3dc1_node_getcurv(nodes+i+3,&curv2);
          for(int i=0; i<2; i++)
            assert(m3dc1_double_isequal(normal1[i], normal2[i]));
          assert(m3dc1_double_isequal(curv1,curv2));
        }
      }
    }
  }
#endif
  double petsc_fill_time=0.0;
  double epetra_fill_time=0.0;
  for(int ielm = 0; ielm < num_elem; ielm++)
  {
    tf1 = MPI_Wtime();
    tf2 = MPI_Wtime();
    int nodes[256];
    int size_alloc=256;
    int num_nodes=-1;
    m3dc1_ent_getadj (&elem_dim, &ielm, &vertex_dim, nodes, &size_alloc, &num_nodes);
    for(int rowVar=0; rowVar< num_values; rowVar++)
    {
      for(int colVar=0; colVar< num_values; colVar++)
      {
         vector<double> block_tmp = block;
         if(rowVar!=colVar)
         {
           for(int i=0; i<block_tmp.size(); i++) block_tmp.at(i)*=0.5/num_values;
         }
         t1 = MPI_Wtime();
         m3dc1_matrix_insertblock(&matrix_solve, &ielm, &rowVar, &colVar, &block_tmp[0]);
         m3dc1_matrix_insertblock(&matrix_mult, &ielm, &rowVar, &colVar, &block_tmp[0]);
        petsc_fill_time+=MPI_Wtime()-t1;
        
#ifdef M3DC1_TRILINOS
         double epetra_t1 = MPI_Wtime();
         m3dc1_epetra_addblock(&matrix_solve, &ielm, &rowVar, &colVar, &block_tmp[0]);
         m3dc1_epetra_addblock(&matrix_mult, &ielm, &rowVar, &colVar, &block_tmp[0]);
         epetra_fill_time+=MPI_Wtime()-epetra_t1;
#endif
      }
    }
  }

  if(!PCU_Comm_Self()) cout<<" assemble/freeze matrix ..."<<endl;
  m3dc1_matrix_freeze(&matrix_mult);
  m3dc1_matrix_freeze(&matrix_solve); 
  assert(!m3dc1_epetra_freeze(&matrix_mult));
  assert(!m3dc1_epetra_freeze(&matrix_solve)); 

  // calculate c field
  if(!PCU_Comm_Self()) cout<<" do matrix-vector multiply Ab=c ..."<<endl;
  t3 = MPI_Wtime();
  m3dc1_matrix_multiply(&matrix_mult, &b_field, &c_field); 
  t4 = MPI_Wtime();
#ifdef M3DC1_TRILINOS
  double epetra_t3 = MPI_Wtime();
  m3dc1_epetra_multiply(&matrix_mult, &epetra_b_field, &epetra_c_field);  
  double epetra_t4 = MPI_Wtime();
  if (m3dc1_field_compare(&b_field, &epetra_b_field) == M3DC1_SUCCESS)
    if(!PCU_Comm_Self()) cout<<"> petsc_b==epetra_b verified\n";
  // verify petsc_c == epetra_c
  if (m3dc1_field_compare(&c_field, &epetra_c_field) == M3DC1_SUCCESS)
    if(!PCU_Comm_Self()) cout<<"> petsc_c==epetra_c verified\n";
#endif
  if(!PCU_Comm_Self()) cout<<" solve Ax=c with petsc ..."<<endl;
//  if (argc<4) 
//    ierr=m3dc1_solver_amesos(&matrix_solve, &x_field, &c_field, "superludist");
//  else
//    ierr=m3dc1_solver_amesos(&matrix_solve, &x_field, &c_field, argv[4]);
  m3dc1_field_copy(&x_field, &c_field);
  t5 = MPI_Wtime(); 
  m3dc1_matrix_solve(&matrix_solve, &x_field);
  if (!PCU_Comm_Self() && m3dc1_mesh::instance()->num_global_ent[0]<100) m3dc1_matrix_print(&matrix_solve);
  double t6 = MPI_Wtime();

  // verify petsc_x==petsc_b
  if (m3dc1_field_compare(&x_field, &b_field) == M3DC1_SUCCESS)
    if(!PCU_Comm_Self()) cout<<"> petsc_x==petsc_b verified\n";

  if(!PCU_Comm_Self()) cout<<" solve Ax=c with aztec..."<<endl;
  double aztec_t5 = MPI_Wtime();
  // Modified by kk
  // m3dc1_solver_aztec(&matrix_solve, &aztec_x_field, &epetra_c_field);
  int num_iter = 100;
  double tol = 1e-14;
  m3dc1_solver_aztec(&matrix_solve, &aztec_x_field, &epetra_c_field,
		     &num_iter, &tol, "gmres", "none");  
  double aztec_t6 = MPI_Wtime();

  // verify aztec_x==epetra_b
  if (m3dc1_field_compare(&aztec_x_field, &epetra_b_field) == M3DC1_SUCCESS)
    if(!PCU_Comm_Self()) cout<<"> aztec_x==epetra_b verified\n";

  // verify petsc_x==aztec_x
  if (m3dc1_field_compare(&x_field, &aztec_x_field) == M3DC1_SUCCESS)
    if(!PCU_Comm_Self()) cout<<"> petsc_x==aztec_x verified\n";

  if(!PCU_Comm_Self())
  {
    cout<<"> Petsc time: fill matrix "<<petsc_fill_time<<", mult "<<t4-t3<<", solve "<<t6-t5<<endl; 
    cout<<"> Epetra/Aztec time: fill matrix "<<epetra_fill_time<<", mult "<<epetra_t4-epetra_t3<<", solve "<<aztec_t6-aztec_t5<<endl; 
//    cout<<"> Amesos/SuperLU time: solve "<<t6-t5<<endl;
  }

  m3dc1_matrix_delete(&matrix_mult);
  m3dc1_matrix_delete(&matrix_solve);
  m3dc1_field_delete(&x_field);
  m3dc1_field_delete(&b_field);
  m3dc1_field_delete(&c_field);
#ifdef M3DC1_TRILINOS
  m3dc1_field_delete(&aztec_x_field);
  m3dc1_field_delete(&epetra_b_field);
  m3dc1_field_delete(&epetra_c_field);
  m3dc1_epetra_delete(&matrix_mult); 
  m3dc1_epetra_delete(&matrix_solve);
#endif
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}


