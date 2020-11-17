/****************************************************************************** 

  (c) 2005-2020 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "name_convert.h"
#include "pumi.h"
#include <iostream>
#include <assert.h>
#include "m3dc1_mesh.h" // debugging purpose
#include <parma.h>
#include "PCU.h"
#include "petscksp.h"
#include <iostream>
#include <assert.h>

extern "C" int setPETScMat(int matrixid, Mat * A) {};
extern "C" int setPETScKSP(int matrixid, KSP * ksp, Mat * A){};

using namespace std;
static char help[] = "testing solver functions; \n first do mat-vec product A*b=c; solve Ax=c; compare x and b\n\n";

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff);
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif

#ifdef __bgq__
#include <spi/include/kernel/memory.h>

static double get_peak()
{
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return heap;
}

#elif defined (__linux__)

static double get_peak()
{
  return mallinfo().arena;
}

#else

static double get_peak()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}

#endif

static void print_stats(const char* name, double value)
{
  double min, max, avg;
  min = value;
  PCU_Min_Doubles(&min, 1);
  max = value;
  PCU_Max_Doubles(&max, 1);
  avg = value;
  PCU_Add_Doubles(&avg, 1);
  avg /= PCU_Comm_Peers();
  double imb = max / avg;
  if (!PCU_Comm_Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

#if defined(__linux__)

static double get_chunks()
{
  struct mallinfo m = mallinfo();
  return m.uordblks + m.hblkhd;
}

#else
static double get_chunks()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}
#endif

double t1, t2, t3, t4, t5, t6, t7;
int num_values=1, vertex_dim=0, scalar_type=0;
int num_dofs, num_dofs_node;
int b_field=1, c_field=2, x_field=3;

void test_matrix(int, int);

int main( int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();

  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  PetscLogDouble mem;
  if (argc<4 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes real(0)/complex(1) "<<endl;
    return M3DC1_FAILURE;
  }


  int num_plane=1;
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && PCU_Comm_Peers()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);
  }

  if (m3dc1_model_load(argv[1])) // model loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  m3dc1_model_print();

#ifdef PETSC_USE_COMPLEX
  scalar_type=1; // complex
#endif

  if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  //int three=3;
  //m3dc1_mesh_write("geoId", &three);

  printStats(m3dc1_mesh::instance()->mesh);
  int ent_dim=2;
  if (num_plane>1)
  {
    int zero=0;
    ent_dim=3;
    m3dc1_mesh_build3d(&zero, &zero, &zero);
  }
 
   int num_ent = 10;
  int* ent_ids=new int[num_ent];
  int* num_adj_ent=new int[num_ent];

  // M3DC1 assumption - local ID is continuous 
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(ent_dim-1);
  apf::MeshEntity* e;
  int i=0;
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    ent_ids[i++] =  get_ent_localid(m3dc1_mesh::instance()->mesh,e);
    if (i==10) break;
  }
  m3dc1_mesh::instance()->mesh->end(it);
  m3dc1_ent_getnumglobaladj (&ent_dim, ent_ids, &num_ent, &ent_dim, num_adj_ent);
/*
  for (int p=0; p<PCU_Comm_Peers();++p)
  {
    if (p==PCU_Comm_Self())   
      for (int i=0; i<10; ++i)
        std::cout<<"("<<p<<") # adj element of elem "<<i<<" = "<<num_adj_ent[i]<<"\n";
   MPI_Barrier(MPI_COMM_WORLD) ;
  }
    
*/

  int num_layer=2;
  m3dc1_ghost_create(&num_layer);
  pumi_mesh_verify(m3dc1_mesh::instance()->mesh, false);

   // set/get field dof values
  int num_vertex = m3dc1_mesh::instance()->mesh->count(0);
  int num_own_vertex = m3dc1_mesh::instance()->num_own_ent[0];

  int value_type[] = {scalar_type,scalar_type};

  num_dofs=6;
  if (num_plane>1) num_dofs=12;
  num_dofs_node = num_values * num_dofs;

  m3dc1_field_create (&b_field, "b_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&c_field, "c_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&x_field, "x_field", &num_values, value_type, &num_dofs);
  m3dc1_field_printcompnorm(&b_field, "b_field init info");
//  PetscMemoryGetCurrentUsage(&mem);
//  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
//  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<"* set b field ..."<<endl;
  // fill b field
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    // 2D mesh, z component =0
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &b_field, &num_dofs_node, &dofs.at(0));
  }
  m3dc1_field_printcompnorm(&b_field, "b_field after set info");
//  PetscMemoryGetCurrentUsage(&mem);
//  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
//  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<"* set matrix ..."<<endl; 
  t1 = MPI_Wtime();
  // fill matrix 
  // the matrix is diagnal dominant; thus should be positive definite
  int matrix_mult=1, matrix_solve=2;
  int matrix_mult_type = M3DC1_MULTIPLY;
  int matrix_solve_type = M3DC1_SOLVE;
  m3dc1_matrix_create(&matrix_mult, &matrix_mult_type, value_type, &b_field);
  m3dc1_matrix_create(&matrix_solve, &matrix_solve_type, value_type, &b_field);

  //check geo
  if(num_plane >1)
  {
    int num_elem = m3dc1_mesh::instance()->mesh->count(3);
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
        if (is_bdy)
        {
          int geom_class_dim,geom_class_id;
          m3dc1_ent_getgeomclass (&adj_dim, nodes+i, &geom_class_dim, &geom_class_id);
          m3dc1_node_getnormvec(nodes+i, normal1);
          m3dc1_node_getcurv(nodes+i,&curv1);
//          if (geom_class_dim==0)
//             cout<<"* node classified on geometric vertex with normal "<<normal1[0]<<" "<<normal1[1]<<" curv "<<curv1<<endl;

          m3dc1_node_isongeombdry(nodes+i+3, &is_bdy);
          assert(is_bdy);
          m3dc1_node_getnormvec(nodes+i+3, normal2);
          m3dc1_node_getcurv(nodes+i+3,&curv2);
          for(int i=0; i<2; i++)
            assert(AlmostEqualDoubles(normal1[i], normal2[i], 1e-6, 1e-6));
          assert(AlmostEqualDoubles(curv1,curv2, 1e-6, 1e-6));
        }
      }
    }
  }

  test_matrix(matrix_mult, matrix_solve);

//  m3dc1_matrix_print(&matrix_mult);
//  pumi_sync();

  if(!PCU_Comm_Self())
    cout<<"* time: fill matrix "<<t2-t1<<" assemble "<<t3-t2<<" mult "<<t4-t3
        <<" solve "<<t5-t4<<endl; 

  m3dc1_matrix_delete(&matrix_mult);
  m3dc1_matrix_delete(&matrix_solve);

  // test ghosting
  m3dc1_ghost_delete();
  pumi_mesh_verify(m3dc1_mesh::instance()->mesh, false);

  m3dc1_field_delete(&x_field);
  m3dc1_field_delete(&b_field);
  m3dc1_field_delete(&c_field);

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}

void test_matrix(int matrix_mult, int matrix_solve)
{
  int elem_dim= m3dc1_mesh::instance()->mesh->getDimension();

  double diag_value=2.0, off_diag=1.0;
  int node_elm = 3;
  if (elem_dim==3) node_elm=6;

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

  
  int num_elem = m3dc1_mesh::instance()->mesh->count(elem_dim);
  int num_vertex = m3dc1_mesh::instance()->mesh->count(0);
  int num_own_vertex = m3dc1_mesh::instance()->num_own_ent[0];

  for(int ielm = 0; ielm < num_elem; ielm++)
  {
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
         m3dc1_matrix_insertblock(&matrix_solve, &ielm, &rowVar, &colVar, &block_tmp[0]);
         m3dc1_matrix_insertblock(&matrix_mult, &ielm, &rowVar, &colVar, &block_tmp[0]);
      }
    }
  }
  t2 = MPI_Wtime();
  //PetscMemoryGetCurrentUsage(&mem);
  //PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  //PetscSynchronizedFlush(MPI_COMM_WORLD,NULL);
  if(!PCU_Comm_Self()) cout<<"* assemble matrix ..."<<endl;
  m3dc1_matrix_assemble(&matrix_mult);
  m3dc1_matrix_assemble(&matrix_solve); 
  //m3dc1_matrix_write(&matrix_mult, "matrixMult.m");
  //m3dc1_matrix_write(&matrix_solve, "matrixSolve.m");
  // print out memory usage from petsc
  //PetscMemoryGetCurrentUsage(&mem);
  //PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  //PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  t3 = MPI_Wtime();
  // calculate c field
  //m3dc1_field_print(&b_field);

  if(!PCU_Comm_Self()) cout<<"* do matrix-vector multiply ..."<<endl;
  m3dc1_matrix_multiply(&matrix_mult, &b_field, &c_field); 
  //m3dc1_field_print(&c_field);
  t4 = MPI_Wtime();
  // let's test field operations here
  m3dc1_field_copy(&x_field, &c_field);
  double val[]={2.};
  int realtype=0;
  m3dc1_field_mult(&x_field, val, &realtype);
  m3dc1_field_add(&x_field, &c_field);
  for(int i=0; i<num_vertex; i++)
  {
    vector<double> dofs1(num_dofs_node*(1+scalar_type)), dofs2(num_dofs_node*(1+scalar_type));
    int num_dofs_t;
    m3dc1_ent_getdofdata(&vertex_dim, &i, &c_field, &num_dofs_t, &dofs1.at(0));
    assert(num_dofs_t==num_dofs_node);
    m3dc1_ent_getdofdata(&vertex_dim, &i, &x_field, &num_dofs_t, &dofs2.at(0));
    assert(num_dofs_t==num_dofs_node);
    for(int i=0; i<num_dofs_node*((1+scalar_type)); i++)
      assert(AlmostEqualDoubles(dofs2.at(i), (val[0]+1)*dofs1.at(i), 1e-6, 1e-6));
  }
  // copy c field to x field
  m3dc1_field_copy(&x_field, &c_field);
  //m3dc1_field_print(&x_field);
  if(!PCU_Comm_Self()) cout<<"* solve ..."<<endl;
  // solve Ax=c
  int solver_type = 0;    // PETSc direct solver
  double solver_tol = 1e-6;
  m3dc1_matrix_solve(&matrix_solve, &x_field); //, &solver_type, &solver_tol);
  //m3dc1_field_print(&x_field);
  t5 = MPI_Wtime();
  // verify x=b
  for(int inode=0; inode<num_vertex; inode++)
  {
    vector<double> dofs_x(num_dofs_node*(1+scalar_type)), dofs_b(num_dofs_node*(1+scalar_type));
    m3dc1_ent_getdofdata(&vertex_dim, &inode, &b_field, &num_dofs_node, &dofs_b.at(0));
    m3dc1_ent_getdofdata(&vertex_dim, &inode, &x_field, &num_dofs_node, &dofs_x.at(0));
    for(int idof=0; idof<num_dofs_node*(1+scalar_type); idof++)
      assert(AlmostEqualDoubles(dofs_b.at(idof),dofs_x.at(idof), 1e-3, 1e-3));
  }
}

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff)
{
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

