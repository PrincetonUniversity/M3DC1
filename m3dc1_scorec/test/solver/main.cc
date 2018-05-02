/******************************************************************************

  (c) 2005-2016 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h" // debugging purpose
#include "m3dc1_field.h"
#include <PCU.h>
#include <pumi.h>
#include <parma.h>
#include <cassert>
#include <fenv.h>
#include <iostream>
#include <stdlib.h>
using namespace std;
//static char help[] = "testing solver functions; \n first do mat-vec product A*b=c; solve Ax=c; compare x and b\n\n";
bool AlmostEqualDoubles(double A, double B, double maxDiff, double maxRelDiff);
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif

/*
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
*/

double t1, t2, t3, t4, t5, t6, t7;
int num_values = 1;
int vrt_dim = 0;
int scalar_type = 0;
int num_dofs = 0;
int dofs_per_nd = 0;
int b_field = 1;
int c_field = 2;
int x_field = 3;

void test_matrix(int, int);

int main(int argc, char * argv[])
{
  // abort on invalid floating point operations and values
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  m3dc1_scorec_init(&argc,&argv);
  if (argc<4 && !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes real(0)/complex(1) "<<endl;
    return M3DC1_FAILURE;
  }
  int num_plane = 1;
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane > 1 && PCU_Comm_Peers() % num_plane==0)
      m3dc1_model_setnumplane(&num_plane);
  }
  m3dc1_model_load(argv[1]);
  m3dc1_model_print();
#ifdef PETSC_USE_COMPLEX
  scalar_type = 1; // complex
#endif
  m3dc1_mesh_load(argv[2]);
  printStats(m3dc1_mesh::instance()->get_mesh());
  if (num_plane > 1)
  {
    int zero = 0;
    m3dc1_mesh_build3d(&zero, &zero, &zero);
  }
  //check geo
  if(num_plane > 1)
  {
    int num_elem = m3dc1_mesh::instance()->get_mesh()->count(3);
    for(int eid = 0; eid < num_elem; eid++)
    {
      int ent_dim = 3;
      int vrt_dim = 0;
      int num_adj_vrts = 0;
      m3dc1_ent_getnumadj(&ent_dim, &eid, &vrt_dim, &num_adj_vrts);
      int * nds = new int[num_adj_vrts];
      m3dc1_ent_getadj(&ent_dim, &eid, &vrt_dim, &nds[0], &num_adj_vrts, &num_adj_vrts);
      double normal1[] = {0.0,0.0,0.0};
      double normal2[] = {0.0,0.0,0.0};
      double curv1 = 0.0;
      double curv2 = 0.0;
      for(int ii = 0; ii < 3; ii++)
      {
        int is_bdy=0;
        m3dc1_node_isongeombdry(&nds[ii], &is_bdy);
        if (is_bdy)
        {
          int geom_class_dim = -1;
          int geom_class_id = -1;
          m3dc1_ent_getgeomclass(&vrt_dim, &nds[ii], &geom_class_dim, &geom_class_id);
          m3dc1_node_getnormvec(&nds[ii], normal1);
          m3dc1_node_getcurv(&nds[ii],&curv1);
          if(geom_class_dim==0)
            std::cout << "* node classified on geometric vertex with normal " << normal1[0] << " " << normal1[1] << " curv " << curv1 << std::endl;
          m3dc1_node_isongeombdry(&nds[ii+3], &is_bdy);
          assert(is_bdy);
          m3dc1_node_getnormvec(&nds[ii+3], normal2);
          m3dc1_node_getcurv(&nds[ii+3],&curv2);
          for(int jj = 0; jj < 2; jj++)
            assert(AlmostEqualDoubles(normal1[jj], normal2[jj], 1e-6, 1e-6));
          assert(AlmostEqualDoubles(curv1,curv2, 1e-6, 1e-6));
        }
      }
    }
  }
  // check field I/O
  if (argc>5)
  {
    int field_13;
    m3dc1_field_load(&field_13, argv[4]);
    m3dc1_field* mf = m3dc1_mesh::instance()->get_field(field_13);
    write_field(m3dc1_mesh::instance()->get_mesh(), mf, argv[5], 0);
    m3dc1_field_delete(&field_13);
  }

  //int num_layer=2;
  //m3dc1_ghost_create(&num_layer);
  //pumi_mesh_verify(m3dc1_mesh::instance()->get_mesh(), false);

   // set/get field dof values
  int num_vertex = m3dc1_mesh::instance()->get_mesh()->count(0);
  //int num_own_vertex = m3dc1_mesh::instance()->num_own_ent[0];

  int value_type[] = {scalar_type,scalar_type};

  num_dofs = 6;
  if (num_plane>1) num_dofs = 12;
  dofs_per_nd = num_values * num_dofs;

  m3dc1_field_create (&b_field, "b_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&c_field, "c_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&x_field, "x_field", &num_values, value_type, &num_dofs);
  m3dc1_field_printcompnorm(&b_field, "b_field init info");

  if(!PCU_Comm_Self()) cout<<"* set b field ..."<<endl;
  // fill b field
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    // 2D mesh, z component =0
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(dofs_per_nd*(1+scalar_type));
    for(int i=0; i<dofs_per_nd*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vrt_dim, &inode, &b_field, &dofs_per_nd, &dofs.at(0));
  }
  double fac[2] = {2.0, 2.0};
  m3dc1_field_mult(&b_field, fac, value_type);

  m3dc1_field_printcompnorm(&b_field, "b_field after set info");
  fac[0] = 0.5;   fac[1] = 0.5;
  m3dc1_field_mult(&b_field, fac, value_type);
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

  test_matrix(matrix_mult, matrix_solve);

//  m3dc1_matrix_print(&matrix_mult);
//  pumi_sync();

  t6 = MPI_Wtime();
  m3dc1_matrix_reset(&matrix_mult);
  m3dc1_matrix_reset(&matrix_solve);
  t7 = MPI_Wtime();

  if(!PCU_Comm_Self())
    cout<<"* time: fill matrix "<<t2-t1<<" assemble "<<t3-t2<<" mult "<<t4-t3
        <<" solve "<<t5-t4<<" reset "<<t7-t6<<endl;

  test_matrix(matrix_mult, matrix_solve);

  m3dc1_matrix_delete(&matrix_mult);
  m3dc1_matrix_delete(&matrix_solve);

  // test ghosting
  // m3dc1_ghost_delete();
  //pumi_mesh_verify(m3dc1_mesh::instance()->get_mesh(), false);

  m3dc1_field_delete(&x_field);
  m3dc1_field_delete(&b_field);
  m3dc1_field_delete(&c_field);

  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}

void test_matrix(int mat_mlt, int mat_slv)
{
  int elt_dim = m3dc1_mesh::instance()->get_mesh()->getDimension();
  assert(elt_dim == 3 || elt_dim == 2);
  int nds_per_elt = (elt_dim == 3 ? 6 : 3);
  int dofs_per_elt = dofs_per_nd * nds_per_elt;
  std::vector<double> elt_diag_blk(dofs_per_elt * dofs_per_elt * (1+scalar_type),0);
  for(int ii = 0; ii < dofs_per_elt; ii++)
  {
    for(int jj = 0; jj < dofs_per_elt; jj++)
    {
      double diag_val = (ii == jj ? 2.0 : 1.0); // 2 on main diag, 1 else
      if(!scalar_type)
        elt_diag_blk[ii * dofs_per_elt + jj] = diag_val;
      else
      {
        elt_diag_blk[2 * (ii * dofs_per_elt + jj)] = diag_val;
        elt_diag_blk[2 * (ii * dofs_per_elt + jj) + 1] = 1.0;
      }
    }
  }
  int num_lcl_elt = m3dc1_mesh::instance()->get_mesh()->count(elt_dim);
  int num_lcl_vrt = m3dc1_mesh::instance()->get_mesh()->count(0);
  if(!PCU_Comm_Self())
    std::cout << "* matrix elemental assembly ..." << std::endl;
  std::vector<double> & blk = elt_diag_blk; //(nd1 == nd2 ? elt_diag_blk : elt_off_blk);
  for(int elt_id = 0; elt_id < num_lcl_elt; ++elt_id)
    m3dc1_matrix_insertentblocks(&mat_mlt, &elt_dim, &elt_id, &blk[0]);//&nd1, &nd2, &blk[0]);
  for(int elt_id = 0; elt_id < num_lcl_elt; ++elt_id)
    m3dc1_matrix_insertentblocks(&mat_slv, &elt_dim, &elt_id, &blk[0]);//&nd1, &nd2, &blk[0]);
  t2 = MPI_Wtime();
  if(!PCU_Comm_Self())
  {
    std::cout << "* matrix parallel assembly ... " << std::endl
              << "** multiplication matrix" << std::endl;
  }
  m3dc1_matrix_assemble(&mat_mlt);
  m3dc1_matrix_assemble(&mat_slv);
  if(!PCU_Comm_Self())
    std::cout << "** solution matrix" << std::endl;
  t3 = MPI_Wtime();
  // calculate c field
  if(!PCU_Comm_Self())
    std::cout << "* do matrix-vector multiply ..." << std::endl;
  m3dc1_matrix_multiply(&mat_mlt, &b_field, &c_field);
  //m3dc1_field_print(&c_field);
  t4 = MPI_Wtime();
  // let's test field operations here
  m3dc1_field_copy(&x_field, &c_field);
  double val[]={2.};
  int realtype=0;
  m3dc1_field_mult(&x_field, val, &realtype);
  m3dc1_field_add(&x_field, &c_field);
  for(int ii = 0; ii < num_lcl_vrt; ii++)
  {
    std::vector<double> dofs1(dofs_per_nd*(1+scalar_type));
    std::vector<double> dofs2(dofs_per_nd*(1+scalar_type));
    int num_dofs_t;
    m3dc1_ent_getdofdata(&vrt_dim, &ii, &c_field, &num_dofs_t, &dofs1.at(0));
    assert(num_dofs_t==dofs_per_nd);
    m3dc1_ent_getdofdata(&vrt_dim, &ii, &x_field, &num_dofs_t, &dofs2.at(0));
    assert(num_dofs_t==dofs_per_nd);
    for(int jj = 0; jj < dofs_per_nd*((1+scalar_type)); jj++)
      assert(AlmostEqualDoubles(dofs2.at(jj), (val[0]+1)*dofs1.at(jj), 1e-6, 1e-6));
  }
  // copy c field to x field
  m3dc1_field_copy(&x_field, &c_field);
  //m3dc1_field_print(&x_field);
  if(!PCU_Comm_Self()) cout<<"* solve ..."<<endl;
  // solve Ax=c
  m3dc1_matrix_solve(&mat_slv, &x_field);
  t5 = MPI_Wtime();
  // verify Ax=b
  // set c = Ax
  m3dc1_matrix_multiply(&mat_slv, &x_field, &c_field);
  // compare c == b
  for(int vrt = 0; vrt < num_lcl_vrt; vrt++)
  {
    std::vector<double> dofs_c(dofs_per_nd*(1+scalar_type));
    std::vector<double> dofs_b(dofs_per_nd*(1+scalar_type));
    m3dc1_ent_getdofdata(&vrt_dim, &vrt, &c_field, &dofs_per_nd, &dofs_c[0]);
    m3dc1_ent_getdofdata(&vrt_dim, &vrt, &b_field, &dofs_per_nd, &dofs_b[0]);
    for(int idof=0; idof<dofs_per_nd*(1+scalar_type); idof++)
      assert(AlmostEqualDoubles(dofs_c.at(idof),dofs_b.at(idof), 1e-3, 1e-3));
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
