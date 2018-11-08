/******************************************************************************
  (c) 2005-2016 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h" // debugging purpose
#include "m3dc1_field.h"
#include "m3dc1_matrix.h"
#include <PCU.h>
#include <pumi.h>
#include <parma.h>
#include <cassert>
#include <fenv.h>
#include <iostream>
#include <stdlib.h>
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif
int b_fld = 1;
int c_fld = 2;
int x_fld = 3;
template <typename scalar>
bool close(scalar a, scalar b, scalar eps = 1e-8, scalar rel_eps = 1e-8)
{
  // http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
  // Check if the numbers are really close -- needed
  // when comparing numbers near zero.
  scalar diff = fabs(a - b);
  if (diff <= eps)
    return true;
  a = fabs(a);
  b = fabs(b);
  scalar largest = (b > a) ? b : a;
  if (diff <= largest * rel_eps)
    return true;
  return false;
}
void init_field(int fid)
{
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  int num_blks = 0;
  int dofs_per_blk = 0;
  int agg_scp = -1;
  m3dc1_field_getinfo(&fid,&num_blks,&dofs_per_blk,&agg_scp);
  int dofs_per_nd = dofs_per_blk * num_blks;
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  for(int vtx_idx = 0; vtx_idx < num_lcl_vtx; ++vtx_idx)
  {
    double xyz[3];
    m3dc1_node_getcoord(&vtx_idx, xyz);
    std::vector<double> dofs(dofs_per_nd * sz);
    for(int ii = 0; ii < dofs_per_nd * sz; ii++)
      dofs.at(ii) = xyz[ii%3];
    m3dc1_ent_setdofdata(&vtx_dim, &vtx_idx, &fid, &dofs_per_nd, &dofs.at(0));
  }
}
void test_field_ops(int fid1, int fid2)
{
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  int num_blks = 0;
  int dofs_per_blk = 0;
  int agg_scp = -1;
  m3dc1_field_getinfo(&fid1,&num_blks,&dofs_per_blk,&agg_scp);
  int dofs_per_nd = dofs_per_blk * num_blks;
  double scale = 2.0;
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  // let's test field operations here
  // f2 = f1
  m3dc1_field_copy(&fid2, &fid1);
  // f2 *= 2
  m3dc1_field_mult(&fid2, &scale);
  // f2 += f1
  m3dc1_field_add(&fid2, &fid1);
  // verify f2 == 3*f1
  for(int ii = 0; ii < num_lcl_vtx; ii++)
  {
    std::vector<double> dofs1(dofs_per_nd * sz);
    std::vector<double> dofs2(dofs_per_nd * sz);
    int num_dofs_t;
    m3dc1_ent_getdofdata(&vtx_dim, &ii, &fid1, &num_dofs_t, &dofs1.at(0));
    assert(num_dofs_t==dofs_per_nd);
    m3dc1_ent_getdofdata(&vtx_dim, &ii, &fid2, &num_dofs_t, &dofs2.at(0));
    assert(num_dofs_t==dofs_per_nd);
    for(int jj = 0; jj < dofs_per_nd * sz; jj++)
      assert(close<double>(dofs2.at(jj), (scale + 1) * dofs1.at(jj), 1e-6, 1e-6));
  }
}
void test_matrix(int mat_id)
{
  if(!PCU_Comm_Self())
    std::cout << "* init matrix " << mat_id << " ..." << std::endl;
  double t[7] = {0.0};
  int elt_dim = m3dc1_mesh::instance()->get_mesh()->getDimension();
  int vtx_dim = 0;
  int lhs_fld_id = -1;
  m3dc1_matrix_getfieldid(&mat_id,&lhs_fld_id);
  int blks_per_nd = 0;
  int dofs_per_blk = 0;
  int agg_scp = -1;
  m3dc1_field_getinfo(&lhs_fld_id,&blks_per_nd,&dofs_per_blk,&agg_scp);
  assert(elt_dim == 3 || elt_dim == 2);
  int nds_per_elt = (elt_dim == 3 ? 6 : 3);
  int dofs_per_nd = dofs_per_blk * blks_per_nd;
  int dofs_per_elt = dofs_per_nd * nds_per_elt;
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  int dbls_per_elt = dofs_per_elt * sz;
  std::vector<double> elt_diag_blk(dbls_per_elt * dbls_per_elt,0);
  for(int ii = 0; ii < dofs_per_elt; ii++)
  {
    for(int jj = 0; jj < dofs_per_elt; jj++)
    {
      double diag_val = (ii == jj ? 2.0 : 1.0); // 2 on main diag, 1 else
#ifndef PETSC_USE_COMPLEX
      elt_diag_blk[ii * dofs_per_elt + jj] = diag_val;
#else
      elt_diag_blk[2 * (ii * dofs_per_elt + jj)] = diag_val;
      elt_diag_blk[2 * (ii * dofs_per_elt + jj) + 1] = 1.0;
#endif
    }
  }
  int num_lcl_elt = m3dc1_mesh::instance()->get_mesh()->count(elt_dim);
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  if(!PCU_Comm_Self())
    std::cout << "* matrix elemental assembly ..." << std::endl;
  t[0] = MPI_Wtime();
  std::vector<double> & blk = elt_diag_blk;
  for(int elt_id = 0; elt_id < num_lcl_elt; ++elt_id)
    m3dc1_matrix_insertentblocks(&mat_id, &elt_dim, &elt_id, &blk[0]);
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(mat_id);
  mat->printAssemblyInfo();
  if(!PCU_Comm_Self())
    std::cout << "* matrix parallel assembly ... " << std::endl;
  t[1] = MPI_Wtime();
  m3dc1_matrix_assemble(&mat_id);
  t[2] = MPI_Wtime();
  if(!PCU_Comm_Self())
    std::cout << "  assembled in " << t[2] - t[1] << std::endl;
  mat->printNNZStats();
  // calculate c field
  if(!PCU_Comm_Self())
    std::cout << "* do matrix-vector multiply ..." << std::endl;
  // c = Ab
  m3dc1_matrix_multiply(&mat_id, &b_fld, &c_fld);
  t[3] = MPI_Wtime();
  // x = c
  m3dc1_field_copy(&x_fld, &c_fld);
  if(!PCU_Comm_Self())
    std::cout << "* solve ..." << std::endl;
  // solve Ax = b
  // x = (A^-1)b
  m3dc1_matrix_solve(&mat_id, &x_fld);
  t[4] = MPI_Wtime();
  // verify b == Ax
  // set c = Ax
  m3dc1_matrix_multiply(&mat_id, &x_fld, &c_fld);
  // compare c = Ax and b
  for(int vtx = 0; vtx < num_lcl_vtx; vtx++)
  {
    std::vector<double> dofs_c(dofs_per_nd);
    std::vector<double> dofs_b(dofs_per_nd);
    m3dc1_ent_getdofdata(&vtx_dim, &vtx, &c_fld, &dofs_per_nd, &dofs_c[0]);
    m3dc1_ent_getdofdata(&vtx_dim, &vtx, &b_fld, &dofs_per_nd, &dofs_b[0]);
    for(int ii=0; ii < dofs_per_nd; ii++)
      if(!close<double>(dofs_c.at(ii),dofs_b.at(ii), 1e-3, 1e-3))
      {
        std::cout << "b field and c field not 'close' "
                  << dofs_c.at(ii) << " c-field value and "
                  << dofs_b.at(ii) << " b-field value" << std::endl;
        assert(false);
      }
  }
  t[5] = MPI_Wtime();
  m3dc1_matrix_reset(&mat_id);
  t[6] = MPI_Wtime();
  if(!PCU_Comm_Self())
    std::cout << "* timings for matrix " << mat_id      << ":\n"
              << "\tfill matrix        " << t[1] - t[0] << "\n"
              << "\tassemble           " << t[2] - t[1] << "\n"
              << "\tmult               " << t[3] - t[2] << "\n"
              << "\tsolve              " << t[4] - t[3] << "\n"
              << "\treset              " << t[6] - t[5] << std::endl;
}
int main(int argc, char * argv[])
{
  // abort on invalid floating point operations and values
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  m3dc1_scorec_init(&argc,&argv);
  if (argc != 4 && !PCU_Comm_Self())
  {
    std::cout << "Usage: ./" << argv[0] << " [model_filename] [mesh_filename] [number_of_planes] " << std::endl;
    return M3DC1_FAILURE;
  }
  int num_plane = atoi(argv[3]);
  //int num_dof_blks = atoi(argv[5]);
  int scalar_tp = M3DC1_REAL;
#ifdef PETSC_USE_COMPLEX
  scalar_tp = M3DC1_COMPLEX;
#endif
  if (num_plane > 1 && PCU_Comm_Peers() % num_plane == 0)
    m3dc1_model_setnumplane(&num_plane);
  m3dc1_model_load(argv[1]);
  m3dc1_model_print();
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
      int vtx_dim = 0;
      int num_adj_vtxs = 0;
      m3dc1_ent_getnumadj(&ent_dim, &eid, &vtx_dim, &num_adj_vtxs);
      int * nds = new int[num_adj_vtxs];
      m3dc1_ent_getadj(&ent_dim, &eid, &vtx_dim, &nds[0], &num_adj_vtxs, &num_adj_vtxs);
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
          m3dc1_ent_getgeomclass(&vtx_dim, &nds[ii], &geom_class_dim, &geom_class_id);
          m3dc1_node_getnormvec(&nds[ii], normal1);
          m3dc1_node_getcurv(&nds[ii],&curv1);
          if(geom_class_dim==0)
            std::cout << "* node classified on geometric vertex with normal " << normal1[0] << " " << normal1[1] << " curv " << curv1 << std::endl;
          m3dc1_node_isongeombdry(&nds[ii+3], &is_bdy);
          assert(is_bdy);
          m3dc1_node_getnormvec(&nds[ii+3], normal2);
          m3dc1_node_getcurv(&nds[ii+3],&curv2);
          for(int jj = 0; jj < 2; jj++)
            assert(close<double>(normal1[jj], normal2[jj], 1e-6, 1e-6));
          assert(close<double>(curv1,curv2, 1e-6, 1e-6));
        }
      }
    }
  }
  // set/get field dof values
  int num_dofs = 1; //num_plane > 1 ? 12 : 6;
  int num_blks = 12;
  int agg_scp = m3dc1_field::PLANE_AGGREGATION;
  if(!PCU_Comm_Self())
    std::cout << "* creating fields with " << num_dofs
              << " dofs in " << num_blks << " blocks" << std::endl;
  m3dc1_field_create (&b_fld, "b_fld", &num_blks, &num_dofs, &agg_scp);
  m3dc1_field_create (&c_fld, "c_fld", &num_blks, &num_dofs, &agg_scp);
  m3dc1_field_create (&x_fld, "x_fld", &num_blks, &num_dofs, &agg_scp);
  //m3dc1_field_printcompnorm(&b_fld, "b_fld init info");
  if(!PCU_Comm_Self())
    std::cout << "* initialize b field ..." << std::endl;
  init_field(b_fld);
  if(!PCU_Comm_Self())
    std::cout << "* test field ops, set x = 3*b ..." << std::endl;
  test_field_ops(b_fld,x_fld);
  // test local matrix
  int lcl_mat_id = 1;
  int lcl_mat_tp = M3DC1_LOCAL_MAT;
  //m3dc1_matrix_create(&lcl_mat_id, &lcl_mat_tp, &scalar_tp, &b_fld);
  //test_matrix(lcl_mat_id);
  // test parallel matrix
  int par_mat_id = 2;
  int par_mat_tp = M3DC1_PARALLEL_MAT;
  m3dc1_matrix_create(&par_mat_id, &par_mat_tp, &scalar_tp, &b_fld);
  test_matrix(par_mat_id);
  m3dc1_matrix_reset(&lcl_mat_id);
  m3dc1_matrix_reset(&par_mat_id);
  // test again to ensure no issues will arise
  //  due to consecutive uses of the same matrices.
  test_matrix(lcl_mat_id);
  test_matrix(par_mat_id);
  /// cleanup
  m3dc1_matrix_delete(&lcl_mat_id);
  m3dc1_matrix_delete(&par_mat_id);
  m3dc1_field_delete(&x_fld);
  m3dc1_field_delete(&b_fld);
  m3dc1_field_delete(&c_fld);
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}
