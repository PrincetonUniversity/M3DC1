#include "m3dc1_scorec.h"
#include "m3dc1_field.h"
#include "m3dc1_matrix.h"
#include "m3dc1_mesh.h"
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
  double dof = 0.0;
  std::vector<double> dofs(dofs_per_nd * sz,1.0);
  for(int vtx_idx = 0; vtx_idx < num_lcl_vtx; ++vtx_idx)
  {
    int own = 0;
    m3dc1_ent_isowner(&vtx_dim, &vtx_idx,&own);
    if(own)
    {
      for(int dof_idx = 0; dof_idx < dofs_per_nd; ++dof_idx)
      {
        dofs[dof_idx] = dof;
        dof += 1.0;
      }
      m3dc1_ent_setdofdata(&vtx_dim, &vtx_idx, &fid, &dofs_per_nd, &dofs[0]);
    }
  }
  m3dc1_field_sync(&fid);
}
int main(int argc, char * argv[])
{
  // abort on invalid floating point operations and values
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  m3dc1_scorec_init(&argc,&argv);
  if (argc != 3 && !PCU_Comm_Self())
  {
    std::cout << "Usage: ./" << argv[0] << " [model_filename] [mesh_filename]" << std::endl;
    return M3DC1_FAILURE;
  }
  int num_plane = 2;
  int scalar_tp = M3DC1_REAL;
  m3dc1_model_setnumplane(&num_plane);
  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]);
  int zero = 0;
  m3dc1_mesh_build3d(&zero, &zero, &zero);
  // set/get field dof values
  int num_dofs = 1;
  int num_blks = 12;
  int agg_scps[] = {m3dc1_field::LOCAL_AGGREGATION,
                    m3dc1_field::PLANE_AGGREGATION,
                    m3dc1_field::GLOBAL_AGGREGATION};
  int fld_ids[] = {0, 1, 2, 3};
  // reverse num_dofs and num_blocks to create a naturally ordered field
  m3dc1_field_create(&fld_ids[0], "nat", &num_dofs, &num_blks, &agg_scps[0]);
  m3dc1_field_create(&fld_ids[1], "lcl", &num_blks, &num_dofs, &agg_scps[0]);
  m3dc1_field_create(&fld_ids[2], "pln", &num_blks, &num_dofs, &agg_scps[1]);
  m3dc1_field_create(&fld_ids[3], "gbl", &num_blks, &num_dofs, &agg_scps[2]);
  int dofs_per_nd = num_dofs * num_blks;
  if(!PCU_Comm_Self())
    std::cout << "* initialize fields ..." << std::endl;
  for(int ii = 0; ii < 4; ++ii)
    init_field(fld_ids[ii]);
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(fld_ids[0]);
  int mat_id = 1;
  int par_mat_tp = M3DC1_PARALLEL_MAT;
  m3dc1_matrix_create(&mat_id,&par_mat_tp,&scalar_tp,&fld_ids[0]);
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(mat_id);
  Vec * v = mat->getRHS();
  field2Vec(fld,*v);
  /*
  PetscViewer vwr;
  PetscViewerASCIIOpen(M3DC1_COMM_WORLD,"vec.m",&vwr);
  VecView(*v,vwr);
  PetscViewerDestroy(&vwr);
  */
  vec2Field(*v,fld);
  int vtx_dim = 0;
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  for(int fld_idx = 0; fld_idx < 4; ++fld_idx)
  {
    double dof_val = 0.0;
    std::vector<double> dofs(dofs_per_nd,-1.0);
    for(int vtx_idx = 0; vtx_idx < num_lcl_vtx; ++vtx_idx)
    {
      int own = 0;
      m3dc1_ent_isowner(&vtx_dim, &vtx_idx,&own);
      if(own)
      {
        m3dc1_ent_getdofdata(&vtx_dim, &vtx_idx, &fld_ids[fld_idx], &dofs_per_nd, &dofs[0]);
        for(int dof_idx = 0; dof_idx < dofs_per_nd; ++dof_idx)
        {
          assert(dofs[dof_idx] == dof_val);
          dof_val += 1.0;
        }
      }
    }
  }
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}
