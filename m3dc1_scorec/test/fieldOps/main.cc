#include "m3dc1_scorec.h"
#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include <fenv.h>
// "./av[0] [model_filename] [mesh_filename]
template <typename scalar>
bool close(scalar a, scalar b, scalar eps = 1e-8, scalar rel_eps = 1e-8)
{
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
  int dofs_per_nd = 0;
  int agg_scp = -1;
  m3dc1_field_getinfo(&fid,&num_blks,&dofs_per_nd,&agg_scp);
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
  int dofs_per_nd = 0;
  int agg_scp = -1;
  m3dc1_field_getinfo(&fid1,&num_blks,&dofs_per_nd,&agg_scp);
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
    assert(num_dofs_t == dofs_per_nd);
    m3dc1_ent_getdofdata(&vtx_dim, &ii, &fid2, &num_dofs_t, &dofs2.at(0));
    assert(num_dofs_t == dofs_per_nd);
    for(int jj = 0; jj < dofs_per_nd * sz; jj++)
      assert(close<double>(dofs2.at(jj), (scale + 1) * dofs1.at(jj), 1e-6, 1e-6));
  }
}
int main(int ac, char * av[])
{
  feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
  m3dc1_scorec_init(&ac,&av);
  int num_pln = 2;
  m3dc1_model_setnumplane(&num_pln);
  m3dc1_model_load(av[1]);
  m3dc1_model_print();
  m3dc1_mesh_load(av[2]);
  printStats(m3dc1_mesh::instance()->get_mesh());
  int zro = 0;
  m3dc1_mesh_build3d(&zro,&zro,&zro);
  int num_dofs = 12;
  int num_blks = 1;
  int agg_scp = m3dc1_field::GLOBAL_AGGREGATION;
  int fld1_id = 1;
  m3dc1_field_create(&fld1_id,"fld1",&num_blks,&num_dofs,&agg_scp);
  int fld2_id = 2;
  m3dc1_field_create(&fld2_id,"fld2",&num_blks,&num_dofs,&agg_scp);
  init_field(fld1_id);
  test_field_ops(fld1_id,fld2_id);
  return 0;
}
