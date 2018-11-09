/******************************************************************************
  (c) 2005-2018 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include "m3dc1_mesh.h"
#include "m3dc1_field.h"
#include <mpi.h>
#include <PCU.h>
#include "gmi_null.h" // FIXME: should be deleted later on since it's added temporarily for null model
#include <gmi_analytic.h>
#include <map>
#include "apfMDS.h"
#include "Expression.h"
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "pumi.h"
#ifdef M3DC1_TRILINOS
#include "m3dc1_ls.h"
#endif
#include <alloca.h>
double begin_mem = 0.0;
double begin_time = 0.0;
MPI_Comm M3DC1_COMM_WORLD = MPI_COMM_WORLD;
// helper routines
void synchronize_field(apf::Field* f);
bool m3dc1_double_isequal(double A, double B)
{
  double maxDiff = 1e-5;
  double maxRelDiff = 1e-5;
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
void m3dc1_scorec_setcomm_c(MPI_Comm cm)
{
  M3DC1_COMM_WORLD = cm;
}
void m3dc1_scorec_setcomm_f(MPI_Fint * cm)
{
  M3DC1_COMM_WORLD = MPI_Comm_f2c(*cm);
}
void m3dc1_scorec_init(int * argc, char ** argv[])
{
  int is_init = 0;
  MPI_Initialized(&is_init);
  if(!is_init)
    MPI_Init(argc,argv);
  pumi_start();
  las_init(argc,argv,M3DC1_COMM_WORLD);
  begin_time=MPI_Wtime();
}
void m3dc1_scorec_finalize()
{
  if (m3dc1_mesh::instance()->get_mesh()->getDimension()==3)
    m3dc1_ghost_delete();
  pumi_mesh_deleteGlobalID(m3dc1_mesh::instance()->get_mesh());  // delete global id
  m3dc1_mesh::instance()->clean(); // delete tag, field and internal data
  pumi_mesh_delete(m3dc1_mesh::instance()->get_mesh());
  if (!pumi_rank()) std::cout<<"\n* [M3D-C1] time: "<<MPI_Wtime()-begin_time<<" (sec)\n";
  PetscFinalize();
  pumi_finalize();
}
// plane functions
/*
void m3dc1_plane_setnum(int * num_plane)
{
  if (m3dc1_mesh::instance()->get_mesh()->getDimension()==3)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: #plane should be set before constructing 3D mesh\n";
    return;
  }
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane))
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: invalid #planes - mod(#ranks, #planes) should be 0\n";
    return;
  }
  m3dc1_model::instance()->set_num_plane(*num_plane);
}
*/
void m3dc1_plane_getnum(int * num_plane)
{
  *num_plane = m3dc1_model::instance()->num_plane;
}
void m3dc1_plane_getid(int* plane_id)
{
  *plane_id = m3dc1_model::instance()->local_planeid;
}
void m3dc1_plane_setphirange(double* min_val, double* max_val)
{
  if (m3dc1_mesh::instance()->get_mesh()->getDimension()==3)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: phi should be set before constructing 3D mesh\n";
    return;
  }
  m3dc1_model::instance()->set_phi(*min_val, *max_val);
}
void m3dc1_plane_setphi(int* planeid, double* phi)
{
  if (m3dc1_mesh::instance()->get_mesh()->getDimension()==3)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: phi should be set before constructing 3D mesh\n";
    return;
  }
  m3dc1_model::instance()->set_phi(*planeid, *phi);
}
void m3dc1_plane_getphi(int* planeid, double* phi)
{
  m3dc1_model::instance()->get_phi(*planeid, phi);
}
/** model functions */
void m3dc1_model_getplaneid(int* plane_id)
{
  *plane_id = m3dc1_model::instance()->local_planeid;
}
void m3dc1_model_getmincoord(double* x_min, double* y_min)
{
  *x_min = m3dc1_model::instance()->boundingBox[0];
  *y_min = m3dc1_model::instance()->boundingBox[1];
}
void m3dc1_model_getmaxcoord(double* x_max, double* y_max)
{
  *x_max = m3dc1_model::instance()->boundingBox[2];
  *y_max = m3dc1_model::instance()->boundingBox[3];
}
void m3dc1_model_load(const char * model_file)
{
  FILE *test_in = fopen (model_file,"r");
  if (!test_in)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: model file \""<<model_file<<"\" doesn't exist\n";
    return;
  }
  else
    fclose(test_in);
  pGeom g = pumi_geom_load(model_file, "analytic");
  m3dc1_model::instance()->model = g->getGmi();
  m3dc1_model::instance()->load_analytic_model(model_file);
  pumi_geom_freeze(g);
  m3dc1_model::instance()->caculateBoundingBox();
  // save the num of geo ent on the oringal plane
  m3dc1_model::instance()->numEntOrig[0]=m3dc1_model::instance()->model->n[0];
  //if (m3dc1_model::instance()->model->n[1]==1) assert(m3dc1_model::instance()->numEntOrig[0]==0); // for a smooth loop, there is no geo vtx
  m3dc1_model::instance()->numEntOrig[1]=m3dc1_model::instance()->model->n[1];
  m3dc1_model::instance()->numEntOrig[2]=m3dc1_model::instance()->model->n[2];
}
void m3dc1_model_print()
{
  if (PCU_Comm_Self() || m3dc1_model::instance()->local_planeid)
    return;
  gmi_iter* gf_it = gmi_begin(m3dc1_model::instance()->model, 2);
  gmi_ent* ge;
  while ((ge = gmi_next(m3dc1_model::instance()->model, gf_it)))
  {
    gmi_set* gf_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())    std::cout<<"model face id "<<gmi_tag(m3dc1_model::instance()->model, ge)
                                      <<" - # adj model edges: "<<gf_edges->n<<"\n";
    if (gf_edges->n)
    {
      if (!PCU_Comm_Self())      std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gf_edges->n; ++i)
        if (!PCU_Comm_Self())        std::cout<<gmi_tag(m3dc1_model::instance()->model,  gf_edges->e[i])<<" ";
      if (!PCU_Comm_Self())      std::cout<<"\n";
    }
    gmi_free_set(gf_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gf_it);
// verify geom edge
  gmi_iter* ge_it = gmi_begin(m3dc1_model::instance()->model, 1);
  while ((ge = gmi_next(m3dc1_model::instance()->model, ge_it)))
  {
    M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,ge);
    if (!pn)
    {
      if (!PCU_Comm_Self())       std::cout<<"["<<PCU_Comm_Self()<<"] model edge "<<gmi_tag(m3dc1_model::instance()->model, ge)
                                           <<" - failed with gmi_analytic_data retrieval\n";
    }
  }
  gmi_end(m3dc1_model::instance()->model, ge_it);
// verify gv_edges
  gmi_iter* gv_it = gmi_begin(m3dc1_model::instance()->model, 0);
  while ((ge = gmi_next(m3dc1_model::instance()->model, gv_it)))
  {
    gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())
      std::cout<<"model vertex id "<<gmi_tag(m3dc1_model::instance()->model, ge)
               <<" - # adj model edges: "<<gv_edges->n<<"\n";
    if (gv_edges->n)
    {
      if (!PCU_Comm_Self())
        std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gv_edges->n; ++i)
        if (!PCU_Comm_Self()) std::cout<<gmi_tag(m3dc1_model::instance()->model,  gv_edges->e[i])<<" ";
      if (!PCU_Comm_Self()) std::cout<<"\n";
    }
    gmi_free_set(gv_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gv_it);
}
void m3dc1_model_setnumplane(int* num_plane)
{
  if (*num_plane < 1 || PCU_Comm_Peers() % (*num_plane))
    return;
  m3dc1_model::instance()->set_num_plane(*num_plane);
}
void m3dc1_model_getnumplane(int* num_plane)
{
  *num_plane = m3dc1_model::instance()->num_plane;
}
/** mesh functions */
#include <parma.h>
void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim)
{
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    m->setDoubleTag(e, tag, &w);
  m->end(it);
}
apf::MeshTag* setWeights(apf::Mesh* m)
{
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}
void clearTags(apf::Mesh* m, apf::MeshTag* t)
{
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}
#include <sstream>
void m3dc1_mesh_load(const char * mesh_file)
{
  std::string in(mesh_file);
  size_t p = in.rfind('.');
  std::string base = in.substr(0,p);
  std::string ext = in.substr(p+1);
  std::stringstream s;
  s << base << pumi_rank()<<"." << ext;
  std::string partFile = s.str();
  FILE* test_in = fopen (partFile.c_str(),"r");
  if (!test_in)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: mesh file \""<<mesh_file<<"\" doesn't exist\n";
    return;
  }
  else
    fclose(test_in);
  if (m3dc1_model::instance()->local_planeid == 0) // master plane
  {
    m3dc1_mesh::instance()->load_mesh(mesh_file);
    /* vertex load balancing */
    //Parma_PrintPtnStats(m3dc1_mesh::instance()->mesh, "initial");
    // clean-up tag, field and numbering loaded from file
    apf::Mesh2 * mesh = m3dc1_mesh::instance()->get_mesh();
    while (mesh->countFields())
    {
      apf::Field* f = mesh->getField(0);
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": field "<<getName(f)<<" deleted\n";
      destroyField(f);
    }
    while(mesh->countNumberings())
    {
      apf::Numbering* n = mesh->getNumbering(0);
      if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
      destroyNumbering(n);
    }
    apf::DynamicArray<apf::MeshTag*> tags;
    mesh->getTags(tags);
    for (size_t ii=0; ii<tags.getSize(); ++ii)
    {
      if (mesh->findTag("norm_curv")==tags[ii]) continue;
      for (int idim=0; idim<4; idim++)
        apf::removeTagFromDimension(mesh, tags[ii], idim);
      mesh->destroyTag(tags[ii]);
    }
  }
  else // non-master plane
    m3dc1_mesh::instance()->create_mesh();
  m3dc1_mesh::instance()->initialize();
  if (m3dc1_model::instance()->num_plane==1) // 2D problem
  {
    compute_globalid(m3dc1_mesh::instance()->get_mesh(), 0);
    compute_globalid(m3dc1_mesh::instance()->get_mesh(), 2);
  }
}
Ghosting* getGhostingPlan(apf::Mesh2* m)
{
  int mdim=pumi_mesh_getDim(m);
  Ghosting* plan = new Ghosting(m, mdim);
  apf::MeshIterator* it = m->begin(mdim);
  apf::MeshEntity* e;
  int pid=m3dc1_model::instance()->prev_plane_partid;
  while ((e = m->iterate(it)))
    plan->send(e, pid);
  m->end(it);
  return plan;
}
void create_backward_plane_ghosting(apf::Mesh2* m)
{
  Ghosting* ghosting_plan = getGhostingPlan(m);
  // ghosting_plan destroyed after ghosting is complete
  pumi_ghost_create(m, ghosting_plan);
  for (int i=0; i<m->countFields();++i)
    synchronize_field(m->getField(i)); // synch fields based on m3dc1's ownership rule
}
void m3dc1_mesh_build3d (int* num_field, int* field_id,
                         int* num_dofs_per_value)
{
  // switch COMM to GLOBAL COMM
  MPI_Comm groupComm = PCU_Get_Comm();
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);
  // initialize phi value and construct 3d
  int num_plane = m3dc1_model::instance()->num_plane;
  if (num_plane==1)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 ERROR] "<<__func__<<" failed: #plane must be greater than 1\n";
    return;
  }
  if (!m3dc1_model::instance()->phi)
  {
    if (!PCU_Comm_Self())
      std::cout<<"[M3D-C1 INFO] "<<__func__
               <<": setting phi range (0.0, "<<2.0*M3DC1_PI/num_plane*(num_plane-1)<<")\n";
    m3dc1_model::instance()->set_phi(0.0, 2.0*M3DC1_PI/num_plane*(num_plane-1));
  }
  m3dc1_model::instance()->setupCommGroupsPlane();
  // returns error if num_plane>1
  pumi_mesh_deleteGlobalID(m3dc1_mesh::instance()->get_mesh());
  m3dc1_mesh::instance()->build3d(*num_field, field_id, num_dofs_per_value);
  // update global ID
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), 0);
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), 3);
  // added as per Bill's request -- 03/14/18
  create_backward_plane_ghosting(m3dc1_mesh::instance()->get_mesh());
#ifdef DEBUG
  int nlocal = m3dc1_mesh::instance()->get_mesh()->count(3);
  int nghost = nlocal - m3dc1_mesh::instance()->get_local_count(3);
  if (!pumi_rank()) std::cout<<"# elements "<<nlocal<<", #ghost elements "<<nghost<<"\n";
  printStats(m3dc1_mesh::instance()->get_mesh());
#endif
}
/* ghosting functions */
void m3dc1_ghost_create (int* num_layer )
{
  pMesh m = m3dc1_mesh::instance()->get_mesh();
  // create layer
  pumi_ghost_createLayer (m, 0, pumi_mesh_getDim(m), *num_layer, 1);
  for (int i=0; i<m->countFields();++i)
    synchronize_field(m->getField(i)); // synch fields based on m3dc1's ownership rule
  if (!pumi_rank()) std::cout<<"[M3D-C1 INFO] "<<* num_layer<<" ghost layer(s) created\n";
}
void m3dc1_ghost_delete()
{
  pumi_ghost_delete (m3dc1_mesh::instance()->get_mesh());
  if (!pumi_rank()) std::cout<<"[M3D-C1 INFO] ghost layer(s) deleted\n";
}
void m3dc1_mesh_getnument (int* /* in*/ edim, int* num_ent)
{
  *num_ent = m3dc1_mesh::instance()->get_mesh()->count(*edim);
}
void m3dc1_mesh_getnumownent (int* /* in*/ edim, int* num_ent)
{
  *num_ent = m3dc1_mesh::instance()->get_own_count(*edim);
}
void m3dc1_mesh_getnumglobalent (int* /* in*/ edim, int* num_ent)
{
  *num_ent = m3dc1_mesh::instance()->get_global_count(*edim);
}
void m3dc1_mesh_getnumghostent (int* /* in*/ edim, int* num_ent)
{
  if (*edim<0 || *edim > 3)
    return;
  *num_ent = m3dc1_mesh::instance()->get_mesh()->count(*edim) -
    m3dc1_mesh::instance()->get_local_count(*edim);
}
void m3dc1_mesh_search(int* initial_simplex,
                       double* final_position,
                       int* final_simplex)
{
  bool located = false;
  apf::MeshEntity* e = NULL;
  apf::MeshEntity* simplex = NULL;
  apf::Mesh2* m = m3dc1_mesh::instance()->get_mesh();
  apf::Adjacent adjacent;
  int edge_curr_index, edge_prev_index;
  int simplex_dim = m->getDimension();
  int vertex_dim = 0, edge_dim = 1;
  int bridge_dim = simplex_dim - 1;
  int count = 0, max_count = m->count(simplex_dim);
  double tol = 1e-15;
  simplex = apf::getMdsEntity(m, simplex_dim, *initial_simplex);
  while (not located)
  {
    int simplex_index = apf::getMdsIndex(m, simplex);
    apf::Downward vertices;
    apf::Vector3 v_coord[3];
    apf::Matrix<3, 3> A;
    int nv = m->getDownward(simplex, vertex_dim, vertices);
    for (int j = 0; j < nv; ++j) {
      m->getPoint(vertices[j], 0, v_coord[j]);
      // Note: second argument is 0 for linear meshes
    }
    // Compute (linear) barycentric coordinates of final position
    // with respect to current simplex
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 3; ++k)
        A[j][k] = v_coord[k][j];
    for (int j = 0; j < 3; ++j)
      A[2][j] = 1;
    apf::Matrix<3, 3> Ainv = apf::invert(A);
    apf::Vector3 b_coords; // = Ainv * final_position
    for (int j = 0; j < 3; ++j) {
      b_coords[j] = (Ainv[j][0] * final_position[0] +
                     Ainv[j][1] * final_position[1] +
                     Ainv[j][2]);
    }
    // If all positive for current simplex, exit.
    if (((b_coords[0] >= -tol) && (b_coords[0] <= (1 + tol))) &&
        ((b_coords[1] >= -tol) && (b_coords[1] <= (1 + tol))) &&
        ((b_coords[2] >= -tol) && (b_coords[2] <= (1 + tol)))) {
      located = true;
      *final_simplex = simplex_index;
      break;
    }
    // Otherwise, check which coordinates are negative
    //bool b_negative[3] = {false, false, false};
    int bneg_index[2] = {0, 0};
    int bneg_count = 0;
    for (int j = 0; j < 3; ++j)
      if (b_coords[j] < 0)
      {
        //b_negative[j] = true;
        bneg_index[bneg_count] = j;
        ++bneg_count;
      }
    // Obtain the index of most negative coordinate
    assert(bneg_count > 0 && bneg_count < 3);
    // Ensure bneg_index[0] is the index of vertex whose corresponding
    // barycentric coordinate is most negative; ties automatically resolved
    // as a result.
    if (bneg_count == 2)
      if (fabs(b_coords[bneg_index[0]]) <
          fabs(b_coords[bneg_index[1]])) {
        int tmp_index = bneg_index[0];
        bneg_index[0] = bneg_index[1];
        bneg_index[1] = tmp_index;
      }
    // Determine edge opposite to this vertex
    apf::MeshEntity* edge_vertices[2];
    int edge_count = 0;
    int edge_type = 1;    // Mesh entity type; see APF documentation.
    for (int j = 0; j < 3; ++j)
      if (j != bneg_index[0]) {
        edge_vertices[edge_count] = vertices[j];
        ++edge_count;
      }
    // Find neighboring simplex sharing this edge
    e = apf::findElement(m, edge_type, edge_vertices);
    edge_curr_index = apf::getMdsIndex(m, e);
    // If current edge choice is same as previous edge, pick edge
    // opposite to second least (actual, not absolute, valued) barycentric
    // coordinate.
    if (edge_curr_index == edge_prev_index) {
      edge_count = 0;
      for (int j = 0; j < 3; ++j)
        if (j != bneg_index[1]) {
          edge_vertices[edge_count] = vertices[j];
          ++edge_count;
        }
      e = apf::findElement(m, edge_type, edge_vertices);
      edge_curr_index = apf::getMdsIndex(m, e);
    }
    apf::getBridgeAdjacent(m, e,
                           bridge_dim, simplex_dim,
                           adjacent);
    if (adjacent.getSize() == 2) {
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
        int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
        if (new_simplex_index != simplex_index)
          simplex = adjacent[j];
      }
    }
    else {
      apf::Downward edges;
      int ne = m->getDownward(simplex, edge_dim, edges);
      for (int jj = 0; jj < ne; ++jj) {
        int edge_tmp_index = apf::getMdsIndex(m, edges[jj]);
        if ((edge_tmp_index != edge_curr_index) &&
            (edge_tmp_index != edge_prev_index)) {
          e = edges[jj];
          break;
        }
      }
      edge_curr_index = apf::getMdsIndex(m, e);
      apf::getBridgeAdjacent(m, e,
                             bridge_dim, simplex_dim,
                             adjacent);
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
        int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
        if (new_simplex_index != simplex_index)
          simplex = adjacent[j];
      }
    }
    // Keep track of edge via which we entered the current simplex
    edge_prev_index = edge_curr_index;
    ++count;
    if (count == max_count){
      // std::cout << "\nError: hit maximum number of simplices to look for.";
      *final_simplex = -2;
      return;
    }
  }
}
void m3dc1_mesh_write(const char* filename, int* option)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->get_mesh();
  apf::MeshEntity* e;
  if (*option==0 ||*option==3)
  {
    int dim=2, num_ent=m3dc1_mesh::instance()->get_mesh()->count(2);
    vector<double> geoId (num_ent);
    apf:: MeshIterator* it = mesh->begin(dim);
    while ((e = mesh->iterate(it)))
    {
      int eid = getMdsIndex(m3dc1_mesh::instance()->get_mesh(), e);
      int geom_class_dim,geom_class_id;
      m3dc1_ent_getgeomclass (&dim, &eid, &geom_class_dim, &geom_class_id);
      geoId.at(eid)=geom_class_id;
    }
    mesh->end(it);
    apf::writeVtkFiles(filename,m3dc1_mesh::instance()->get_mesh());
    int one=1;
    if (*option==3) output_face_data (&one, &geoId[0], "geoId");
    /*apf::removeTagFromDimension(mesh, tag, dim);
      mesh->destroyTag(tag);*/
  }
  else
  {
    char filename_buff[256];
    sprintf(filename_buff, "%s.smb",filename);
    int fid = 12; // FIXME: why is the field ID hard-coded?
    double dofBuff[1024];
    m3dc1_field * fld = m3dc1_mesh::instance()->get_field(fid);
    apf::Field * f = fld->getCoreField();
    int numDof = countComponents(f);
    apf::MeshTag* tag = mesh->createDoubleTag("field12", numDof);
    apf:: MeshIterator* it = mesh->begin(0);
    while ((e = mesh->iterate(it)))
    {
      get_ent_dofdata(fld, e, dofBuff);
      mesh->setDoubleTag(e, tag, dofBuff);
    }
    mesh->end(it);
    m3dc1_mesh::instance()->get_mesh()->writeNative(filename_buff);
    apf::removeTagFromDimension(mesh, tag, 0);
    mesh->destroyTag(tag);
  }
}
void m3dc1_field_verify()
{
  m3dc1_mesh::instance()->verify_fields();
}
/* mesh entity functions */
void m3dc1_ent_getglobalid (int * edim, int * eid, int * global_eid)
{
  apf::MeshEntity * e = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  assert(e);
  *global_eid = get_ent_globalid(m3dc1_mesh::instance()->get_mesh(), e);
}
void m3dc1_ent_getgeomclass (int*  edim, int*  eid,
                             int*  geom_class_dim, int*  geom_class_id)
{
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  assert(ent);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->get_mesh()->toModel(ent));
  *geom_class_dim = gmi_dim(m3dc1_model::instance()->model,gent);
  *geom_class_id = gmi_tag(m3dc1_model::instance()->model,gent);
  // if 3D mesh, need to return the classification on the original plane
  if ( m3dc1_mesh::instance()->get_mesh()->getDimension() ==3 )
  {
    int numEntOrig[3];
    memcpy( numEntOrig, m3dc1_model::instance()->numEntOrig, sizeof(numEntOrig));
    *geom_class_id-=1;
    switch (*geom_class_dim)
    {
    case 3: *geom_class_id%=numEntOrig[2]; break;
    case 2: *geom_class_id%=(numEntOrig[1]+numEntOrig[2]); break;
    case 1:  *geom_class_id%=(numEntOrig[0]+numEntOrig[1]); break;
    case 0: *geom_class_id%=(numEntOrig[0]);
    }
    *geom_class_id+=1;
  }
}
void m3dc1_ent_getadj (int*  edim, int*  eid,
                       int*  adj_dim, int*  adj_ent,
                       int*  adj_ent_allocated_size, int*  num_adj_ent)
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  if (!e || *adj_dim==*edim)
    return;
  if (*adj_dim>*edim) // upward
  {
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->get_mesh()->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
    if (*adj_ent_allocated_size<*num_adj_ent)
      return;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m3dc1_mesh::instance()->get_mesh(), adjacent[i]);
  }
  else if (*adj_dim<*edim)
  {
    apf::Downward downward;
    *num_adj_ent = m3dc1_mesh::instance()->get_mesh()->getDownward(e, *adj_dim, downward);
    if (*adj_ent_allocated_size<*num_adj_ent)
      return;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m3dc1_mesh::instance()->get_mesh(), downward[i]);
    //adjust the order to work with m3dc1
    if (m3dc1_mesh::instance()->get_mesh()->getDimension()==3 && *edim==3 &&*adj_dim==0 &&adj_ent[0]>adj_ent[3])
    {
      int buff[3];
      memcpy(buff, adj_ent, 3*sizeof(int));
      memcpy(adj_ent, adj_ent+3, 3*sizeof(int));
      memcpy(adj_ent+3, buff, 3*sizeof(int));
    }
  }
}
void m3dc1_ent_getnumadj (int*  edim, int*  eid,
                          int*  adj_dim, int*  num_adj_ent)
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  if (!e || *adj_dim==*edim)
    return;
  if (*adj_dim>*edim) // upward
  {
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->get_mesh()->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
  }
  else if (*adj_dim<*edim)
  {
    apf::Downward downward;
    *num_adj_ent = m3dc1_mesh::instance()->get_mesh()->getDownward(e, *adj_dim, downward);
  }
}
void m3dc1_ent_getownpartid (int*  edim, int*  eid,
                             int*  owning_partid)
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  assert(e);
  *owning_partid = get_ent_ownpartid(e);
}

void m3dc1_ent_isowner(int * edim,
                       int * eid,
                       int * ismine)
{
  apf::Mesh2 * m = m3dc1_mesh::instance()->get_mesh();
  apf::MeshEntity * e = getMdsEntity(m, *edim, *eid);
  assert(e);
  *ismine = (is_ent_original(e) ? 1 : 0);
}
void m3dc1_ent_isghost(int*  edim, int*  eid, int* isghost)
{
  pMeshEnt e =getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  if (pumi_ment_isGhost(e))
    *isghost=1;
  else
    *isghost=0;
}
// node-specific functions
void m3dc1_node_getcoord (int*  node_id, double*  coord)
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, *node_id);
  assert(e);
  apf::Vector3 xyz;
  m3dc1_mesh::instance()->get_mesh()->getPoint(e, 0, xyz);
  for (int i=0; i<3; ++i)
    coord[i] = xyz[i];
}
void m3dc1_node_getglobalid (int * edim, int*  eid, int*  global_eid)
{
  (void)edim;
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, *eid);
  assert(e);
  *global_eid = get_ent_globalid(m3dc1_mesh::instance()->get_mesh(), e);
}
void get_gv_edges(gmi_ent* gvertex, std::vector<gmi_ent*>& gedges)
{
  gedges.clear();
  gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, gvertex, 1);
  for (int i=0; i<gv_edges->n; ++i)
    gedges.push_back(gv_edges->e[i]);
  gmi_free_set(gv_edges);
  assert(gedges.size()>=1);
}
void m3dc1_node_getnormvec(int * vtx, double * xyzt)
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, *vtx);
  xyzt[2] = 0.0;
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->get_mesh()->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if (gType !=  1 && gType !=  0)
  {
    xyzt[0] = xyzt[1] = 0.0;
    return;
  }
  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->get_mesh()->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->get_mesh()->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->get_mesh()->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    xyzt[0] = norm_curv[0];
    xyzt[1] = norm_curv[1];
    return;
  }
  else
  { // if norm/curv is not attached, evaluate
    apf::Vector3 param(0,0,0);
    m3dc1_mesh::instance()->get_mesh()->getParam(vt,param);
    // geo node avage on the connected edges
    if (gType == 0) // node is on the
    {
      apf::Vector3 vcd_t;
      double vcd[3];
      m3dc1_mesh::instance()->get_mesh()->getPoint(vt, 0, vcd_t);
      for (int ii=0; ii<3; ++ii)
        vcd[ii]=vcd_t[ii];
      std::vector<gmi_ent*> gEdges;
      get_gv_edges(gent, gEdges);
      int numEdgePlane=0;
      double normalvec[3]={0.,0.,0.};
      xyzt[0]=xyzt[1]=xyzt[2]=0;
      if (gEdges.size()<2)
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: #adjEdge of gVertex="<<gEdges.size()<<" (it should be minimum 2) \n";
      assert(gEdges.size()>=2);
      for (size_t ii=0;ii<gEdges.size();++ii)
      {
        gmi_ent* pe = gEdges.at(ii);
        double cd[3]={0,0,0};
        double paraRange[2];
        gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
        M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
        if (!pn) continue;
        numEdgePlane++;
        M3DC1::evalCoord(paraRange[0],cd, pn);
        if (checkSamePoint2D(vcd,cd))
        {
          M3DC1::evalNormalVector(pn[0],pn[1], paraRange[0], normalvec);
        }
        else
        {
          evalNormalVector(pn[0],pn[1], paraRange[1], normalvec);
        }
        xyzt[0]+=normalvec[0];
        xyzt[1]+=normalvec[1];
        xyzt[2]+=normalvec[2];
      }
      if (numEdgePlane!=2)
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: numEdgePlane="
                 <<numEdgePlane<<" (it should be 2) \n";
      assert(numEdgePlane==2);
      double arclen=sqrt(xyzt[0]*xyzt[0]+xyzt[1]*xyzt[1]+xyzt[2]*xyzt[2]);
      assert(arclen>0);
      xyzt[0]=xyzt[0]/arclen;
      xyzt[1]=xyzt[1]/arclen;
      xyzt[2]=xyzt[2]/arclen;
    }
    else
    {
      apf::Vector3 param(0,0,0);
      m3dc1_mesh::instance()->get_mesh()->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalNormalVector(pn[0],pn[1], param[0], xyzt);
    }
  }
}
void m3dc1_node_getcurv (int*  vtx, double*  curv)
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, *vtx);
  assert(vt);
  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->get_mesh()->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->get_mesh()->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->get_mesh()->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    *curv = norm_curv[2];
    return;
  }
  *curv=0.0;
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->get_mesh()->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if (gType==0)
  {
    apf::Vector3 vcd_t;
    double vcd[3];
    m3dc1_mesh::instance()->get_mesh()->getPoint(vt, 0, vcd_t);
    for (int ii=0; ii<3; ++ii)
      vcd[ii]=vcd_t[ii];
    std::vector<gmi_ent*> gEdges;
    get_gv_edges(gent, gEdges);
    int numEdgesPlane=0;
    double curv_tmp;
    for (size_t ii = 0; ii < gEdges.size(); ++ii)
    {
      gmi_ent* pe = gEdges.at(ii);
      double cd[3]={0,0,0};
      double paraRange[2];
      gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
      if (!pn) continue;
      numEdgesPlane++;
      evalCoord(paraRange[0], cd, pn);
      if (checkSamePoint2D(vcd,cd))
      {
        evalCurvature(pn[0],pn[1], paraRange[0], &curv_tmp);
      }
      else
      {
        evalCurvature(pn[0],pn[1], paraRange[1], &curv_tmp);
      }
      *curv+=curv_tmp;
    }
    assert(numEdgesPlane==2);
    *curv/=numEdgesPlane;
  }
  else if (gType==1)
  {
    apf::Vector3 param(0,0,0);
    m3dc1_mesh::instance()->get_mesh()->getParam(vt,param);
    M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
    evalCurvature(pn[0],pn[1], param[0], curv);
  }
}
void m3dc1_node_isongeombdry (int*  vtx, int*  on_geom_bdry)
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, *vtx);
  assert(vt);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->get_mesh()->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  *on_geom_bdry=(gType==0||gType==1);
}
//=========================================================================
void write_node(apf::Mesh2* m, const char* filename, int start_index)
{
  char node_filename[256];
  sprintf(node_filename,"%s-%d", filename, PCU_Comm_Self());
  FILE * np =fopen(node_filename, "w");
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    apf::Vector3 xyz;
    m->getPoint(e, 0, xyz);
    fprintf(np, "%d\t%lf\t%lf\t%lf\n", get_ent_globalid(m,e)+start_index, xyz[0],xyz[1],xyz[2]);
  } // while
  m->end(it);
  fclose(np);
}
void m3dc1_node_write (const char* filename, int* start_index)
{
  write_node(m3dc1_mesh::instance()->get_mesh(), filename, *start_index);
}
// region-specific function
void m3dc1_region_getoriginalface(int*  elm, int*  fac)
{
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 3, *elm);
  apf::Downward downward;
  int num_adj_ent = m3dc1_mesh::instance()->get_mesh()->getDownward(ent, 2, downward);
  assert(num_adj_ent==5);
  int triFace[2];
  int counter=0;
  for (int i=0; i<num_adj_ent; ++i)
  {
    apf::Downward downward2;
    int num_edge = m3dc1_mesh::instance()->get_mesh()->getDownward(downward[i], 1, downward2);
    if (num_edge==3) triFace[counter++]= getMdsIndex(m3dc1_mesh::instance()->get_mesh(),downward[i]);
  }
  assert(counter==2);
  *fac = std::min(triFace[0],triFace[1]);
}
void m3dc1_field_create(FieldID * fid,
                        const char * fnm,
                        int * blks_per_nd,
                        int * dofs_per_blk,
                        int * agg)
{
  auto fld = new m3dc1_field(*fid,
                             fnm,
                             *blks_per_nd,
                             *dofs_per_blk,
                             (m3dc1_field::aggregation_scope)(*agg));
  m3dc1_mesh::instance()->add_field(*fid,fld);
}
void m3dc1_field_delete (FieldID * fid)
{
  m3dc1_mesh::instance()->destroy_field(*fid);
}
const char * m3dc1_field_getname(FieldID * fid)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  const std::string & fnm = fld->getName();
  return fnm.c_str();
}
void m3dc1_field_getinfo(FieldID * fid,
                         int * blks_per_nd,
                         int * dofs_per_blk,
                         int * agg_scp)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  *blks_per_nd = fld->getBlocksPerNode();
  *dofs_per_blk = fld->getDofsPerBlock();
  *agg_scp = fld->getAggregationScope();
}
void m3dc1_field_exist(FieldID * fid, int * exists)
{
  *exists = m3dc1_mesh::instance()->field_exists(*fid);
}
void synchronize_field(apf::Field* f)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->get_mesh();
  apf::MeshEntity* e;
  int n = countComponents(f);
  double* sender_data = new double[n];
  double* dof_data = new double[n];
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(e) || (!m->isShared(e)&&!m->isGhosted(e)))
      continue;
    getComponents(f, e, 0, dof_data);
    if (m->isShared(e))
    {
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(apf::Copies,remotes,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
      }
    }
    if (m->isGhosted(e))
    {
      apf::Copies ghosts;
      m->getGhosts(e,ghosts);
      APF_ITERATE(apf::Copies,ghosts,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
      }
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      // FIXME based on new m3dc1_field
      setComponents(f, r, 0, sender_data);
    }
  delete [] dof_data;
  delete [] sender_data;
}
void m3dc1_field_sync(FieldID *  fid)
{
  synchronize_field(m3dc1_mesh::instance()->get_field(*fid)->getCoreField());
}
// send non-owned copies' dof to owner copy and add them up
void accumulate_field(apf::Field* f)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->get_mesh();
  apf::MeshEntity* e;
  int own_partid, n = countComponents(f);
  double* dof_data = new double[n];
  double* sender_data = new double[n];
  apf::MeshEntity* own_e;
  apf::MeshEntity* r;
  std::map<apf::MeshEntity*, std::vector<double> > save_map;
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid=get_ent_ownpartid(e);
    if (own_partid==PCU_Comm_Self() || pumi_ment_isGhost(e)) continue;
    assert(m->isShared(e));
    own_e = get_ent_owncopy(e);
    getComponents(f, e, 0, &(dof_data[0]));
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(dof_data[0]),n*sizeof(double));
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while (! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      for (int i = 0; i < n; ++i)
        save_map[r].push_back(sender_data[i]);
    }
  for (std::map<apf::MeshEntity*, std::vector<double> >::iterator mit=save_map.begin(); mit!=save_map.end(); ++mit)
  {
    e = mit->first;
    getComponents(f, e, 0, dof_data);
    int num_data = mit->second.size()/n;
    for (int i=0; i<num_data;++i)
    {
      for (int j=0; j<n; ++j)
        dof_data[j] += mit->second[i*n+j];
    }
    // FIXME based on new m3dc1_field
    setComponents(f, e, 0, dof_data);
  }
  delete [] dof_data;
  delete [] sender_data;
}
void m3dc1_field_sum(FieldID *  fid)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  apf::Field * cfld = fld->getCoreField();
  accumulate_field(cfld);
  synchronize_field(cfld);
}
// todo : replace with a FrozenFieldOp / FieldOp
void m3dc1_field_sumsq(FieldID *  fid, double *  sum)
{
  *sum = 0.0;
  m3dc1_mesh * m3msh = m3dc1_mesh::instance();
  apf::Mesh2 * msh = m3msh->get_mesh();
  m3dc1_field * m3fld = m3msh->get_field(*fid);
  apf::Field * fld = m3fld->getCoreField();
  int num_dof = countComponents(fld);
  double * dof_data = new double[num_dof];
  apf::MeshEntity * e = NULL;
  apf::MeshIterator * it = msh->begin(0);
  while ((e = msh->iterate(it)))
  {
    if(msh->isOwned(e))
      continue;
    getComponents(fld, e, 0, dof_data);
    for(int ii = 0; ii < num_dof; ++ii)
      *sum += dof_data[ii] * dof_data[ii];
  }
  msh->end(it);
  delete [] dof_data;
}
// todo : calculate and reduce the local sum in lieu of calculating the entire sum over the communicator
void m3dc1_field_sum_plane (FieldID *  fid)
{
  int value_type = 0;
#ifdef PETSC_USE_COMPLEX
  value_type = 1;
#endif
  m3dc1_mesh * m3msh = m3dc1_mesh::instance();
  apf::Mesh2 * msh = m3msh->get_mesh();
  int num_lcl_vtx = msh->count(0);
  m3dc1_field * m3fld = m3msh->get_field(*fid);
  int num_dof = m3fld->getDofsPerNode() * num_lcl_vtx;
  MPI_Comm pcm = m3dc1_model::instance()->getMPICommPlane();
  int data_size = num_dof * (value_type+1);
  double * fld_arr = NULL;
  m3dc1_field_getdataptr(fid, &fld_arr);
  double * sendbuf = new double [data_size];
  m3dc1_field_retrieve(fid, sendbuf, &num_dof);
  MPI_Allreduce(sendbuf,fld_arr,data_size,MPI_DOUBLE,MPI_SUM,pcm);
  synchronize_field(m3fld->getCoreField());
  delete [] sendbuf;
}
void m3dc1_field_printcompnorm(FieldID*  fid, const char * info)
{
  double * pts = NULL;
  m3dc1_field_getdataptr(fid, &pts);
  int num_dof = 0;
  m3dc1_field_getnumlocaldof(fid, &num_dof);
  m3dc1_field * m3fld = m3dc1_mesh::instance()->get_field(*fid);
  if (m3fld->getBlocksPerNode()>1)
  {
    if(!PCU_Comm_Self()) std::cout<<"[M3D-C1 WARNING] "<<__func__<<": not implemented for multiple valued DOF's\n";
    return;
  }
  apf::Field* f = m3fld->getCoreField();
  int dof_per_node = countComponents(f);
  int num_comp=dof_per_node/C1TRIDOFNODE;
  vector<double> norms(dof_per_node/C1TRIDOFNODE);
  int j=0;
  for(int i=0; i<num_dof/C1TRIDOFNODE; ++i)
  {
    for(int k=0; k<6; k++)
      norms.at(j)+=pts[i*C1TRIDOFNODE+k]*pts[i*C1TRIDOFNODE+k];
    ++j;
    j%=num_comp;
  }
  vector<double> buff=norms;
  MPI_Allreduce(&buff[0],&norms[0], num_comp, MPI_DOUBLE,MPI_SUM,M3DC1_COMM_WORLD);
  int psize;
  MPI_Comm_size(M3DC1_COMM_WORLD,&psize);
  if (PCU_Comm_Self() == psize-1)
  {
    std::cout<< "norm of vec "<<info;
    for(int i=0; i<num_comp; ++i)
      std::cout<<" "<<std::sqrt(norms[i]);
    std::cout<<std::endl;
  }
}
void m3dc1_field_getlocaldofid (FieldID * fid, int * start_dof_id, int * end_dof_id_plus_one)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int dof_per_nd = fld->getDofsPerNode();
  *start_dof_id = 0;
  *end_dof_id_plus_one = dof_per_nd * m3dc1_mesh::instance()->get_mesh()->count(0); // assumes only verts have nodes
}
void m3dc1_field_getglobaldofid(FieldID * fid,
                                int * start_dof_id,
                                int * end_dof_id_plus_one)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int dof_per_nd = fld->getDofsPerNode();
  *start_dof_id = 0;
  *end_dof_id_plus_one = dof_per_nd * m3dc1_mesh::instance()->get_global_count(0); // assumes only verts have nodes
}
void m3dc1_field_getowndofid (FieldID * fid, int * start_dof_id, int * end_dof_id_plus_one)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int num_own_ent = m3dc1_mesh::instance()->get_own_count(0); // assumes only verts have nodes
  int dof_per_nd = fld->getDofsPerNode();
  int start_id = num_own_ent;
  PCU_Exscan_Ints(&start_id,1);
  *start_dof_id = start_id * dof_per_nd;
  *end_dof_id_plus_one = *start_dof_id + num_own_ent * dof_per_nd;
}
void m3dc1_field_getnumglobaldof(FieldID * fid, int * num_global_dof)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  *num_global_dof = (msh->get_global_count(0))*fld->getDofsPerNode();
}
void m3dc1_field_getnumlocaldof(FieldID * fid, int * num_local_dof)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  *num_local_dof = (msh->get_mesh()->count(0)) * fld->getDofsPerNode();
}
void m3dc1_field_getnumowndof (FieldID * fid, int * num_own_dof)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  *num_own_dof = (msh->get_own_count(0)) * fld->getDofsPerNode();
}
void m3dc1_field_getnumghostdof (FieldID* fid, int*  num_ghost_dof)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  int num_vtx = msh->get_mesh()->count(0) - msh->get_local_count(0);
  *num_ghost_dof = num_vtx * fld->getDofsPerNode();
}
void m3dc1_field_getdataptr (FieldID* fid, double** pts)
{
  apf::Field * f = m3dc1_mesh::instance()->get_field(*fid)->getCoreField();
  if (!isFrozen(f)) freeze(f);
  *pts=getArrayData(f);
}
// add field2 to field1
void m3dc1_field_add(FieldID*  fid1, FieldID* fid2)
{
  m3dc1_field * fld1 = m3dc1_mesh::instance()->get_field(*fid1);
  m3dc1_field * fld2 = m3dc1_mesh::instance()->get_field(*fid2);
  int dof_per_nd_1 = fld1->getDofsPerNode();
  int dof_per_nd_2 = fld2->getDofsPerNode();
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  std::vector<double> dofs1(dof_per_nd_1 * sz);
  std::vector<double> dofs2(dof_per_nd_2 * sz);
  int dofMin = std::min(dof_per_nd_1,dof_per_nd_2);
  int num_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  int dofPerEntDummy[2];
  for (int vtx=0; vtx<num_vtx; ++vtx)
  {
    m3dc1_ent_getdofdata (&vtx_dim, &vtx, fid1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vtx_dim, &vtx, fid2, dofPerEntDummy+1, &dofs2[0]);
    for (int ii = 0; ii < dofMin * sz; ++ii)
      dofs1.at(ii) += dofs2.at(ii);
    m3dc1_ent_setdofdata (&vtx_dim, &vtx, fid1, &dof_per_nd_1, &dofs1[0]);
  }
}
void m3dc1_field_mult(FieldID * fid, double * scale)
{
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int dof_per_nd = fld->getDofsPerNode();
  double * dofs = new double[dof_per_nd]();
  for (int vtx = 0; vtx < num_lcl_vtx; ++vtx)
  {
    m3dc1_ent_getdofdata(&vtx_dim, &vtx, fid, &dof_per_nd, dofs);
#ifndef PETSC_USE_COMPLEX
    for (int ii = 0; ii < dof_per_nd; ++ii)
      dofs[ii] *= *scale;
#else
    for (int ii = 0; ii < dof_per_nd; ++ii)
    {
      double rl = dofs[2 * ii];
      double img = dofs[2 * ii + 1];
      dofs[2 * ii] = rl * scale[0] - img * scale[1];
      dofs[2 * ii + 1] = img * scale[0] + rl * scale[1];
    }
#endif
    m3dc1_ent_setdofdata(&vtx_dim, &vtx, fid, &dof_per_nd, &dofs[0]);
  }
  delete [] dofs;
}
void m3dc1_field_assign(FieldID * fid, double * val)
{
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int dof_per_nd = fld->getDofsPerNode();
  int num_vtx=m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  std::vector<double> dofs(dof_per_nd * sz, val[0]);
#ifdef PETSC_USE_COMPLEX
  for (int ii = 0; ii < dof_per_nd; ++ii)
    dofs[2 * ii + 1] = val[1];
#endif
  for (int vtx = 0; vtx < num_vtx; ++vtx)
    m3dc1_ent_setdofdata(&vtx_dim, &vtx, fid, &dof_per_nd, &dofs[0]);
}
void m3dc1_field_copy(FieldID * fid1, FieldID * fid2)
{
  m3dc1_field * fld1 = m3dc1_mesh::instance()->get_field(*fid1);
  m3dc1_field * fld2 = m3dc1_mesh::instance()->get_field(*fid1);
  int dof_per_nd_1 = fld1->getDofsPerNode();
  int dof_per_nd_2 = fld2->getDofsPerNode();
  assert(dof_per_nd_1 = dof_per_nd_2);
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  std::vector<double> dofs(dof_per_nd_1 * sz);
  int num_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  int num_dof = -1;
  for(int vtx = 0; vtx < num_vtx; ++vtx)
  {
    m3dc1_ent_getdofdata(&vtx_dim, &vtx, fid2, &num_dof, &dofs[0]);
    m3dc1_ent_setdofdata(&vtx_dim, &vtx, fid1, &num_dof, &dofs[0]);
  }
}
void m3dc1_field_retrieve (FieldID*  fid, double* /*out*/ data, int* size)
{
  int scalar_type=0;
#ifdef PETSC_USE_COMPLEX
  scalar_type=1;
#endif
  double* pts=NULL;
  m3dc1_field_getdataptr(fid, &pts);
  memcpy(data, pts, *size*(1+scalar_type)*sizeof(double));
}
void m3dc1_field_set (FieldID*  fid, double* /*in*/ data, int* size)
{
  int scalar_type=0;
#ifdef PETSC_USE_COMPLEX
  scalar_type=1;
#endif
  double* pts=NULL;
  m3dc1_field_getdataptr (fid, &pts);
  memcpy(pts, data, *size*(1+scalar_type)*sizeof(double));
}
void m3dc1_field_insert(FieldID * fid,
                        int * lcl_dofs,
                        int * num_dofs,
                        double * vals,
                        int * op)
{
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  int num_lcl_dof = 0;
  m3dc1_field_getnumlocaldof(fid, &num_lcl_dof);
  double * dat = NULL;
  m3dc1_field_getdataptr(fid, &dat);
  if(*op == 0) // set value
    for(int ii = 0; ii < *num_dofs; ++ii)
      for(int tp = 0; tp < sz; ++tp)
        dat[lcl_dofs[ii + tp]] = vals[ii + tp];
  else
    for(int ii = 0; ii < *num_dofs * sz; ++ii)
      for(int tp = 0; tp < sz; ++tp)
        dat[lcl_dofs[ii + tp]] = vals[ii + tp];
}
void m3dc1_field_load(FieldID * fid, const char * filename)
{
  if (!PCU_Comm_Self())
    std::cout << "[M3D-C1 INFO] " << __func__ << "(field id " << *fid << ", file \"" << filename << "\")\n";
  load_field(m3dc1_mesh::instance()->get_mesh(), *fid, filename);
}
void m3dc1_field_write(FieldID * fid, const char * filename, int * start_index)
{
  if (!PCU_Comm_Self())
    std::cout << "[M3D-C1 INFO] " << __func__ << "(field id " << *fid << ", file \"" << filename << "\")\n";
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  write_field(msh->get_mesh(), fld, filename, *start_index);
}
void m3dc1_field_print(FieldID* fid)
{
  /*
  if (!m->findField("node global id field"))
  {
    if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<" failed as global node id not found\n";
    return;
  }
  */
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  int num_dof = fld->getDofsPerNode();
  double * dofs = new double[num_dof];
  apf::Mesh2 * m = m3dc1_mesh::instance()->get_mesh();
  apf::MeshEntity * e = NULL;
  apf::MeshIterator * it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    get_ent_dofdata(fld, e, &dofs[0]);
    std::cout << "[p" << PCU_Comm_Self() << "] field " << fld->getName()
              << "/ent " << get_ent_globalid(m, e)
              <<": [" << dofs[0];
    for(int ii = 1; ii < num_dof; ++ii)
      std::cout << ", " << dofs[ii];
    std::cout << "]" << std::endl;
  }
  m->end(it);
  delete [] dofs;
}
void m3dc1_field_compare(FieldID * fid1, FieldID * fid2)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld1 = msh->get_field(*fid1);
  m3dc1_field * fld2 = msh->get_field(*fid2);
  apf::Field * f1 = fld1->getCoreField();
  apf::Field * f2 = fld2->getCoreField();
  double * field_data_1 = getArrayData(f1);
  double * field_data_2 = getArrayData(f2);
  int num_dof_1 = fld1->getDofsPerNode();
  int num_dof_2 = fld2->getDofsPerNode();
  assert(num_dof_1 == num_dof_2);
  int num_lcl_vtx = msh->get_mesh()->count(0);
  for(int ii = 0; ii < num_dof_1 * num_lcl_vtx; ++ii)
  {
    if(!m3dc1_double_isequal(field_data_1[ii], field_data_2[ii]))
    {
      std::cout << "[M3D-C1] : "
                << apf::getName(f1) << "[" << ii << "] =" << field_data_1[ii]
                << ", "
                << apf::getName(f2) << "[" << ii << "] =" << field_data_2[ii]
                << std::endl;
      break;
    }
  }
}
void m3dc1_ent_getlocaldofid(int * edim, int * eid, FieldID * fid, int * dof_id, int * dof_cnt)
{
  (void)edim;
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  get_ent_localdofid(fld, *eid, dof_id, dof_cnt);
}
void m3dc1_ent_getglobaldofid(int * edim,
                              int * eid,
                              FieldID * fid,
                              int * dof_id,
                              int * dof_cnt)
{
  (void)edim;
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  get_ent_globaldofid(fld, *eid, dof_id, dof_cnt);
}
void m3dc1_ent_getnumdof(int * edim,
                         int * eid,
                         FieldID * fid,
                         int * num_dof)
{
  (void)edim;
  (void)eid;
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  *num_dof = fld->getDofsPerNode();
}
void m3dc1_ent_setdofdata(int * edim,
                          int * eid,
                          FieldID * fid,
                          int * num_dof,
                          double * dof_data)
{
  apf::MeshEntity * e = apf::getMdsEntity(m3dc1_mesh::instance()->get_mesh(), *edim, *eid);
  assert(*edim==0 && e);
  m3dc1_field * fld = m3dc1_mesh::instance()->get_field(*fid);
  assert(*num_dof == fld->getDofsPerNode());
  set_ent_dofdata(fld, e, dof_data);
}
void m3dc1_ent_getdofdata(int * edim,
                          int * eid,
                          FieldID * fid,
                          int * num_dof,
                          double * dof_data)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  apf::MeshEntity * e = apf::getMdsEntity(msh->get_mesh(), *edim, *eid);
  m3dc1_field * fld = msh->get_field(*fid);
  *num_dof = fld->getDofsPerNode( );
  // todo: check node count on ent type
  get_ent_dofdata(fld, e, dof_data);
}
#ifdef M3DC1_PETSC
// matrix and solver functions
void m3dc1_matrix_create(int * mid, int * matrix_type, int * scalar_type, FieldID * fid)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  m3dc1_field * fld = msh->get_field(*fid);
  MPI_Comm cm = *matrix_type == M3DC1_LOCAL_MAT ? MPI_COMM_SELF :
                *matrix_type == M3DC1_PARALLEL_MAT ? M3DC1_COMM_WORLD : MPI_COMM_NULL;
  m3dc1_matrix * mat = new m3dc1_matrix(*mid,*scalar_type,msh,fld,cm);
  m3dc1_solver::instance()->add_matrix(*mid, (m3dc1_matrix*)mat);
}
void m3dc1_matrix_assemble(int * mid)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist.");
  mat->fix();
}

void m3dc1_matrix_delete(int * mid)
{
  m3dc1_solver::instance()->destroy_matrix(*mid);
}
void m3dc1_matrix_reset(int* mid)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist.");
  mat->zero();
}
void m3dc1_matrix_insert(int * mid, int * row, int * col, int * scalar_type, double * val)
{
  (void) scalar_type;
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist.");
  mat->set_values(1,row,1,col,val);
}
void m3dc1_matrix_add(int * mid, int * row, int * col, int * scalar_type, double * val)
{
  (void) scalar_type;
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist.");
  mat->add_values(1,row,1,col,val);
}
void m3dc1_matrix_setbc(int * mid, int * row)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist");
  double one = 1.0;
  mat->zero_rows(1,row);
  mat->set_values(1,row,1,row,&one);
}
void m3dc1_matrix_setlaplacebc(int * mid, int * rw, int * cl_cnt, int * cols, double * vals)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with speficied ID does not exist");
  mat->zero_rows(1,rw);
  mat->set_values(1,rw,*cl_cnt,cols,vals);
}
void m3dc1_matrix_solve(int * mid, FieldID * rhs)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist");
  if(!mat->is_fixed())
    mat->fix();
  m3dc1_field * rhsf = m3dc1_mesh::instance()->get_field(*rhs);
  mat->solve(rhsf);
}
void m3dc1_matrix_multiply(int * mid, FieldID * in, FieldID * out)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist");
  m3dc1_field * inf = m3dc1_mesh::instance()->get_field(*in);
  m3dc1_field * outf = m3dc1_mesh::instance()->get_field(*out);
  mat->multiply(inf,outf);
}
void m3dc1_matrix_getfieldid(int * mid, int * fid)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  m3dc1_field * fld = mat->get_field();
  *fid = fld->getId();
}
void m3dc1_matrix_getnumiter(int * mid, int * iter_num)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  assert(mat && "[M3D-C1 Error] Matrix with specified ID does not exist");
  *iter_num = mat->solver_iteration_count();
}
// insert the blocks associated with a single node of the given element
// mid is the matrix identifier
// dim is the dimension of the element (0,1,2,3)
// eid is the element id of dimension dim
// nd1 is the first node, in canonical ordering on the element using the field order used to create the matrix
// nd2 is the second node, in the canonical ordering on the element using the field order used to create the matrix
// vals is the logically 2d-array of values for the node's blocks
void m3dc1_matrix_insertnodeblocks(int * mid,
                                   int * dim,
                                   int * eid,
                                   int * nd1,
                                   int * nd2,
                                   double * vals)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  insert_node_blocks(mat,*dim,*eid,*nd1,*nd2,vals);
}
void m3dc1_matrix_insertentblocks(int * mid,
                                  int * edim,
                                  int * eid,
                                  double * vals)
{
  m3dc1_matrix * mat = m3dc1_solver::instance()->get_matrix(*mid);
  insert_element_blocks(mat,*edim,*eid,vals);
}
void m3dc1_matrix_write(int* mid, const char * filename)
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*mid);
  if (!filename)
    return m3dc1_matrix_print(mid);
#ifdef DEBUG
  if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<": matrix id "<<*mid<<", file "<<filename<<"\n";
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*mid<<" does not exist\n";
    return;
  }
#endif
  mat->write(filename);
}
void m3dc1_matrix_print(int* mid)
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*mid);
#ifdef DEBUG
  if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<": matrix id "<<*mid<<"\n";
  if (!mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*mid<<" does not exist\n";
    return;
  }
#endif
  int row = -1;
  int csize = 0;
  size_t sum_csize = 0;
  size_t index = 0;
  vector<int> rows;
  vector<int> n_cols;
  vector<int> cols;
  vector<double> vals;
  mat->get_values(rows, n_cols, cols, vals);
  for(size_t ii = 0; ii < rows.size(); ++ii)
    sum_csize += n_cols[ii];
  assert(vals.size() == sum_csize);
  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": printing matrix "<<*mid<<"\n";
  for (size_t ii = 0; ii < rows.size(); ++ii)
  {
    row = rows[ii];
    csize = n_cols[ii];
    for (int j=0; j<csize; ++j)
    {
      std::cout<<"["<<PCU_Comm_Self()<<"]\t"<<row<<"\t"<<cols[index]<<"\t"<<vals[index]<<"\n";
      ++index;
    }
  }
  assert(index == vals.size());
}
#else
#include "m3dc1_ls.h"
#endif // #ifdef M3DC1_PETSC
int adapt_time = 0;
int adapt_by_field (int * fid, double* psi0, double * psil)
{
  m3dc1_mesh * msh = m3dc1_mesh::instance();
  FILE * fp = fopen("sizefieldParam", "r");
  if(!fp)
  {
    std::cout<<" file sizefieldParam not found "<<std::endl;
    throw 1;
  }
  double param[13];
  set<int> field_keep;
  field_keep.insert(*fid);
  apf::Mesh2* mesh = m3dc1_mesh::instance()->get_mesh();
  for(int i=0; i<13; ++i)
    fscanf(fp, "%lf ", &param[i]);
  fclose(fp);
  apf::Field * psiField = msh->get_field(*fid)->getCoreField();
  // delete all the matrix
#ifdef M3DC1_TRILINOS
  while (m3dc1_ls::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_epetra*>::iterator mat_it = m3dc1_ls::instance()->matrix_container->begin();
    mat_it->second->destroy();
    delete mat_it->second;
    m3dc1_ls::instance()->matrix_container->erase(mat_it);
  }
#endif
#ifdef M3DC1_PETSC
  m3dc1_solver::instance()->destroy_matrices();
#endif
  int valueType = 0;
#ifdef PETSC_USE_COMPLEX
  valueType = 1;
#endif
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);
  double mmax[2];
  double mmin[2];
  m3dc1_model_getmaxcoord(mmax,mmax+1);
  m3dc1_model_getmincoord(mmin,mmin+1);
  ReducedQuinticImplicit shape;
  std::vector<m3dc1_field*> fields;
  std::vector<apf::Field*> apf_fields;
  msh->retrieve_fields(std::back_inserter(fields));
  for(auto fld_it = fields.begin(); fld_it != fields.end(); ++fld_it)
  {
    m3dc1_field * fld = *fld_it;
#ifdef PETSC_USE_COMPLEX
    //group_complex_dof(fld, 1);
#endif
    apf::Field * apf_fld = fld->getCoreField();
    if (isFrozen(apf_fld))
      unfreeze(apf_fld);
    if(!PCU_Comm_Self())
      std::cout << "Solution transfer: add field " << apf::getName(apf_fld) << std::endl;
    apf_fields.push_back(apf_fld);
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    apf::destroyNumbering(n);
  }
  ReducedQuinticTransfer slnTrans(mesh,apf_fields,&shape);
  ma::Input * in = ma::configure(mesh,&sf,&slnTrans);
  //ma::Input* in = ma::configure(mesh,&sfv);
  in->maximumIterations = 9;
  //char filename[256];
  //sprintf(filename,"before%d",adapt_time);
  //apf::writeVtkFiles(filename,mesh);
  //set<int> field_keep;
  //field_keep.insert(*fid);
  //m3dc1_mesh ::instance ()->clean(field_keep);
  in->shouldRunPostZoltan = true;
  ma::adapt(in);
  reorderMdsMesh(mesh);
  m3dc1_mesh::instance()->initialize();
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), 0);
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), m3dc1_mesh::instance()->get_mesh()->getDimension());
  msh->retrieve_fields(std::back_inserter(fields));
  for(auto fld_it = fields.begin(); fld_it != fields.end(); ++fld_it)
  {
    m3dc1_field * fld = *fld_it;
#ifdef PETSC_USE_COMPLEX
    //group_complex_dof(fld, 0);
#endif
    apf::Field * field = fld->getCoreField();
    if(!isFrozen(field))
      freeze(field);
    synchronize_field(field);
  }
  return M3DC1_SUCCESS;
}
double absSize[2]={0,1}, relSize[2]={0.3, 1.5};
int set_mesh_size_bound (double* abs_size, double* rel_size)
{
  for (int i=0; i<2; ++i)
  {
    absSize[i] =  abs_size[i];
    relSize[i] = rel_size[i];
  }
  return M3DC1_SUCCESS;
}
void smooth_size_field (apf::Field* sizeField)
{
  int numVert = m3dc1_mesh::instance()->get_mesh()->count(0);
  for(int ii = 0; ii < numVert; ++ii)
  {
    vector<apf::MeshEntity*> nodes;
    apf::MeshEntity* e = apf:: getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, ii);
    assert(e);
    nodes.push_back(e);
    double sizeOrg = 0.0;
    getComponents(sizeField, e, 0, &sizeOrg);
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->get_mesh()->getAdjacent(e,1,adjacent);
    for(size_t jj = 0; jj < adjacent.getSize(); ++jj)
    {
      apf::Downward downward;
      m3dc1_mesh::instance()->get_mesh()->getDownward(adjacent[jj], 0, downward);
      nodes.push_back(downward[0]==e?downward[1]:downward[0]);
    }
    double size = 0.0;
    for(size_t jj = 0; jj < nodes.size(); ++jj)
    {
      double buff = 0.0;
      getComponents(sizeField, nodes[jj], 0, &buff);
      size += 1.0/buff;
    }
    size /= nodes.size();
    size = 1.0 / size;
    setComponents(sizeField, e, 0, &size);
  }
}
/*
void group_complex_dof (m3dc1_field * fld, int option)
{
  int num_dof = fld->getDofsPerNode();
  int num_dof_double = num_dof * 2;
  std::vector<double> dofs(num_dof_double);
  std::vector<double> newdofs(num_dof_double);
  apf::Mesh * msh = m3dc1_mesh::instance()->get_mesh();
  int num_lcl_vtx = msh->count(0);
  for (int ii = 0; ii < num_lcl_vtx; ++ii)
  {
    apf::MeshEntity * e = apf::getMdsEntity(msh, 0, ii);
    get_ent_dofdata(field, e, &(dofs[0]));
    for (int jj = 0; jj < num_dof/6; ++jj)
    {
      if (option)
      {
        for (int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+k)=dofs.at(2*j*6+2*k);
          newdofs.at(2*j*6+6+k)=dofs.at(2*j*6+2*k+1);
        }
      }
      else
      {
        for (int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+2*k)=dofs.at(2*j*6+k);
          newdofs.at(2*j*6+2*k+1)=dofs.at(2*j*6+6+k);
        }
      }
    }
    // FIXME based on new m3dc1_field
    set_ent_dofdata(field, e, &(newdofs[0]));
  }
}
*/
double p=4;
int set_adapt_p (double* pp)
{
  p = *pp;
  return M3DC1_SUCCESS;
}
int adapt_by_error_field (double * errorData, double * errorAimed, int * max_adapt_node, int * option)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->get_mesh();
  apf::Field* sizeField = createPackedField(m3dc1_mesh::instance()->get_mesh(), "size_field", 1);
  SizeFieldError sf (m3dc1_mesh::instance()->get_mesh(), sizeField); //, *errorAimed);
  int entDim=0, num_lcl_vtx=0;
  m3dc1_mesh_getnument (&entDim, &num_lcl_vtx);
  // first sum error ^ (2d/(2p+d))
  double d=2;
  double errorSum=0;
  if(*option)
  {
    for(int i=0; i<num_lcl_vtx; ++i)
    {
      if(is_ent_original(getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, i)))
        errorSum+=pow(errorData[i],d/(p+d/2.0));
    }
    double errorSumBuff=errorSum;
    MPI_Allreduce(&errorSumBuff, &errorSum, 1, MPI_DOUBLE, MPI_SUM, M3DC1_COMM_WORLD);
    errorSum = *errorAimed*(*errorAimed)/errorSum;
    errorSum = pow(errorSum,1./(2.*p));
  }
  else errorSum=pow(*errorAimed,1./(p+d/2.));
  double size_estimate=0;
  for(int i=0; i<num_lcl_vtx; ++i)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, i);
    assert(e);
    if(!is_ent_original(e)) continue;
    double targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    size_estimate+=max(1.,1./targetSize/targetSize);
  }
  double size_estimate_buff=size_estimate;
  MPI_Allreduce(&size_estimate_buff, &size_estimate, 1, MPI_DOUBLE, MPI_SUM, M3DC1_COMM_WORLD);
  int numNodeGlobl=0, dim=0;
  m3dc1_mesh_getnumglobalent(&dim, &numNodeGlobl);
  cout<<" num_lcl_vtx "<<numNodeGlobl<<" size_estimate "<<size_estimate;
  if(size_estimate>*max_adapt_node) errorSum*=sqrt(size_estimate/(*max_adapt_node));
  for(int i=0; i<num_lcl_vtx; ++i)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->get_mesh(), 0, i);
    assert(e);
    double size = sf.getSize(e);
    assert(errorData[i]==errorData[i]);
    double targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    if(targetSize>relSize[1]) targetSize=relSize[1]; // not too much coarsening
    if(targetSize<relSize[0]) targetSize=relSize[0]; // not too much refining
    targetSize*=size;
    if(targetSize>absSize[1]) targetSize=absSize[1];
    if(targetSize<absSize[0]) targetSize=absSize[0];
    setComponents(sizeField, e, 0, &targetSize);
  }
  smooth_size_field(sizeField);
  smooth_size_field(sizeField);
  synchronize_field(sizeField);
  m3dc1_solver::instance()->destroy_matrices();
  ReducedQuinticImplicit shape;
  std::vector<apf::Field*> fields;
  fields.push_back(sizeField);
  std::vector<m3dc1_field*> m3_flds;
  m3dc1_mesh::instance()->retrieve_fields(std::back_inserter(m3_flds));
  for(auto fld_it = m3_flds.begin(); fld_it != m3_flds.end(); ++fld_it)
  {
    m3dc1_field * fld = *fld_it;
#ifdef PETSC_USE_COMPLEX
    //group_complex_dof(fld, 1);
#endif
    apf::Field * field = fld->getCoreField();
    if (isFrozen(field))
      unfreeze(field);
    fields.push_back(field);
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
    apf::destroyNumbering(n);
  }
  char filename[256];
  sprintf(filename,"before%d",adapt_time);
  //apf::writeVtkFiles(filename,mesh);
  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
  in->maximumIterations = 5;
  in->shouldRunPostZoltan = true;
  ma::adapt(in);
  reorderMdsMesh(mesh);
  m3dc1_mesh::instance()->initialize();
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), 0);
  compute_globalid(m3dc1_mesh::instance()->get_mesh(), m3dc1_mesh::instance()->get_mesh()->getDimension());
  for(auto fld_it = m3_flds.begin(); fld_it != m3_flds.end(); ++fld_it)
  {
    m3dc1_field * fld = *fld_it;
#ifdef PETSC_USE_COMPLEX
//    group_complex_dof(fld, 0);
#endif
    apf::Field* field = fld->getCoreField();
    if(!isFrozen(field))
      freeze(field);
    synchronize_field(field);
  }
  destroyField(sizeField);
  return M3DC1_SUCCESS;
}
int sum_edge_data (double* data, int* size)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->get_mesh();
  int num_edge=m3dc1_mesh::instance()->get_mesh()->count(1), edg_dim=1;
  PCU_Comm_Begin();
  for (int i=0; i<num_edge; ++i)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    int own_partid=get_ent_ownpartid(e);
    apf::MeshEntity* own_e = get_ent_owncopy(e);
    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(data[i*(*size)]),(*size)*sizeof(double));
  }
  PCU_Comm_Send();
  double* receive_buff = new double [*size];
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; ++i)
        data[iedge*(*size)+i]+=receive_buff[i];
    }
  PCU_Comm_Begin();
  for (int i=0; i<num_edge; ++i)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    if (!is_ent_original(e) || !m->isShared(e))
      continue;
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; ++i)
        data[iedge*(*size)+i]=receive_buff[i];
    }
   delete [] receive_buff;
   return M3DC1_SUCCESS;
}
int get_node_error_from_elm (double* elm_data, int* size, double* nod_data)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->get_mesh();
  int num_node=m3dc1_mesh::instance()->get_mesh()->count(0);
  int num_elm=m3dc1_mesh::instance()->get_mesh()->count(2);
  int nod_dim=0;
  PCU_Comm_Begin();
  double * buff = new double[*size];
  double * area = new double[num_node];
  memset(&area[0],0.0,sizeof(double)*num_node);
  for (int ii = 0; ii < (*size)*num_elm; ++ii)
    for (int jj = 0; jj < *size; ++jj)
    {
      assert(elm_data[(*size)*ii+jj]==elm_data[(*size)*ii+jj]);
      assert(elm_data[(*size)*ii+jj]>=0);
    }
  for (int ii=0; ii<num_node; ++ii)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, ii);
    int own_partid=get_ent_ownpartid(e);
    apf::MeshEntity* own_e = get_ent_owncopy(e);
    apf::Adjacent adjacent;
    m->getAdjacent(e,2,adjacent);
    for (size_t jj = 0; jj < adjacent.getSize(); ++jj)
    {
       apf::MeshElement * me = createMeshElement(m, adjacent[jj]);
       double s = apf::measure(me);
       int ielm = getMdsIndex(m, adjacent[jj]);
       assert(ielm>=0 &&ielm<num_elm);
       for (int kk = 0; kk < *size; kk++)
          nod_data[ii*(*size)+kk] += s*elm_data[(*size)*ielm+kk];
       area[ii] += s;
       destroyMeshElement(me);
    }
    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_COMM_PACK(own_partid, area[ii]);
    PCU_Comm_Pack(own_partid, &(nod_data[(*size)*ii]), sizeof(double)*(*size));
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      int vtx = getMdsIndex(m, node);
      double s;
      PCU_COMM_UNPACK(s);
      area[vtx]+=s;
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      for (int i = 0; i < *size; ++i)
        nod_data[vtx*(*size)+i]+=buff[i];
    }
  for (int i=0; i<num_node; ++i)
  {
    for (int j=0; j<*size; ++j)
      nod_data[i*(*size)+j]/=area[i];
  }
  PCU_Comm_Begin();
  for (int i=0; i<num_node; ++i)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    if (!is_ent_original(e) || !m->isShared(e))
      continue;
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(nod_data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      int vtx = getMdsIndex(m, node);
      for (int i = 0; i < *size; ++i)
        nod_data[vtx*(*size)+i]=buff[i];
    }
  delete []buff;
  delete []area;
  return M3DC1_SUCCESS;
}
void m3dc1_field_max (FieldID * fid, double * max_val, double * min_val)
{
  m3dc1_field* fld = m3dc1_mesh::instance()->get_field(*fid);
  int dof_per_nd = fld->getDofsPerNode();
  int sz = 1;
#ifdef PETSC_USE_COMPLEX
  sz = 2;
#endif
  int num_lcl_vtx = m3dc1_mesh::instance()->get_mesh()->count(0);
  int vtx_dim = 0;
  std::vector<double> max_cmp(dof_per_nd * sz, -1e30);
  std::vector<double> min_cmp(dof_per_nd * sz,  1e30);
  std::vector<double> dofs(dof_per_nd);
  for (int vtx = 0; vtx < num_lcl_vtx; ++vtx)
  {
    m3dc1_ent_getdofdata(&vtx_dim, &vtx, fid, &dof_per_nd, &dofs[0]);
    for (int ii = 0; ii < dof_per_nd; ++ii)
    {
      if(max_cmp[ii] < dofs[ii])
        max_cmp[ii] = dofs[ii];
      if(min_cmp[ii] > dofs[ii])
        min_cmp[ii] = dofs[ii];
    }
  }
  MPI_Allreduce(&(max_cmp[0]),
                max_val,
                dof_per_nd,
                MPI_DOUBLE,
                MPI_MAX,
                M3DC1_COMM_WORLD);
  MPI_Allreduce(&(min_cmp[0]),
                min_val,
                dof_per_nd,
                MPI_DOUBLE,
                MPI_MIN,
                M3DC1_COMM_WORLD);
}

