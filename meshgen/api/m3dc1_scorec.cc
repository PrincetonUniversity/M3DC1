/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include "m3dc1_mesh.h"
#include "m3dc1_ghost.h"
#include "m3dc1_field.h"
#include <mpi.h>
#include <PCU.h>
#include <gmi_analytic.h>
#include <map>
#include <cstring>
#include "apfMDS.h"
#include "apfNumbering.h"
#include "Expression.h"
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "apfMesh2.h"
#include "apfShape.h"

bool m3dc1_double_isequal(double A, double B)
{
  double maxDiff = 1e-8;
  double maxRelDiff = 1e-8;
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

bool is_ent_ghost(apf::Mesh2* m, apf::MeshEntity* e, int ent_dim);

//*******************************************************
void destroy_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  if (m->findField("node own partid field"))
    destroyField(m->findField("node own partid field"));
  if (m->findField("node global id field"))
    destroyField(m->findField("node global id field"));
}

//*******************************************************
void generate_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
//  if (!PCU_Comm_Self())  std::cout<<"[M3D-C1 INFO] ***** GENERATING GLOBAL NODE NUMBERING ***** \n"; 
  destroy_node_global_numbering(m);
  double id[1];
  double own_partid[1];
  apf::Field* node_ownpid_f = createPackedField(m, "node own partid field", 1);
  apf::freeze(node_ownpid_f);
  apf::Field* node_globalid_f = createPackedField(m, "node global id field", 1);
  apf::freeze(node_globalid_f);
  // generate global node_id
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0];
  apf::MeshEntity* e;
  PCU_Exscan_Ints(&num_own_ent,1);
  int start=num_own_ent;
  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid[0]=(double)get_ent_ownpartid(m,e); 
    setComponents(node_ownpid_f, e, 0, own_partid);    
    if ((int)(own_partid[0])!=PCU_Comm_Self()) continue;
    id[0] = (double) start;
    setComponents(node_globalid_f, e, 0, id);
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&start,sizeof(int));
    }
    ++start;
  }
  m->end(it);
  PCU_Comm_Send();

  int value;
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&value,sizeof(int));
      id[0] = (double) value;
      setComponents(node_globalid_f, r, 0, id);
    }
}


//*******************************************************
void destroy_elem_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  if (m->findField("elem global id field"))
    destroyField(m->findField("elem global id field"));
}

//*******************************************************
void generate_elem_global_numbering(apf::Mesh2* m)
//*******************************************************
{
//  if (!PCU_Comm_Self())  std::cout<<"[M3D-C1 INFO] ***** GENERATING GLOBAL ELEMENT NUMBERING ***** \n"; 
  destroy_elem_global_numbering(m);

  double id[1];
  apf::Field* elem_globalid_f = createPackedField(m, "elem global id field", 1,apf::getConstant(m->getDimension()));
  apf::freeze(elem_globalid_f);

  // generate global elem_id
  int own_partid, mesh_dim = m->getDimension();
  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[mesh_dim];
  apf::MeshEntity* e;
  PCU_Exscan_Ints(&num_own_ent,1);
  int start=num_own_ent;

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(mesh_dim);
  while ((e = m->iterate(it)))
  {
    if (get_ent_ownpartid(m,e)!=PCU_Comm_Self()) continue;
    id[0] = (double) start;
    setComponents(elem_globalid_f, e, 0, id);
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&start,sizeof(int));
    }
    ++start;
  }
  m->end(it);
  PCU_Comm_Send();

  int value;
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&value,sizeof(int));
      id[0] = (double) value;
      setComponents(elem_globalid_f, r, 0, id);
    }
}

void destroy_global_numbering(apf::Mesh2* m)
{
  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[NODE_GLB_ORDER]->get_field();
  std::string name = getName(f);
  destroyNumbering(m->findNumbering(name.c_str()));
  destroyField(f);
  delete (*m3dc1_mesh::instance()->field_container)[NODE_GLB_ORDER];
  m3dc1_mesh::instance()->field_container->erase(NODE_GLB_ORDER);

  destroy_node_global_numbering(m);
  destroy_elem_global_numbering(m);
}

void generate_global_numbering(apf::Mesh2* m)
{
  get_global_numbering();
  generate_node_global_numbering(m);
  generate_elem_global_numbering(m);
}

#include <sys/resource.h>
double start_time;
double _gettime()
{
  struct rusage ruse_now;
  getrusage(RUSAGE_SELF, &ruse_now);
  return double(ruse_now.ru_utime.tv_sec) + double(ruse_now.ru_utime.tv_usec)/1000000.0;
}

//*******************************************************
int m3dc1_scorec_init()
//*******************************************************
{ 
  PCU_Comm_Init();
  start_time=_gettime();
  m3dc1_model :: instance()-> model = gmi_make_analytic();
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_scorec_finalize()
//*******************************************************
{ 
  // destroy global numbering and associated field
  // delete existing numbering
  apf::Numbering* old_n = m3dc1_mesh::instance()->mesh->findNumbering("m3dc1_global_node_id");
  if (old_n) destroyNumbering(old_n);
  int node_glb_order=NODE_GLB_ORDER; 
  m3dc1_field_delete (&node_glb_order);

  //actually destroy badly designed singletons.
  //this was chosen over statically allocating the objects
  //and having the runtime deallocate them in order to avoid
  //possible issues linking to FORTRAN.
  //feel free to make them static objects and see if that works
  m3dc1_mesh::destroy();
  m3dc1_model::destroy();
  // print elapsed time and increased heap memory
  if (!PCU_Comm_Self()) std::cout<<"[PUMI INFO] elapsed time "<<_gettime()-start_time<<"(sec)\n";

  PCU_Comm_Free();
  return M3DC1_SUCCESS; 
}
/** plane functions */

//*******************************************************
int m3dc1_plane_setnum(int* num_plane)
//*******************************************************
{
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) return M3DC1_FAILURE;
  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_plane_getnum(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_plane_getid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphirange(double* min_val, double* max_val)
//*******************************************************
{
  m3dc1_model::instance()->set_phi(*min_val, *max_val);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphi(int* planeid, double* phi)
//*******************************************************
{
  m3dc1_model::instance()->set_phi(*planeid, *phi);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_getphi(int* planeid, double* phi)
//*******************************************************
{
  m3dc1_model::instance()->get_phi(*planeid, phi);
  return M3DC1_SUCCESS; 
}

/** model functions */
//*******************************************************
int m3dc1_model_getplaneid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
}

//*******************************************************
int m3dc1_model_getedge (int*  /* out */  left_edge, 
                         int*  /* out */  right_edge,
                         int*  /* out */  bottom_edge, 
                         int*  /* out */  top_edge)
//*******************************************************
{
  std::cout<<"m3dc1_model_getedge not implemented"; throw 1;
}

//*******************************************************
int m3dc1_model_getmincoord(double* x_min, double* y_min)
//*******************************************************
{
  *x_min = m3dc1_model::instance()->boundingBox[0];
  *y_min = m3dc1_model::instance()->boundingBox[1];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getmaxcoord(double* x_max, double* y_max)
//*******************************************************
{
  double mincoord[3],maxcoord[3];
  *x_max = m3dc1_model::instance()->boundingBox[2];
  *y_max = m3dc1_model::instance()->boundingBox[3];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_load(const char* /* in */ model_file)
//*******************************************************
{  
  std::string str_model_name(model_file);
  m3dc1_model::instance()->load_analytic_model(model_file); 
  m3dc1_model::instance()->caculateBoundingBox();
  // save the num of geo ent on the oringal plane
  m3dc1_model::instance()->numEntOrig[0]=m3dc1_model::instance()->model->n[0];
  m3dc1_model::instance()->numEntOrig[1]=m3dc1_model::instance()->model->n[1];
  m3dc1_model::instance()->numEntOrig[2]=m3dc1_model::instance()->model->n[2];

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_print()
//*******************************************************
{
  if (PCU_Comm_Self() || m3dc1_model::instance()->local_planeid) 
    return M3DC1_SUCCESS;

  double min[3], max[3];
  gmi_iter* gf_it = gmi_begin(m3dc1_model::instance()->model, 2);
  gmi_ent* ge;
  int count = 0;
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

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_setnumplane(int* num_plane)
//*******************************************************
{
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) return M3DC1_FAILURE;
  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getnumplane(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_setpbc (int* /* in */ x_pbc, 
                        int* /* in */ y_pbc)
//*******************************************************
{
  m3dc1_model::instance()->xperiodic = *x_pbc;
  m3dc1_model::instance()->yperiodic = *y_pbc;
  return  M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getpbc (int* /* out */ x_pbc, 
                        int* /* out */ y_pbc)
//*******************************************************
{
  *x_pbc = m3dc1_model::instance()->xperiodic;
  *y_pbc = m3dc1_model::instance()->yperiodic;
  return  M3DC1_SUCCESS;
}

/** mesh functions */
#include <parma.h>

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    m->setDoubleTag(e, tag, &w);
  m->end(it);
}

apf::MeshTag* setWeights(apf::Mesh* m) {
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}

void clearTags(apf::Mesh* m, apf::MeshTag* t) {
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}

#include <unistd.h>
//*******************************************************
int m3dc1_mesh_load(const char* mesh_file)
//*******************************************************
{
  if (m3dc1_model::instance()->local_planeid == 0) 
  {
    m3dc1_mesh::instance()->mesh = apf::loadMdsMesh(m3dc1_model::instance()->model, mesh_file);
    apf::disownMdsModel(m3dc1_mesh::instance()->mesh);
    /* vertex load balancing */
//    Parma_PrintPtnStats(m3dc1_mesh::instance()->mesh, "initial");
    // will not work not non-man geo 
    //m3dc1_mesh::instance()->mesh->verify();
    apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
    // remove field and numbering loaded from the file
    while(mesh->countFields())
      destroyField(mesh->getField(0));
    while(mesh->countNumberings())
      destroyNumbering(mesh->getNumbering(0));

    // remove tags loaded from the file except for "normal curvature"
    apf::DynamicArray<apf::MeshTag*> tags;
    mesh->getTags(tags);
    for(int i=0; i<tags.getSize(); i++)
    {
      if (mesh->findTag("norm_curv")==tags[i]) continue;
      for(int idim=0; idim<4; idim++)
        apf::removeTagFromDimension(mesh, tags[i], idim);
      mesh->destroyTag(tags[i]);
    }
  }
  else {
    m3dc1_mesh::instance()->mesh = apf::makeEmptyMdsMesh(m3dc1_model::instance()->model, 2, false);
  }
  m3dc1_mesh::instance()->initialize();
  // Set up global entity numbering for nodes and elements
  if (m3dc1_model::instance()->num_plane==1) // 2D
    generate_global_numbering(m3dc1_mesh::instance()->mesh);

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_mesh_build3d (int* num_field, int* field_id,  
                        int* num_dofs_per_value)
//*******************************************************
{ 
  // switch COMM to GLOBAL COMM
  MPI_Comm groupComm = PCU_Get_Comm();
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);

  // initialize phi value and construct 3d
  int num_plane = m3dc1_model::instance()->num_plane;
  m3dc1_model::instance()->set_phi(0.0, 2.0*M3DC1_PI/num_plane*(num_plane-1));
  m3dc1_model::instance()->setupCommGroupsPlane();
  m3dc1_mesh::instance()->build3d(*num_field, field_id, num_dofs_per_value);
  generate_global_numbering(m3dc1_mesh::instance()->mesh);
  // now construct 3d mesh
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  if (!m3dc1_ghost::instance()->num_local_ent[0])   
    *num_ent = m3dc1_mesh::instance()->num_local_ent[*ent_dim];
  else
    *num_ent = m3dc1_ghost::instance()->num_local_ent[*ent_dim];
    
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;

  *num_ent = m3dc1_mesh::instance()->num_own_ent[*ent_dim];
   
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_mesh::instance()->num_global_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_setordering (int* option)
//*******************************************************
{
  if (*option < M3DC1_NO_ORDER || *option > M3DC1_ADJ_SOLVER_ORDER)
    return M3DC1_FAILURE;

  m3dc1_mesh::instance()->ordering_opt = *option;
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getordering (int* option)
//*******************************************************
{
  *option = m3dc1_mesh::instance()->ordering_opt;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_mesh_print()
//*******************************************************
{
  apf::Mesh2* m = NULL; 
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else 
    m = m3dc1_ghost::instance()->mesh;
  apf::MeshEntity* e;

  // print vertex
  apf::Vector3 xyz;
  apf::Copies remotes;
  for (int pid=0; pid<PCU_Comm_Peers();++pid)
  {
    if (pid==PCU_Comm_Self())
    {
      apf:: MeshIterator* it = m->begin(0);
      while ((e = m->iterate(it)))
      {
        m->getPoint(e, 0, xyz);
        if (is_ent_ghost(m, e, 0)) 
          std::cout<<"(p"<<PCU_Comm_Self()<<") ghost vtx "<<get_ent_globalid(m,e,0)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", coord ["<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]
                   <<"]\n";

       else 
         std::cout<<"(p"<<PCU_Comm_Self()<<") vtx "<<get_ent_globalid(m,e,0)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", coord ["<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]
                   <<"]\n";
        if (m->isShared(e))
        {
          m->getRemotes(e,remotes);
          APF_ITERATE(apf::Copies,remotes,it)
          std::cout<<"(p"<<PCU_Comm_Self()<<") vtx "<<get_ent_globalid(m,e,0)
                   <<" remote ("<<it->first<<","<< it->second<<")\n";
        }        
      } // while
      m->end(it);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  } // for

  int num_vtx, mesh_dim = m->getDimension();
  apf::Downward downward;

  for (int pid=0; pid<PCU_Comm_Peers();++pid)
  {
    if (pid==PCU_Comm_Self())
    {
      apf:: MeshIterator* it = m->begin(mesh_dim);
      while ((e = m->iterate(it)))
      {
        num_vtx = m->getDownward(e, 0, downward);
        int* vtx_id=new int[num_vtx];
        for (int i=0; i<num_vtx; ++i)
          vtx_id[i] = get_ent_globalid(m,downward[i],0);
        if (is_ent_ghost(m, e, mesh_dim))
        {
          if (num_vtx==3)
            std::cout<<"(p"<<PCU_Comm_Self()<<") ghost elem "<<get_ent_globalid(m,e,mesh_dim)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", vtx ("<<vtx_id[0]<<", "<<vtx_id[1]
                   <<", "<<vtx_id[2]<<")\n";
          else if (num_vtx==4)
            std::cout<<"(p"<<PCU_Comm_Self()<<") ghost elem "<<get_ent_globalid(m,e,mesh_dim)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", vtx ("<<vtx_id[0]<<", "<<vtx_id[1]
                   <<", "<<vtx_id[2]<<", "<<vtx_id[3]<<")\n";
        }
        else
        {
          if (num_vtx==3)
            std::cout<<"(p"<<PCU_Comm_Self()<<") elem "<<get_ent_globalid(m,e,mesh_dim)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", vtx ("<<vtx_id[0]<<", "<<vtx_id[1]
                   <<", "<<vtx_id[2]<<")\n";
          else if (num_vtx==4)
            std::cout<<"(p"<<PCU_Comm_Self()<<") elem "<<get_ent_globalid(m,e,mesh_dim)
                   <<": owner "<<get_ent_ownpartid(m, e)<<", vtx ("<<vtx_id[0]<<", "<<vtx_id[1]
                   <<", "<<vtx_id[2]<<", "<<vtx_id[3]<<")\n";
        }
        if (m->isShared(e))
        {
          m->getRemotes(e,remotes);
          APF_ITERATE(apf::Copies,remotes,it)
          std::cout<<"(p"<<PCU_Comm_Self()<<") elem "<<get_ent_globalid(m,e,mesh_dim)
                   <<" remote ("<<it->first<<","<< it->second<<")\n";
        }  
        delete [] vtx_id;      
      } // while
      m->end(it);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  } // for
}


//*******************************************************
int m3dc1_mesh_write(const char* filename, int *option)
//*******************************************************
{
  apf::Mesh2* m = NULL; 
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else 
    m = m3dc1_ghost::instance()->mesh;
  apf::MeshEntity* e;

  if(*option==0 ||*option==3)
  {  
    int dim=2, num_ent, ent_id, geom_class_dim,geom_class_id;
    m3dc1_mesh_getnument (&dim, &num_ent);
    vector<double> geoId (num_ent);
    apf:: MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      ent_id = getMdsIndex(m, e);
      m3dc1_ent_getgeomclass (&dim, &ent_id, &geom_class_dim, &geom_class_id);
      geoId.at(ent_id)=geom_class_id;
    }
    m->end(it);
    apf::writeVtkFiles(filename,m);     
    int one=1;
    if(*option==3) output_face_data (&one, &geoId[0], "geoId");
  }
  else
  {
    char filename_buff[256];
    sprintf(filename_buff, "%s.smb",filename);
    int fieldID=12;
    double dofBuff[1024];
    m3dc1_field* mf = NULL;
    if (!m3dc1_ghost::instance()->num_local_ent[0])  
      mf = (*(m3dc1_mesh::instance()->field_container))[fieldID];
    else 
      mf = (*(m3dc1_ghost::instance()->field_container))[fieldID];
    apf::Field* f = mf ->get_field();
    int numDof = countComponents(f);
    apf::MeshTag* tag = m->createDoubleTag("field12", numDof);

    const int dim=0;
    apf:: MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      apf::getComponents(f, e, 0, dofBuff);
      m->setDoubleTag(e,tag, dofBuff);
    }
    m->end(it);
    m->writeNative(filename_buff);
    apf::removeTagFromDimension(m, tag, dim);
    m->destroyTag(tag);
  }
}

/* mesh entity functions */
//*******************************************************
int m3dc1_node_getglobalid (int* /* in */ ent_id, int* /* out */ global_ent_id)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;

  apf::MeshEntity* e = apf::getMdsEntity(m, 0, *ent_id);
  *global_ent_id=get_ent_globalid(m, e, 0);
  
#ifdef DEBUG  // debug
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    assert(*global_ent_id==getNumber(get_global_numbering(), e, 0, 0)); 
#endif
  return M3DC1_SUCCESS;
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
    fprintf(np, "%d\t%lf\t%lf\t%lf\n", get_ent_globalid(m,e,0)+start_index, xyz[0],xyz[1],xyz[2]);
  } // while
  m->end(it);
  fclose(np);
}

//*******************************************************
int m3dc1_node_write (const char* filename, int* start_index)
//*******************************************************
{
  apf::Mesh2* m = NULL;
 
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else 
    m = m3dc1_ghost::instance()->mesh;
  write_node(m, filename, *start_index);
  return M3DC1_SUCCESS;
}


// Retrieve global entity id

//*******************************************************
int m3dc1_ent_getglobalid (int* /* in */ ent_dim,
                           int* /* in */ ent_id,
			   int* /* out */ global_ent_id)
//*******************************************************
{
  if (*ent_dim>0 && *ent_dim<m3dc1_mesh::instance()->mesh->getDimension())
  {
    cout<<"[M3D-C1 ERROR] "<<__func__<<": not supported for entity dimension "<<*ent_dim<<"\n";
    return M3DC1_FAILURE;
  }

  apf::Mesh2* m=NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;
    
  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  assert(e);
  *global_ent_id=get_ent_globalid(m, e, *ent_dim);
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id, 
            int* /* out */ geom_class_dim, int* /* out */ geom_class_id)
//*******************************************************
{
  apf::Mesh2* m=NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;
    
   apf::MeshEntity* ent = getMdsEntity(m, *ent_dim, *ent_id);
   gmi_ent* gent = (gmi_ent*)(m->toModel(ent));

  *geom_class_dim = gmi_dim(m3dc1_model::instance()->model,gent);
  *geom_class_id = gmi_tag(m3dc1_model::instance()->model,gent);
  // if 3D mesh, need to return the classification on the original plane
  if (m->getDimension() ==3 )
  {
    int numEntOrig[3];
    int numPlane = m3dc1_model::instance()->num_plane;
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

//*******************************************************
int m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                      int* /* in */ adj_dim, int* /* out */ adj_ent, 
                      int* /* in */ adj_ent_allocated_size, int* /* out */ num_adj_ent)
//*******************************************************
{
  apf::Mesh2* m=NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;
    
   apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);

  if (!e || *adj_dim==*ent_dim)
    return M3DC1_FAILURE;

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    m->getAdjacent(e,*adj_dim,adjacent);
      
    *num_adj_ent = adjacent.getSize();
    if (*adj_ent_allocated_size<*num_adj_ent)
      return M3DC1_FAILURE;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m, adjacent[i]);
  }
  else if (*adj_dim<*ent_dim) 
  {
    apf::Downward downward;
    *num_adj_ent = m->getDownward(e, *adj_dim, downward);
    if (*adj_ent_allocated_size<*num_adj_ent)
      return M3DC1_FAILURE;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m, downward[i]);
      
    //adjust the order to work with m3dc1
    if (m->getDimension()==3 && *ent_dim==3 &&*adj_dim==0 &&adj_ent[0]>adj_ent[3])
    {
      int buff[3];
      memcpy(buff, adj_ent, 3*sizeof(int));
      memcpy(adj_ent, adj_ent+3, 3*sizeof(int));
      memcpy(adj_ent+3, buff, 3*sizeof(int));	
    }
  }
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                         int* /* in */ adj_dim, int* /* out */ num_adj_ent)
//*******************************************************
{
  apf::Mesh2* m=NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;
    
   apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);

  if (!e || *adj_dim==*ent_dim)
    return M3DC1_FAILURE;

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    m->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
  }
  else if (*adj_dim<*ent_dim) 
  {
    apf::Downward downward;
    *num_adj_ent = m->getDownward(e, *adj_dim, downward);
  }
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ owning_partid)
//*******************************************************
{  
  if (*ent_dim>0 && *ent_dim<m3dc1_mesh::instance()->mesh->getDimension())
  {
    std::cout<<"[M3D-C1 ERROR] "<<__func__<<" not supported for entity type "<<*ent_dim<<"\n";
    return M3DC1_FAILURE;
  }

  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;    

  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  assert(e);
  *owning_partid = get_ent_ownpartid(m, e);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_isowned (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ ismine)
//*******************************************************
{
  *ent_id -= 1; //index change from Fortran to C

  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;    

  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;

  if (is_ent_owned(m,e)) 
     *ismine = 1;   // 
  else
     *ismine = 0; 
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_isghost (int* /* in */ ent_dim,
		       int* /* in */ ent_id, 
		       int* /* out */ isghost)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    *isghost = 0;
    return M3DC1_SUCCESS;
  }

  apf::Mesh2* m = m3dc1_ghost::instance()->mesh;    
  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  assert(e);
  assert (*ent_dim==0 || *ent_dim==m->getDimension());

  if (is_ent_ghost(m, e, *ent_dim))
    *isghost=1;
  else 
    *isghost=0;
  return M3DC1_SUCCESS;
}


// node-specific functions
//*******************************************************
int m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;    

  apf::MeshEntity* e = getMdsEntity(m, 0, *node_id);
  if (!e)
    return M3DC1_FAILURE;

  apf::Vector3 xyz;
  m->getPoint(e, 0, xyz);
  for (int i=0; i<3; ++i)
    coord[i] = xyz[i]; 
  return M3DC1_SUCCESS;
}

//*******************************************************
void get_gv_edges(gmi_ent* gvertex, std::vector<gmi_ent*>& gedges)
//*******************************************************
{
  gedges.clear();
  gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, gvertex, 1);
  for (int i=0; i<gv_edges->n; ++i)
    gedges.push_back(gv_edges->e[i]);
  gmi_free_set(gv_edges);
  assert(gedges.size()>=1);
}

//*******************************************************
int m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyzt)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh; 

  apf::MeshEntity* vt = getMdsEntity(m, 0, *node_id);
  assert(vt);
  xyzt[2]=0.0;
  //cout<<"nodnormalvec_ "<<*iNode<<" "<<vt<<endl;
  gmi_ent* gent= (gmi_ent*)(m->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if(gType !=  1 && gType !=  0)
  {
    xyzt[0] = xyzt[1] = 0.0;
    return M3DC1_SUCCESS;
  }

  apf::MeshTag* norm_curv_tag = m->findTag("norm_curv");
  if (norm_curv_tag &&m->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    xyzt[0] = norm_curv[0]; 
    xyzt[1] = norm_curv[1];
    return M3DC1_SUCCESS;
  }
  else
  { // if norm/curv is not attached, evaluate
    apf::Vector3 param(0,0,0);
    m->getParam(vt,param);
    // geo node avage on the connected edges
    if (gType == 0) // node is on the 
    {
      apf::Vector3 vcd_t;
      double vcd[3];
      m->getPoint(vt, 0, vcd_t);
      for(int i=0; i<3; i++) 
        vcd[i]=vcd_t[i];
      std::vector<gmi_ent*> gEdges;
      get_gv_edges(gent, gEdges);
      int numEdgePlane=0;
      double normalvec[3]={0.,0.,0.};
      xyzt[0]=xyzt[1]=xyzt[2]=0;
      if (gEdges.size()<2) 
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: #adjEdge of gVertex="<<gEdges.size()<<" (it should be minimum 2) \n";
      assert(gEdges.size()>=2);
      for(int i=0;i<gEdges.size();i++)
      {
        gmi_ent* pe = gEdges.at(i);
        double cd[3]={0,0,0};
        double paraRange[2];
        gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
        M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
        if(!pn) continue;
        numEdgePlane++;
        M3DC1::evalCoord(paraRange[0],cd, pn);
        if(checkSamePoint2D(vcd,cd))
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
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: numEdgePlane="<<numEdgePlane<<" (it should be 2) \n";
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
      m->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalNormalVector(pn[0],pn[1], param[0], xyzt);
    }
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh; 

  apf::MeshEntity* vt = getMdsEntity(m, 0, *node_id);
  assert(vt);

  apf::MeshTag* norm_curv_tag = m->findTag("norm_curv");
  if (norm_curv_tag && m->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    *curv = norm_curv[2]; 
    return M3DC1_SUCCESS;
  }

  *curv=0.0;
  gmi_ent* gent= (gmi_ent*)(m->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if(gType==0)
  {
    apf::Vector3 vcd_t;
    double vcd[3];
    m->getPoint(vt, 0, vcd_t);
    for(int i=0; i<3; i++)
      vcd[i]=vcd_t[i];
    std::vector<gmi_ent*> gEdges;
    get_gv_edges(gent, gEdges);
    int numEdgesPlane=0;
    double curv_tmp;
    for(int i=0;i<gEdges.size();i++)
    {
      gmi_ent* pe = gEdges.at(i);
      double cd[3]={0,0,0};
      double paraRange[2];
      gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
      if(!pn) continue;
      numEdgesPlane++;
      evalCoord(paraRange[0], cd, pn);
      if(checkSamePoint2D(vcd,cd))
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
  else if(gType==1)
  {
      apf::Vector3 param(0,0,0);
      m->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalCurvature(pn[0],pn[1], param[0], curv);
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh; 

  apf::MeshEntity* vt = getMdsEntity(m, 0, *node_id);
  assert(vt);
  gmi_ent* gent= (gmi_ent*)(m->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  *on_geom_bdry=(gType==0||gType==1);
  return M3DC1_SUCCESS;
}

// region-specific function
//*******************************************************
int m3dc1_region_getoriginalface( int * /* in */ elm, int * /* out */ fac)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh; 

  apf::MeshEntity* ent = getMdsEntity(m, 3, *elm);
  apf::Downward downward;
  int num_adj_ent = m->getDownward(ent, 2, downward);
  assert(num_adj_ent==5);
  int triFace[2];
  int counter=0;
  for(int i=0; i<num_adj_ent; i++)
  {
    int num_adj_ent;
    apf::Downward downward2;
    int num_edge= m->getDownward(downward[i], 1, downward2);
    if(num_edge==3) triFace[counter++]= getMdsIndex(m,downward[i]);    
  }
  assert(counter==2);
  *fac = std::min(triFace[0],triFace[1]);
  return M3DC1_SUCCESS;
}

/** field manangement */
int fieldIdMax=0;
//*******************************************************
int m3dc1_field_getnewid ( FieldID* /*out*/field_id )
//*******************************************************
{
  *field_id = fieldIdMax+1;
  return M3DC1_SUCCESS;
}

// *scalar_type is either M3DC1_REAL or M3DC1_COMPLEX
//*******************************************************
int m3dc1_field_create (FieldID* /*in*/ field_id, 
        const char* /* in */ field_name, int* /*in*/ num_values, 
        int* /*in*/ scalar_type, int* /*in*/ num_dofs_per_value)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  apf::Field* f = NULL;
  int components = (*num_values)*(*scalar_type+1)*(*num_dofs_per_value);

  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
  {
    if (!m3dc1_mesh::instance()->field_container)
      m3dc1_mesh::instance()->field_container=new std::map<FieldID, m3dc1_field*>;
    m = m3dc1_mesh::instance()->mesh;
    f = createPackedField(m, field_name, components);
    m3dc1_mesh::instance()->field_container->insert(
          std::map<FieldID, m3dc1_field*>::value_type(*field_id,
          new m3dc1_field(*field_id, f, *num_values, *scalar_type,*num_dofs_per_value)));

  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container)
      m3dc1_ghost::instance()->field_container=new std::map<FieldID, m3dc1_field*>;
    m = m3dc1_ghost::instance()->mesh;
    f = createPackedField(m, field_name, components);
    m3dc1_ghost::instance()->field_container->insert(
          std::map<FieldID, m3dc1_field*>::value_type(*field_id,
          new m3dc1_field(*field_id, f, *num_values, *scalar_type, *num_dofs_per_value)));
  }

  apf::freeze(f); // switch dof data from tag to array

#ifdef DEBUG
  if (!PCU_Comm_Self() && *field_id!=NODE_GLB_ORDER) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<", #values "<<*num_values<<", #dofs "<<countComponents(f)<<"\n";
#endif

  if (*field_id>fieldIdMax) fieldIdMax=*field_id;
  double val[2]={0,0};
  m3dc1_field_assign(field_id, val, scalar_type);
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_field_delete (FieldID* /*in*/ field_id)
//*******************************************************
{
  apf::Field* f = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_ghost::instance()->field_container)[*field_id]->get_field();
  }    

#ifdef DEBUG
  if (!PCU_Comm_Self() && *field_id!=NODE_GLB_ORDER) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<"\n";
#endif

  destroyField(f);

  // remove f from field container
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    delete (*m3dc1_mesh::instance()->field_container)[*field_id];
    m3dc1_mesh::instance()->field_container->erase(*field_id);
  }
  else {
    delete (*m3dc1_ghost::instance()->field_container)[*field_id];
    m3dc1_ghost::instance()->field_container->erase(*field_id);
  }

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getinfo(FieldID* /*in*/ field_id, 
                        char* /* out*/ field_name, int* num_values, 
                        int* scalar_type, int* total_num_dof)
//*******************************************************
{
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else 
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    

  assert(mf);

  apf::Field* f = mf ->get_field();
  strcpy(field_name, getName(f));
  *num_values = mf -> get_num_value();
  *scalar_type = mf ->get_value_type();
  *total_num_dof = countComponents(f);
  if (*scalar_type) *total_num_dof/=2;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_exist(FieldID* field_id, int * exist)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      *exist = 0;
    else
      *exist = 1;
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      *exist = 0;
    else
      *exist = 1;
  }
  
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_field_synchronize(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else 
    m = m3dc1_ghost::instance()->mesh;

  apf::MeshEntity* e;       

  int num_dof, n = countComponents(f);
  double* sender_data = new double[n];
  double* dof_data = new double[n]; 

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_owned(m, e) || !m->isShared(e)) continue;
    getComponents(f, e, 0, dof_data);

    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
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
      setComponents(f, r, 0, sender_data);
    }
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_field_sync (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  apf::Field* f = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_ghost::instance()->field_container)[*field_id]->get_field();
  }    

  m3dc1_field_synchronize(f);

#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

// send non-owned copies' dof to owner copy and add them up
//*******************************************************
void m3dc1_field_accumulate(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else 
    m = m3dc1_ghost::instance()->mesh;
  apf::MeshEntity* e;       

  int num_dof, own_partid, n = countComponents(f);
  double* dof_data = new double[n];
  double* sender_data = new double[n];
  apf::MeshEntity* own_e;
  apf::MeshEntity* r;
  std::map<apf::MeshEntity*, std::vector<double> > save_map;

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid=get_ent_ownpartid(m, e);
    if (own_partid==PCU_Comm_Self()) continue;

    own_e = get_ent_owncopy(m, e);

    getComponents(f, e, 0, &(dof_data[0]));
      
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(dof_data[0]),n*sizeof(double));
  }
  m->end(it);

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
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
    setComponents(f, e, 0, dof_data);
  } 
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_field_sum (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  apf::Field* f = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    f = (*m3dc1_ghost::instance()->field_container)[*field_id]->get_field();
  } 

  m3dc1_field_accumulate(f);
  m3dc1_field_synchronize(f);
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_sumsq (FieldID* /* in */ field_id, double* /* out */ sum)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  apf::Field* f = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    m = m3dc1_mesh::instance()->mesh;
    f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    m = m3dc1_ghost::instance()->mesh;
    f = (*m3dc1_ghost::instance()->field_container)[*field_id]->get_field();
  }    
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  *sum=0.;
  int num_dof = countComponents(f);

  double* dof_data= new double[num_dof];
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_owned(m,e)) continue;
    getComponents(f, e, 0, dof_data);
    for(int i=0; i<num_dof; ++i)
      *sum+=dof_data[i]*dof_data[i];
  }
  m->end(it);
  delete [] dof_data;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_sum_plane (FieldID* /* in */ field_id)
//*******************************************************
{
  MPI_Comm icomm= m3dc1_model::instance()->getMPICommPlane();
  int num_vtx,num_dof=0, vertex_type=0;
  m3dc1_field_getnumlocaldof(field_id, &num_dof);
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);

  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

  int data_size=total_num_dof*num_vtx*(1+value_type);
  assert(total_num_dof*num_vtx==num_dof);
  double * thevec = NULL;
  m3dc1_field_getdataptr(field_id, &thevec);
  double * sendbuf = new double [data_size];
  m3dc1_field_retrieve (field_id, sendbuf, &num_dof);
  MPI_Allreduce (sendbuf,thevec,data_size,MPI_DOUBLE,MPI_SUM,icomm) ;
  m3dc1_field_sync(field_id);
  delete []sendbuf;
}

/** field dof functions */
//*******************************************************
int m3dc1_field_getlocaldofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  }
  else 
  {
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    
  } 
    
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  *start_dof_id=0;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    *end_dof_id_plus_one=num_dof*m3dc1_mesh::instance()->num_local_ent[0];
  else
    *end_dof_id_plus_one=num_dof*m3dc1_ghost::instance()->num_local_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getowndofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  // ghosting doesn't change the ownership therefore this routine is the same with/without ghosting
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0];
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  
  int start_id = num_own_ent;
  PCU_Exscan_Ints(&start_id,1);

  *start_dof_id=start_id*num_dof;
  *end_dof_id_plus_one=*start_dof_id+num_own_ent*num_dof;
  return M3DC1_SUCCESS;
}
 
//******************************************************* 
int m3dc1_field_getglobaldofid ( FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  // ghosting doesn't change the ownership therefore this routine is the same with/without ghosting
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  assert(mf->get_num_value()*mf->get_dof_per_value()==num_dof);  

  *start_dof_id=0;
  *end_dof_id_plus_one=*start_dof_id+num_dof*m3dc1_mesh::instance()->num_global_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumlocaldof (FieldID* field_id, int* /* out */ num_local_dof)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
#ifdef DEBUG
    if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
#endif
    m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
    *num_local_dof = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  }
  else
  {
#ifdef DEBUG
    if (!m3dc1_ghost::instance()->field_container || !m3dc1_ghost::instance()->field_container->count(*field_id))
      return M3DC1_FAILURE;
#endif
    m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
    *num_local_dof = (m3dc1_ghost::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumowndof (FieldID* field_id, int* /* out */ num_own_dof)
//*******************************************************
{
  // ghosting doesn't change the ownership therefore this routine is the same with/without ghosting
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_own_dof = (m3dc1_mesh::instance()->num_own_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumglobaldof (FieldID* field_id, int* /* out */ num_global_dof)
//*******************************************************
{
  // ghosting doesn't change the ownership therefore this routine is the same with/without ghosting
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_global_dof = (m3dc1_mesh::instance()->num_global_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getdataptr (FieldID* field_id, double** pts)
//*******************************************************
{
  apf::Field* f=NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  else
    f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  if (!isFrozen(f)) freeze(f);
  *pts=getArrayData(f);
  return M3DC1_SUCCESS;
}

// add field2 to field1
//*******************************************************
int m3dc1_field_add(FieldID* /*inout*/ field_id1, FieldID* /*in*/ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif

  m3dc1_field *mf1 = NULL, *mf2 = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {  
    mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
    mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  }
  else 
  {
    mf1 = (*(m3dc1_ghost::instance()->field_container))[*field_id1];
    mf2 = (*(m3dc1_ghost::instance()->field_container))[*field_id2];
  }
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value(); 
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)+=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  } 
#ifdef DEBUG
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_mult(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  double dofs[FIXSIZEBUFF], dofsNew[FIXSIZEBUFF];
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  int value_type = mf->get_value_type();
  assert(dofPerEnt<=sizeof(dofs)/sizeof(double)*(1+value_type));
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    if(*scalar_type==0)
    {
      for(int i=0; i<dofPerEnt*(1+value_type); i++)
        dofsNew[i]=*fac*dofs[i];
    }
    else
    {
      for(int i=0; i<dofPerEnt; i++)
      {
        dofsNew[2*i]=fac[0]*dofs[2*i]-fac[1]*dofs[2*i+1];
        dofsNew[2*i+1]=fac[0]*dofs[2*i+1]+fac[1]*dofs[2*i];
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofsNew[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_assign(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()), fac[0]);
  if(*scalar_type)
    for(int i=0; i<dofPerEnt; i++)
      dofs.at(2*i+1)=fac[1];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
  }
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_copy(FieldID* /* out */ field_id1, FieldID* /* in */ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif

  m3dc1_field *mf1 = NULL, *mf2 = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
    mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  }
  else {
    mf1 = (*(m3dc1_ghost::instance()->field_container))[*field_id1];
    mf2 = (*(m3dc1_ghost::instance()->field_container))[*field_id2];
  }
  apf::Field* f1 =  mf1->get_field();
  apf::Field* f2 =  mf2->get_field();
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value();
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_retrieve (FieldID* /* in */ field_id, double * /*out*/ data, int * /* in */size)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(data, pts, *size*(1+value_type)*sizeof(double));
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_set (FieldID* /* in */ field_id, double * /*in*/ data, int * /* in */size)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(pts, data, *size*(1+value_type)*sizeof(double));
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_insert(FieldID* /* in */ field_id, int /* in */ * local_dof, 
         int * /* in */ size, double* /* in */ values, int * type, int * op)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
#ifdef DEBUG
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);
  assert(*local_dof<num_local_dof);
  if(!value_type) assert(!(*type)); // can not insert complex value to real vector
  for(int i=0; i<*size*(1+(*type)); i++)
  {
#ifdef REPLACENANWITHZERO
    if(values[i]!=values[i]) values[i]=0;
#else
    assert(values[i]==values[i]);
#endif
  }
#endif
  std::vector<double> values_convert(*size*(1+value_type),0);
  if(!(*type)&&value_type) // real into complex
  {
    for(int i=0; i<*size; ++i)
      values_convert.at(2*i)=values[i];
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); ++i)
      values_convert.at(i)=values[i];
  }
  double * dataptr;
  int ibegin=*local_dof*(1+value_type);
  m3dc1_field_getdataptr(field_id, &dataptr);
  if(*op==0) // set value
  {
   for(int i=0; i<*size*(1+value_type); ++i)
     dataptr[ibegin+i]=values_convert.at(i);
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); ++i)
      dataptr[ibegin+i]+=values_convert[i];
  }
  return M3DC1_SUCCESS;
}
#define FIELDVALUELIMIT 1e100
bool value_is_nan(double val)
{
  return val!=val ||fabs(val) >FIELDVALUELIMIT;
}

//*******************************************************
int m3dc1_field_isnan(FieldID* /* in */ field_id, int * isnan)
//*******************************************************
{
  *isnan=0;
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEnt;
  double dofs[FIXSIZEBUFF];
  for(int inode=0; inode<num_vtx; ++inode)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    for(int i=0; i<dofPerEnt; i++)
      if(value_is_nan(dofs[i])) 
        *isnan=1;
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
//int m3dc1_field_write(FieldID* field_id, const char* filename)
//*******************************************************
/*{
  if (!filename) 
  {
    m3dc1_field_print(field_id);
    return M3DC1_SUCCESS;
  }
   
  apf::Mesh2* m = NULL;
  apf::Field* f=NULL;
  
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    m = m3dc1_mesh::instance()->mesh;
    f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  }
  else 
  {
    m = m3dc1_ghost::instance()->mesh;
    f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  }

  double* field_data =getArrayData(f);

  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);

  char field_filename[256];
  sprintf(field_filename,"%s-%d-%d",filename, *field_id, PCU_Comm_Self());
  FILE * fp =fopen(field_filename, "w");
  
  apf::MeshEntity* e;
  int num_dof=countComponents(f);
  double* dof_data = new double[num_dof];

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    switch (num_dof)
    {
      case 1: {
        getComponents(f, e, 0, dof_data);
	fprintf(fp, "[p %d] field %d/ent %d: [%lf]\n", PCU_Comm_Self(), *field_id, 
                get_ent_globalid(m,e,0), dof_data[0]);
        break;
      }
      case 2: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf]\n", PCU_Comm_Self(), *field_id, 
                get_ent_globalid(m,e,0), dof_data[0], dof_data[1]);
        break;
      }
      case 3: {
        getComponents(f, e, 0, dof_data);
	fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf]\n", PCU_Comm_Self(), *field_id, 
                get_ent_globalid(m,e,0), dof_data[0], dof_data[1],dof_data[2]);
        break;
      }
      case 4: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf]\n", PCU_Comm_Self(), *field_id,  
                get_ent_globalid(m,e,0), dof_data[0], dof_data[1],dof_data[2],dof_data[3]);
        break;
      }
      case 6: {
        getComponents(f, e, 0, dof_data);
  	fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf, %lf, %lf]\n", PCU_Comm_Self(), *field_id, 
                get_ent_globalid(m,e,0),
                dof_data[0], dof_data[1],dof_data[2],dof_data[3],dof_data[4],dof_data[5]);
        break;
      }
      case 8: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf]\n", PCU_Comm_Self(), *field_id,
               get_ent_globalid(m,e,0), 
               dof_data[0], dof_data[1],dof_data[2],dof_data[3],dof_data[4],dof_data[5],dof_data[6], dof_data[7]);
        break;
      }
      case 12: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf]\n", 
                   PCU_Comm_Self(), *field_id, get_ent_globalid(m,e,0), 
                   dof_data[0],dof_data[1],dof_data[2],dof_data[3],dof_data[4],dof_data[5],dof_data[6],dof_data[7],
                   dof_data[8],dof_data[9],dof_data[10],dof_data[11]);
        break;
      }
      case 18: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf]\n", 
                   PCU_Comm_Self(), *field_id, get_ent_globalid(m,e,0), 
                   dof_data[0],dof_data[1],dof_data[2],dof_data[3],dof_data[4],dof_data[5],dof_data[6],dof_data[7],
                   dof_data[8],dof_data[9],dof_data[10],dof_data[11],dof_data[12],dof_data[13],dof_data[14],
                   dof_data[15], dof_data[16],dof_data[17]);  
        break;
      }
      case 24: {
        getComponents(f, e, 0, dof_data);
        fprintf(fp, "[p %d] field %d/ent %d: [%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf]\n", 
        PCU_Comm_Self(), *field_id, get_ent_globalid(m,e,0), 
        dof_data[0],dof_data[1],dof_data[2],dof_data[3],dof_data[4],dof_data[5],dof_data[6],dof_data[7],
        dof_data[8],dof_data[9],dof_data[10],dof_data[11],dof_data[12],dof_data[13],dof_data[14],dof_data[15],
        dof_data[16],dof_data[17],dof_data[18],dof_data[19],dof_data[20],dof_data[21],dof_data[22],dof_data[23]);
        break;
      }
      default: if (!PCU_Comm_Self()) std::cout<<__func__<<" failed for field "<<getName(f)<<": does support "<<num_dof<<" dofs\n";
    } // switch
  } // while
  m->end(it);
  fclose(fp);  
  return M3DC1_SUCCESS;
}
*/
//=========================================================================
void write_vector(apf::Mesh2* m, apf::Field* f, const char* filename, int start_index)
{
  char field_filename[256];
  sprintf(field_filename,"%s-%d",filename, PCU_Comm_Self());
  FILE * fp =fopen(field_filename, "w");

  apf::MeshEntity* e;
  int num_dof=countComponents(f);
  double* dof_data = new double[num_dof];

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    // write vector
    getComponents(f, e, 0, dof_data);
    for (int i=0; i<num_dof; ++i)
      fprintf(fp, "%d\t%E\n", get_ent_globalid(m,e,0)+start_index, dof_data[i]);
  } // while
  m->end(it);
  fclose(fp);  
}

//*******************************************************
int m3dc1_field_write (FieldID* field_id, const char* filename, int* start_index)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  apf::Field* f=NULL;
  
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    m = m3dc1_mesh::instance()->mesh;
    f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  }
  else 
  {
    m = m3dc1_ghost::instance()->mesh;
    f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  }
  write_vector(m, f, filename, *start_index);
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_field_compare(FieldID* field_id_1, FieldID* field_id_2)
//*******************************************************
{
  int ierr = M3DC1_SUCCESS, it_num_local_ent;
  apf::Field *f_1 = NULL, *f_2 = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
  {
    it_num_local_ent = m3dc1_mesh::instance()->num_local_ent[0];
    f_1 = (*(m3dc1_mesh::instance()->field_container))[*field_id_1]->get_field();
    f_2 = (*(m3dc1_mesh::instance()->field_container))[*field_id_2]->get_field();
  }
  else 
  {
    it_num_local_ent = m3dc1_ghost::instance()->num_local_ent[0];
    f_1 = (*(m3dc1_ghost::instance()->field_container))[*field_id_1]->get_field();
    f_2 = (*(m3dc1_ghost::instance()->field_container))[*field_id_2]->get_field();
  }    
  double* field_data_1 =getArrayData(f_1);
  double* field_data_2 =getArrayData(f_2);

  int num_dof_1=countComponents(f_1);
  int num_dof_2=countComponents(f_2);
  if (num_dof_1!=num_dof_2) 
  {
    if (!PCU_Comm_Self()) 
      cout<<"[M3D-C1 INFO] "<<__func__<<": #dof mismatch "<<getName(f_1)
          <<"- "<<num_dof_1<<", "<<getName(f_2)<<"- "<<num_dof_2<<"\n";
    return M3DC1_FAILURE;
  }

  for (int i=0; i<num_dof_1*it_num_local_ent; ++i)
  {  
    if (!m3dc1_double_isequal(field_data_1[i], field_data_2[i])) 
    {
     cout<<"[M3D-C1 ERROR] "<<__func__<<": "<<getName(f_1)<<"["<<i<<"]="<<field_data_1[i]
          <<", "<<getName(f_2)<<"["<<i<<"]="<<field_data_2[i]<<"\n";
      ierr=M3DC1_FAILURE;
      break;
    }
  }
  int global_ierr;
  MPI_Allreduce(&ierr, &global_ierr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  if (global_ierr==M3DC1_FAILURE)
  {
    if (!PCU_Comm_Self())
      cout<<"[M3D-C1 INFO] "<<__func__<<": dof value mismatch of fields "<<getName(f_1)
          <<" and "<<getName(f_2)<<"\n";
    
    return M3DC1_FAILURE;
  }
  else
    return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_print(FieldID* field_id)
//*******************************************************
{
  apf::Mesh2* m = NULL;
  apf::Field* f=NULL;
  
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
  {
    m = m3dc1_mesh::instance()->mesh;
    f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  }
  else 
  {
    m = m3dc1_ghost::instance()->mesh;
    f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  }

  
  if (!f)
  {
    if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<" failed as field "<<*field_id<<" not found\n";
    return M3DC1_FAILURE;
  }

  double* field_data =getArrayData(f);

  if (!m->findField("node global id field"))
  {
    if (!PCU_Comm_Self()) cout<<"[M3D-C1 INFO] "<<__func__<<" failed as global node id not found\n";
    return M3DC1_FAILURE;
  }

  apf::MeshEntity* e;
  int num_dof=countComponents(f);
  double* dof_data = new double[num_dof];
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    switch (num_dof)
    {
      case 1: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<"]\n";
        break;}
      case 2: {     
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<"]\n";
        break;}
      case 3: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<"]\n";
        break;}
    case 4: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<"]\n";
 
        break; }
      case 6: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<"]\n";
        break; }
      case 8: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<"]\n";
        break; }
      case 12: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<"]\n";
        break; }
      case 18: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<"]\n";
        break; }
      case 24: {
        getComponents(f, e, 0, dof_data);
	std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
		     <<"/ent "<<get_ent_globalid(m, e, 0)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<", "<<dof_data[18]
		     <<", "<<dof_data[19]
		     <<", "<<dof_data[20]
		     <<", "<<dof_data[21]
		     <<", "<<dof_data[22]
		     <<", "<<dof_data[23]
		     <<"]\n";
        break; }
      default: if (!PCU_Comm_Self()) std::cout<<__func__<<" failed for field "
               <<getName(f)<<": does support "<<num_dof<<" dofs\n";
               break;
    } // switch
  } // while
  m->end(it);
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                       int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::MeshEntity* e = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  else
    e = getMdsEntity(m3dc1_ghost::instance()->mesh, *ent_dim, *ent_id);
  
  assert(e);

  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else 
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    
  assert(mf);
   
  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;
  *start_dof_id = *ent_id*dof_per_node;
  *end_dof_id_plus_one = *start_dof_id +dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
         int* /* out */ start_global_dof_id, int* /* out */ end_global_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::Mesh2* m=NULL;
  apf::MeshEntity* e =NULL;
  m3dc1_field * mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    m = m3dc1_mesh::instance()->mesh;
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  }
  else
  {
    m = m3dc1_ghost::instance()->mesh;
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  }
  e = getMdsEntity(m, *ent_dim, *ent_id);
  assert(e);

  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;

  *start_global_dof_id = (get_ent_globalid(m, e, 0))*dof_per_node;
  *end_global_dof_id_plus_one =*start_global_dof_id + dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e = NULL;
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
    e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  }
  else 
  {
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    
    e = getMdsEntity(m3dc1_ghost::instance()->mesh, *ent_dim, *ent_id);
  } 
  assert(e && mf);     
  apf::Field* f = mf->get_field();

#ifdef DEBUG
  assert(countComponents(f)==*num_dof*(1+mf->get_value_type()));
  for(int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
#endif
  setComponents(f, e, 0, dof_data);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e = NULL;  
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0])
  {
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
    e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  }
  else 
  {
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    
    e = getMdsEntity(m3dc1_ghost::instance()->mesh, *ent_dim, *ent_id);
  } 

  assert(e && mf);
  
  apf::Field* f = mf->get_field();

  getComponents(f, e, 0, dof_data);
  *num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
#ifdef DEBUG
  for(int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
  int start_dof_id,end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid(ent_dim, ent_id,field_id, &start_dof_id, &end_dof_id_plus_one);
  double* data;
  m3dc1_field_getdataptr(field_id, &data);
  int start=start_dof_id*(1+mf->get_value_type());
  for( int i=0; i< *num_dof; i++)
    assert(data[start++]==dof_data[i]);
#endif
  return M3DC1_SUCCESS;
}

#ifdef M3DC1_PETSC
/** matrix and solver functions */
std::map<int, int> matHit;
int getMatHit(int id) { return matHit[id];};
void addMatHit(int id) { matHit[id]++; }

//*******************************************************
int m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID *field_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }

  if (*matrix_type==M3DC1_MULTIPLY) // matrix for multiplication
  {
    matrix_mult* new_mat = new matrix_mult(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }
  else 
  {
    matrix_solve* new_mat= new matrix_solve(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_freeze(int* matrix_id) 
//*******************************************************
{
  double t1 = MPI_Wtime();
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->assemble();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_delete(int* matrix_id)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  typedef std::map<int, m3dc1_matrix*> matrix_container_map;
  m3dc1_solver::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  delete mat;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_insert(int* matrix_id, int* row, 
         int* col, int* scalar_type, double* val)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_status()==M3DC1_FIXED)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" is fixed\n";
    return M3DC1_FAILURE;
  }
  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, INSERT_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, INSERT_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_add (int* matrix_id, int* row, int* col, 
                      int* scalar_type, double* val) //globalinsertval_
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_status()==M3DC1_FIXED)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" is fixed\n";
    return M3DC1_FAILURE;
  }
  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, ADD_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, ADD_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_setbc(int* matrix_id, int* row)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  int row_g = start_global_dof_id+*row%total_num_dof;
  (dynamic_cast<matrix_solve*>(mat))->set_bc(row_g);
}

//*******************************************************
int m3dc1_matrix_setlaplacebc(int * matrix_id, int *row,
         int * numVals, int *columns, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_SOLVE)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  std::vector <int> columns_g(*numVals);
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif
  int row_g = start_global_dof_id+*row%total_num_dof;
  for(int i=0; i<*numVals; i++)
  {
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;
  }
  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
}

int m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol, int* solver_type, double* solver_tol) //solveSysEqu_
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<*rhs_sol<<"\n";
  // solver_type: 0 for iterative solver, 1 for direct solver.
  (dynamic_cast<matrix_solve*>(mat))->solve(*rhs_sol, *solver_type, *solver_tol);

  addMatHit(*matrix_id);
}

//*******************************************************
int m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, 
         FieldID* outputvecid) 
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_MULTIPLY)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for solving (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
#ifdef DEBUG_
  int isnan;
  if(!PCU_Comm_Self()) std::cout<<" mult matrix "<<*matrix_id<<" with input vec "<<*inputvecid<<" output vec "<<*outputvecid<<std::endl;
  m3dc1_field_isnan(inputvecid, &isnan);
  assert(isnan==0);
#endif

#ifdef PRINTSOLVEMATRIX
  char filename[256];
  int time=getMatHit(*matrix_id);
  sprintf(filename, "mat%din%d.m",*matrix_id, time);
  int zero=0;
  m3dc1_field_write(inputvecid, filename,&zero);
#endif
  (dynamic_cast<matrix_mult*>(mat))->multiply(*inputvecid, *outputvecid);
#ifdef PRINTSOLVEMATRIX
  sprintf(filename, "mat%dout%d.m",*matrix_id, time);
  m3dc1_field_write(outputvecid, filename,&zero);
#endif
#ifdef DEBUG_
  m3dc1_field_isnan(outputvecid, &isnan);
  assert(isnan==0);
#endif
  addMatHit(*matrix_id);
}

//*******************************************************
int m3dc1_matrix_flush(int* matrix_id)
//*******************************************************
{
  double t1 = MPI_Wtime();
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->flushAssembly();
}

//*******************************************************
int m3dc1_matrix_getiternum(int* matrix_id, int * iter_num)
//*******************************************************
{ 
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  *iter_num = dynamic_cast<matrix_solve*> (mat)->iterNum;
}

//*******************************************************
int m3dc1_matrix_insertblock(int* matrix_id, int * ielm, 
          int* rowIdx, int * columnIdx, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  int field = mat->get_fieldOrdering();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = 2;
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  if (m3dc1_mesh::instance()->mesh->getDimension()==3) ielm_dim =3;
  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;
  int start_global_dof,end_global_dof_id;
  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowIdx<numVar && *columnIdx<numVar);
  int rows[1024], columns[1024];
  assert(sizeof(rows)/sizeof(int)>=dofPerVar*nodes_per_element);
  if(mat->get_type()==0)
  {
    int localFlag=0;
    matrix_mult* mmat = dynamic_cast<matrix_mult*> (mat);
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      if(mmat->is_mat_local()) m3dc1_ent_getlocaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      else m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
  }
  else
  {
    matrix_solve* smat = dynamic_cast<matrix_solve*> (mat);
    int nodeOwner[6];
    int columns_bloc[6], rows_bloc[6];
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnIdx;
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      if(nodeOwner[inode]!=PCU_Comm_Self()&&!m3dc1_solver::instance()->assembleOption)
        smat->add_blockvalues(1, rows_bloc+inode, nodes_per_element, columns_bloc, values+offset);
      else 
        smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      offset+=numValuesNode;
    }
  }
}

//*******************************************************
int m3dc1_matrix_write(int* matrix_id, const char* filename, int* start_index)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (!filename) 
    return m3dc1_matrix_print(matrix_id);

  char matrix_filename[256];
  sprintf(matrix_filename,"%s-%d",filename, PCU_Comm_Self());
  FILE * fp =fopen(matrix_filename, "w");

  int row, col, csize, sum_csize=0, index=0;

  vector<int> rows;
  vector<int> n_cols;
  vector<int> cols;
  vector<double> vals;

  mat->get_values(rows, n_cols, cols, vals);
  for (int i=0; i<rows.size(); ++i)
    sum_csize += n_cols[i];
  assert(vals.size()==sum_csize);

  fprintf(fp, "%d\t%d\t%d\n", rows.size(), n_cols.size(), vals.size());

  for (int i=0; i<rows.size(); ++i)
  {
    row = rows[i];
    csize = n_cols[i];
    for (int j=0; j<csize; ++j)
    {
      fprintf(fp, "%d\t%d\t%E\n", row+*start_index, cols[index]+*start_index,vals[index]);
      ++index;
    }
  }
  fclose(fp);  
  assert(index == vals.size());
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_matrix_print(int* matrix_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  int row, col, csize, sum_csize=0, index=0;

  vector<int> rows;
  vector<int> n_cols;
  vector<int> cols;
  vector<double> vals;

  mat->get_values(rows, n_cols, cols, vals);
  for (int i=0; i<rows.size(); ++i)
    sum_csize += n_cols[i];
  assert(vals.size()==sum_csize);

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": printing matrix "<<*matrix_id<<"\n";

  for (int i=0; i<rows.size(); ++i)
  {
    row = rows[i];
    csize = n_cols[i];
    for (int j=0; j<csize; ++j)
    {
      std::cout<<"["<<PCU_Comm_Self()<<"]\t"<<row<<"\t"<<cols[index]<<"\t"<<vals[index]<<"\n";
      ++index;
    }
  }
  assert(index == vals.size());
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_matrix_setassembleoption(int * op)
//*******************************************************
{
  m3dc1_solver::instance()->assembleOption= *op;
}

#else
#include "m3dc1_ls.h"
#endif // #ifdef M3DC1_PETSC

int adapt_time=0;
void group_complex_dof (apf::Field* field, int option);

int adapt_by_field (int * fieldId, double* psi0, double * psil)
{
  FILE *fp = fopen("sizefieldParam", "r");
  if(!fp)
  {
    std::cout<<" file sizefieldParam not found "<<std::endl;
    throw 1;
  }
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  // read the size field parameters
  double param[13];
  for(int i=0; i<13; i++)
    fscanf(fp, "%lf ", &param[i]);
  fclose(fp);

  // delete global numbering
  destroy_global_numbering(mesh);

  // delete all matrices
#ifdef M3DC1_TRILINOS
  std::map<int, m3dc1_epetra*>::iterator mat_it = m3dc1_ls::instance()->matrix_container->begin();
  while (m3dc1_ls::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_epetra*>::iterator mat_it = m3dc1_ls::instance()->matrix_container->begin();
    mat_it->second->destroy();
    delete mat_it->second;
    m3dc1_ls::instance()->matrix_container->erase(mat_it);
  }
#endif
#ifdef M3DC1_PETSC
  std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }
#endif
  int valueType = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_value_type();
  apf::Field* psiField = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_field();
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);
  double mmax[2], mmin[2];
  m3dc1_model_getmaxcoord(mmax,mmax+1);
  m3dc1_model_getmincoord(mmin,mmin+1);

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    assert(valueType==complexType);
    if(complexType) group_complex_dof(field, 1);
    if(isFrozen(field)) unfreeze(field);
    if(!PCU_Comm_Self()) std::cout<<__func__<<": transferring field "<<apf::getName(field)<<"\n";
    fields.push_back(field);
    it++;
  }

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<__func__<<": deleting numbering "<<getName(n)<<endl;
    apf::destroyNumbering(n);
  }

  char filename[256];
  sprintf(filename,"before-adapt-%d",adapt_time);
  apf::writeVtkFiles(filename,mesh);

  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
#ifdef OLDMA
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
#else
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh,&sf,&slnTrans));
#endif
  in->maximumIterations = 9;
  in->shouldRunPostZoltan = true;
  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] running MeshAdapt...\n";
  ma::adapt(in);
  reorderMdsMesh(mesh);

  sprintf(filename,"after-adapt-%d",adapt_time);
  apf::writeVtkFiles(filename,mesh);
  ++adapt_time;

  m3dc1_mesh::instance()->initialize();
  generate_global_numbering(mesh); // new global numbering

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 0);
    if(!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    m3dc1_field_synchronize(field);
#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
}


double absSize[2]={0,1}, relSize[2]={0.3, 1.5};
int set_mesh_size_bound (double* abs_size, double * rel_size)
{
  for(int i=0; i<2; i++)
  {
    absSize[i] =  abs_size[i];
    relSize[i] = rel_size[i];
  }
}

void smooth_size_field (apf::Field* sizeField)
{
  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  for(int i=0; i<numVert; i++)
  {
    vector<apf::MeshEntity*> nodes;
    apf::MeshEntity* e = apf:: getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    nodes.push_back(e);
    double sizeOrg=0;
    getComponents(sizeField, e, 0, &sizeOrg);
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,1,adjacent);
    for(int i=0; i<adjacent.getSize(); i++)
    {
      apf::Downward downward;
      m3dc1_mesh::instance()->mesh->getDownward(adjacent[i], 0, downward);
      nodes.push_back(downward[0]==e?downward[1]:downward[0]);
    }
    double size=0;
    for(int i=0; i<nodes.size(); i++)
    {
      double buff;
      getComponents(sizeField, nodes[i], 0, &buff);
      size+=1./buff;
    }
    size/=nodes.size();
    size=1./size;
    setComponents(sizeField, e, 0, &size);
  }
}

void group_complex_dof (apf::Field* field, int option)
{
  int num_dof_double = countComponents(field);
  assert(num_dof_double/6%2==0);
  int num_dof = num_dof_double/2;
  vector<double> dofs(num_dof_double);
  vector<double> newdofs(num_dof_double);
  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    getComponents(field, e, 0, &(dofs[0]));
    for(int j=0; j<num_dof/6; j++)
    {
      if(option)
      {
        for(int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+k)=dofs.at(2*j*6+2*k);
          newdofs.at(2*j*6+6+k)=dofs.at(2*j*6+2*k+1);
        }
      }
      else
      {
        for(int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+2*k)=dofs.at(2*j*6+k);
          newdofs.at(2*j*6+2*k+1)=dofs.at(2*j*6+6+k);
        }
      }
    }
    setComponents(field, e, 0, &(newdofs[0]));
  }
}

double p=4;
int set_adapt_p (double * pp) {p=*pp;}

int adapt_by_error_field (double * errorData, double * errorAimed, int * max_adapt_node, int * option)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::Field* sizeField = createPackedField(m3dc1_mesh::instance()->mesh, "size_field", 1);
  SizeFieldError sf (m3dc1_mesh::instance()->mesh, sizeField, *errorAimed);

  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  // first sum error ^ (2d/(2p+d))
  double d=2;
  double errorSum=0;
  if(*option)
  {
    for(int i=0; i<numVert; i++)
    {
      if(is_ent_owned(m3dc1_mesh::instance()->mesh,getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i)))
        errorSum+=pow(errorData[i],d/(p+d/2.0));
    }
    double errorSumBuff=errorSum;
    MPI_Allreduce(&errorSumBuff, &errorSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errorSum = *errorAimed*(*errorAimed)/errorSum;
    errorSum = pow(errorSum,1./(2.*p));
  }
  else errorSum=pow(*errorAimed,1./(p+d/2.));
  double size_estimate=0;
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    if(!is_ent_owned(m3dc1_mesh::instance()->mesh,e)) continue;
    double size = sf.getSize(e);
    double targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    size_estimate+=max(1.,1./targetSize/targetSize);
  }
  double size_estimate_buff=size_estimate;
  MPI_Allreduce(&size_estimate_buff, &size_estimate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int numNodeGlobl=0, dim=0;
  m3dc1_mesh_getnumglobalent(&dim, &numNodeGlobl);
  cout<<" numVert "<<numNodeGlobl<<" size_estimate "<<size_estimate;
  if(size_estimate>*max_adapt_node) errorSum*=sqrt(size_estimate/(*max_adapt_node));
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
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
  // only implemented for one process
  //if(PCU_Comm_Peers()==1)
  {
    smooth_size_field(sizeField);
    smooth_size_field(sizeField);
  }
  m3dc1_field_synchronize(sizeField);

  // delete global numbering
  destroy_global_numbering(mesh);

  // delete all matrices
#ifdef M3DC1_PETSC
  std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }
#endif
  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  fields.push_back(sizeField);
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 1);
    if(isFrozen(field)) unfreeze(field);
    fields.push_back(field);
    it++;
  }

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<__func__<<": deleting numbering "<<getName(n)<<endl;
    apf::destroyNumbering(n);
  }

  char filename[256];
  sprintf(filename,"before-adapt-%d",adapt_time);
  apf::writeVtkFiles(filename,mesh);

  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
#ifdef OLDMA
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
#else
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh,&sf,&slnTrans));
#endif
  in->maximumIterations = 5;
  in->shouldRunPostZoltan = true;
  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] running MeshAdapt...\n";
  ma::adapt(in);
  reorderMdsMesh(mesh);

  sprintf(filename,"after-adapt-%d",adapt_time);
  apf::writeVtkFiles(filename,mesh);
  ++adapt_time;

  m3dc1_mesh::instance()->initialize();
  generate_global_numbering(mesh);

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 0);
    if(!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    m3dc1_field_synchronize(field);
#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
  destroyField(sizeField);
}

static double square(double x) {return x * x;}

int m3dc1_field_printcompnorm(FieldID* /* in */ field_id, const char* info)
{
  /* assumes DOFs are real */
  double* array=NULL;
  m3dc1_field_getdataptr(field_id, &array);
  int dof_in_array;
  m3dc1_field_getnumlocaldof(field_id, &dof_in_array);

  m3dc1_field* mf = NULL;  
  if (!m3dc1_ghost::instance()->num_local_ent[0])
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    

  apf::Field* f = mf ->get_field();
  int dof_per_node = countComponents(f);
  int dof_per_value = C1TRIDOFNODE;
  int values_per_node = dof_per_node / dof_per_value;
  vector<double> norms(values_per_node);
  int values_in_array = dof_in_array / dof_per_value;
  int nnodes = dof_in_array / dof_per_node;
  double* p = array;
  for (int i = 0; i < nnodes; ++i)
    for (int j = 0; j < values_per_node; ++j)
      for (int k = 0; k < dof_per_value; ++k)
        norms.at(j) += square(*p++);
  assert(nnodes * values_per_node * dof_per_value == dof_in_array);
  vector<double> buff = norms;
  MPI_Allreduce(&buff[0], &norms[0], values_per_node, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int psize;
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  if(PCU_Comm_Self() == psize-1)
  {
    std::cout<< "norm of vec "<<info;
    for(int i = 0; i < values_per_node; ++i)
      std::cout<<" "<<std::sqrt(norms[i]);
    std::cout<<std::endl;
  }
}

int sum_edge_data (double * data, int* size)
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;

  int num_edge=0, edg_dim=1;
  m3dc1_mesh_getnument(&edg_dim, &num_edge);
  PCU_Comm_Begin();
  for(int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
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
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]+=receive_buff[i];
    }

  PCU_Comm_Begin();
  for(int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    if (!is_ent_owned(m,e) || !m->isShared(e))
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
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]=receive_buff[i];
    }
   delete []receive_buff;
}

int get_node_error_from_elm (double * elm_data, int * size, double* nod_data)
{
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;    

  int num_node=0, num_elm=0, nod_dim=0, elm_dim=2;
  m3dc1_mesh_getnument(&nod_dim, &num_node);
  m3dc1_mesh_getnument(&elm_dim, &num_elm);
  PCU_Comm_Begin();
  double* buff = new double[*size];
  double* area = new double[num_node];
  for(int i=0; i<num_node; i++)
    area[i]=0.;
  for(int i=0; i<(*size)*num_elm; i++ )
    for(int j=0; j<*size; j++)
    {
      assert(elm_data[(*size)*i+j]==elm_data[(*size)*i+j]);
      assert(elm_data[(*size)*i+j]>=0);
    }
  for(int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
    apf::Adjacent adjacent;
    m->getAdjacent(e,2,adjacent);
    for(int j=0; j<adjacent.getSize(); j++)
    {
       apf::MeshElement* me = createMeshElement(m, adjacent[j]);
       double s = apf::measure(me);
       int ielm = getMdsIndex(m, adjacent[j]);
       assert(ielm>=0 &&ielm<num_elm);
       for(int k=0; k<*size; k++)
          nod_data[i*(*size)+k]+=s*elm_data[(*size)*ielm+k];
       area[i]+=s;
       destroyMeshElement(me);
    }

    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_COMM_PACK(own_partid, area[i]);
    PCU_Comm_Pack(own_partid, &(nod_data[(*size)*i]), sizeof(double)*(*size));
  }

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      int inode = getMdsIndex(m, node);
      double s;
      PCU_COMM_UNPACK(s);
      area[inode]+=s;
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]+=buff[i];
    }

  for(int i=0; i<num_node; i++)
  {
    for(int j=0; j<*size; j++)
      nod_data[i*(*size)+j]/=area[i];
  }
  PCU_Comm_Begin();
  for(int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    if (!is_ent_owned(m,e) || !m->isShared(e))
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
      int inode = getMdsIndex(m, node);
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]=buff[i];
    }
  delete []buff;
  delete []area;  
}
int m3dc1_field_max (FieldID* field_id, double * max_val, double * min_val)
{
  m3dc1_field* mf = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  else
    mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];    

  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if(mf->get_value_type()) dofPerEnt *= 2;
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  std::vector<double> maxVal(dofPerEnt, -1e30), minVal(dofPerEnt,1e30), dofs(dofPerEnt);
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for(int i=0; i<dofPerEnt; i++)
    {
      if(maxVal[i]<dofs[i]) maxVal[i]=dofs[i];
      if(minVal[i]>dofs[i]) minVal[i]=dofs[i];
    }
  }
  MPI_Allreduce(&(maxVal[0]), max_val, dofPerEnt, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(minVal[0]), min_val, dofPerEnt, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return M3DC1_SUCCESS;
}

#ifdef M3DC1_TRILINOS
#include <Epetra_MultiVector.h>
#include <AztecOO.h>
#include <Epetra_Version.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>

void verifyFieldEpetraVector(apf::Field* f, Epetra_MultiVector* x)
{
  double* field_data =getArrayData(f);
  assert(countComponents(f)*m3dc1_mesh::instance()->num_local_ent[0]==x->MyLength());

  for(int i=0; i<x->MyLength(); ++i)
  { 
    assert(!value_is_nan((*x)[0][i]) && !value_is_nan(field_data[i]));

    if (!(m3dc1_double_isequal((*x)[0][i], field_data[i])))
      std::cout<<"[p"<<PCU_Comm_Self()<<"] x["<<i<<"]="<<(*x)[0][i]
                <<", field_data["<<i<<"]="<<field_data[i]<<"\n";
      assert(m3dc1_double_isequal((*x)[0][i], field_data[i]));
  }
}
#endif

int m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }
  m3dc1_ls::instance()->add_matrix(*matrix_id, new m3dc1_epetra(*matrix_id, *matrix_type, *scalar_type, *field_id));
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_delete(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  typedef std::map<int, m3dc1_epetra*> matrix_container_map;
  m3dc1_ls::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  mat->destroy();
  delete mat;
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_insert(int* matrix_id, int* row, int* col, int* scalar_type, double* val)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  assert(*scalar_type==M3DC1_REAL);

  int err = mat->epetra_mat->ReplaceGlobalValues(*row, 1, val, col);
  if (err) {
    err =mat->epetra_mat->InsertGlobalValues(*row, 1, val, col);
    assert(err == 0);
  }
  return M3DC1_SUCCESS;
#endif
}

void print_elem (int elem_id)
{
  int ielm_dim;
  apf::Mesh2* m=NULL;
  
  int num_node_per_element;
  apf::Downward downward;
  if (!m3dc1_ghost::instance()->num_local_ent[0])  
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;

  ielm_dim = (m->getDimension()==2)? 2:3; 
  apf::MeshEntity* e = getMdsEntity(m, ielm_dim, elem_id);
  num_node_per_element = m->getDownward(e, 0, downward);
 
  int *id = new int[num_node_per_element];
  for (int i=0; i<num_node_per_element; ++i)
    id[i] = get_ent_globalid(m, downward[i], 0);

  switch (num_node_per_element)
  { 
    case 3: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<"\n";
            break;
    case 4: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<"\n";
            break;
    case 5: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<" "<<id[4]<<"\n";
            break;
    case 6: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
                      <<id[0]<<" "<<id[1]<<" "<<id[2]<<" "<<id[3]<<" "<<id[4]<<" "<<id[5]<<"\n";
            break;
    default: break;
  }
  delete [] id;
}

#ifdef M3DC1_TRILINOS
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  double val[1];
  int col[1];
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    for(int j=0; j<csize; j++)
    {
      col[0] = columns[j];
      val[0] = values[i*csize+j];
      int ierr = mat->SumIntoGlobalValues(rows[i], 1, val, col);
      if (ierr) 
        ierr =mat->InsertGlobalValues(rows[i], 1, val, col);
      assert(!ierr);
    }
  } // for i
}

// seol -- this does weird thing so shouldn't be used
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values_wrong(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    int ierr = mat->SumIntoGlobalValues(rows[i], csize, &values[i*csize], columns);
    if (ierr) 
      ierr =mat->InsertGlobalValues(rows[i], csize, &values[i*csize], columns);
    assert(!ierr);
  } // for i
}

#endif


int m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (!mat)
    return M3DC1_FAILURE;

  int field = mat->get_field_id();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = 2;
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  if (m3dc1_mesh::instance()->mesh->getDimension()==3) ielm_dim =3;
  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;
  int start_global_dof,end_global_dof_id;
  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowVarIdx<numVar && *columnVarIdx<numVar);
  int* rows = new int[dofPerVar*nodes_per_element];
  int* columns = new int[dofPerVar*nodes_per_element];

  if (mat->matrix_type==M3DC1_MULTIPLY)
  {
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    //FIXME: mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
    epetra_add_values(mat->epetra_mat, dofPerVar*nodes_per_element, 
                      rows,dofPerVar*nodes_per_element, columns, values);     
  }
  else //M3DC1_SOLVE
  {
    int nodeOwner[6];
    int columns_bloc[6], rows_bloc[6];
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowVarIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnVarIdx;
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      // FIXME: smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      epetra_add_values(mat->epetra_mat, dofPerVar, rows+dofPerVar*inode, 
                       dofPerVar*nodes_per_element, columns, values+offset);
      offset += numValuesNode;
    }
  }
  delete [] rows;
  delete [] columns;

  return M3DC1_SUCCESS;
#endif
}


int m3dc1_epetra_setbc(int* matrix_id, int* row)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  global_ordinal_type col[1]; col[0] = row_g;
  double val[1]; val[0]=1.0; 
 
  // MatSetValue(*A, row, row, 1.0, ADD_VALUES);

  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, 1, val, col);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, 1, val, col);
  assert(err == 0);

/* zero out non-diagonal columns
  apf::MeshEntity* e = apf::getMdsEntity(m3dc1_mesh::instance()->mesh, 0, inode);

  int MaxNumIndices = mat->epetra_mat->MaxNumEntries();
  int * Indices_int = 0;
  assert (mat->epetra_mat->RowMap().GlobalIndicesInt());
  Indices_int = new int[MaxNumIndices];
  double * values  = new double[MaxNumIndices];
  int err, NumIndices;
  mat->epetra_mat->ExtractGlobalRowCopy(row_g, MaxNumIndices, NumIndices, values, Indices_int);
  for (int j = 0; j < NumIndices ; j++)
  {
    col[0] = Indices_int[j];
    if (Indices_int[j]!=row_g) // zero out
    {
      val[0]=0.0; 
      err = mat->epetra_mat->ReplaceGlobalValues(row_g, 1, val, col);
      if (err) 
        err =mat->epetra_mat->InsertGlobalValues(row_g, 1, val, col);
      std::cout<<"("<<PCU_Comm_Self()<<") "<<__func__<<": "<<*matrix_id<<" - node "
           <<getNumber(get_global_numbering(), e, 0, 0)<<", ReplaceGlobalValues("<< row_g<<','<<col[0]<<",0.0"<<endl;
    } 
    else
    {
      val[0]=1.0; 
      int err = mat->epetra_mat->SumIntoGlobalValues(row_g, 1, val, col);
      if (err) 
        err =mat->epetra_mat->InsertGlobalValues(row_g, 1, val, col);
      std::cout<<"("<<PCU_Comm_Self()<<") "<<__func__<<": "<<*matrix_id<<" - node "
           <<getNumber(get_global_numbering(), e, 0, 0)<<", SumIntoGlobalValues("<< row_g<<','<<col[0]<<",1.0"<<endl;
    }
    assert(err == 0);
  }

  std::cout<<"("<<PCU_Comm_Self()<<") "<<__func__<<": "<<*matrix_id<<" - node "
           <<getNumber(get_global_numbering(), e, 0, 0)<<", local dof "<<*row<<", global dof "<<row_g<<std::endl;
*/
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  std::vector <global_ordinal_type> columns_g(*numVals);
  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  for(int i=0; i<*numVals; i++)
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;
//  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, *numVals, values, &columns_g[0]);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, *numVals, values, &columns_g[0]);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_print(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),/*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();
  A.Print(cout);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_write(int* matrix_id, const char* filename, int* skip_zero, int* start_index)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (!filename)
    return m3dc1_epetra_print(matrix_id);

  char matrix_filename[256];
  sprintf(matrix_filename,"%s-%d",filename, PCU_Comm_Self());
  if (*skip_zero==0)
    write_matrix(mat->epetra_mat, matrix_filename,false,*start_index);
  else
    write_matrix(mat->epetra_mat, matrix_filename,true,*start_index);

  // assemble matrix
  Epetra_Export exporter(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

  sprintf(matrix_filename,"assembled-%s-%d",filename, PCU_Comm_Self());

  if (*skip_zero==0)
    write_matrix(&A, matrix_filename,false,*start_index);
  else
    write_matrix(&A, matrix_filename,true,*start_index);

  return M3DC1_SUCCESS;
#endif
}

#ifdef M3DC1_TRILINOS
void copyEpetraVec2Field(Epetra_MultiVector* x, apf::Field* f)
{
  int start_global_dofid, num_dof = countComponents(f);
  std::vector<double> dof_data(num_dof);
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;

  int index=0;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_owned(m,e)) continue;
    for(int i=0; i<num_dof; ++i)
      dof_data.at(i)=(*x)[0][index++];
    setComponents(f, e, 0, &(dof_data[0]));
  }
  m->end(it);
  assert(index == num_dof*m3dc1_mesh::instance()->num_own_ent[0]);
  m3dc1_field_synchronize(f);
}
#endif

int m3dc1_solver_aztec(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, int* num_iter, double* tolerance)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  else
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<* x_fieldid<<" (tol "<<*tolerance<<")\n";

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),
			 /*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

//  const char* mfilename = mfilename_str.c_str();
//  EpetraExt::RowMatrixToMatlabFile(mfilename, A);

  // copy field to vec  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  m3dc1_field_synchronize(b_field);
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1);
  b.Export(b_field_vec,exporter,Insert);

  // vector for solution
  Epetra_MultiVector x(*(mat->_owned_map), 1);

  Epetra_LinearProblem problem(&A,&x,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_output,1);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  int overlap = 1;
  solver.SetAztecOption(AZ_overlap, overlap);

  solver.Iterate(*num_iter,*tolerance);
  mat->num_solver_iter = solver.NumIters();
  
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  copyEpetraVec2Field(&x, x_field);
  return M3DC1_SUCCESS;
#endif
}

// solve Ax=b
#ifdef M3DC1_TRILINOS
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
//#include "Amesos2.hpp"
//#include "Amesos2_Version.hpp"
#endif

int m3dc1_solver_amesos(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, const char* solver_name)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: check if amesos2 is available\n";
  return M3DC1_FAILURE;

  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  Epetra_Export exporter(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  // copy field to vec
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1 );

  A.Export(*(mat->epetra_mat),exporter,Add);
  b.Export(b_field_vec,exporter,Insert);

  A.FillComplete();
  A.OptimizeStorage();
  A.MakeDataContiguous();

// using SuperLUDIST

  // Before we do anything, check that the solver is enabled
/*
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if( !Amesos2::query(solver_name) )
  {
    if (!PCU_Comm_Self()) 
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": "<< solver_name << " not enabled.  Exiting..." << std::endl;
    return M3DC1_FAILURE;	// Otherwise CTest will pick it up as failure, which it isn't really
  }
  else
    if (!PCU_Comm_Self())  *fos <<__func__<<": "<<Amesos2::version() <<" with "<<solver_name<< std::endl;
  // Constructor from Factory
  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;

  Teuchos::RCP<MAT> rcp_A = Teuchos::rcp(&A);
//  A->Print(*fos);
//    *fos <<"["<<PCU_Comm_Self()<<"] GlobalNumNonZero="<<rcp_A->NumGlobalNonzeros()<<"\n";
//  int nrows = A->NumGlobalEntries();
//  if (!PCU_Comm_Self()) 
//      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": nrows="<<nrows<<std::endl;

    int numVecs=1;
    Epetra_MultiVector x(*(mat->_owned_map), 1 );
    Teuchos::RCP<MV> X = Teuchos::rcp(&x);
    Teuchos::RCP<MV> B = Teuchos::rcp(&b);

    // copy field to vec
    for (int i=0; i<B->MyLength(); ++i)
      (*B)[0][i] = b_field_data[i];

    // Solve A*Xhat = B for Xhat using the Superlu solver
    Teuchos::RCP<Amesos2::Solver<MAT,MV> > solver;
    try 
    {
      solver = Amesos2::create<MAT,MV>(solver_name, rcp_A, X, B );
    }
    catch (std::invalid_argument e)
    {
      if (!PCU_Comm_Self()) *fos <<"[M3D-C1 ERROR] "<<__func__<<": "<< e.what() << std::endl;
      return M3DC1_FAILURE;
    }

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    //solver->printTiming(*fos);
    //X.Describe(*fos, Teuchos::VERB_EXTREME);

  //  X->Print(*fos); //, Teuchos::VERB_EXTREME);
  // print residual
//  double norm_data[1];
//  x.Norm2(norm_data);
//  *fos << "["<<PCU_Comm_Self()<<"] Norm2 of Ax - b = " << norm_data[0] << std::endl;

  // get solution
  Epetra_Import importer(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_MultiVector sol_x(*(mat->_overlap_map),1);
  sol_x.Import(x, importer, Add);
 
  double** s;
  sol_x.ExtractView(&s);
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  double* x_field_data = getArrayData(x_field);
  for (int i=0; i<sol_x.MyLength(); ++i)
    x_field_data[i] =s[0][i];
*/
  return M3DC1_SUCCESS;

#endif
}

int m3dc1_solver_getnumiter(int* matrix_id, int * num_iter)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  *num_iter = mat->num_solver_iter;
  return M3DC1_SUCCESS;
#endif

}

// local matrix multiplication
// do accumulate(out_field) for global result
int m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_MULTIPLY)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  if (!mat->epetra_mat->Filled())
    mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());

  apf::Field* in_field = (*(m3dc1_mesh::instance()->field_container))[*in_fieldid]->get_field();
  double* in_field_data =getArrayData(in_field);

  Epetra_MultiVector x(mat->epetra_mat->RowMap(), 1);
  // copy field to vec
  for (int i=0; i<x.MyLength(); ++i)
    x[0][i] = in_field_data[i];
  Epetra_MultiVector b(mat->epetra_mat->RowMap(), 1);
  EPETRA_CHK_ERR(mat->epetra_mat->Multiply(false, x, b));
  apf::Field* out_field = (*(m3dc1_mesh::instance()->field_container))[*out_fieldid]->get_field();
  double* out_field_data =getArrayData(out_field);
  b.ExtractCopy(out_field_data, b.MyLength());
  accumulate(out_field);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_freeze(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());
  return M3DC1_SUCCESS;
#endif
}


// Ghosting Functions

#ifdef M3DC1_OMEGA_H
#include "apfOmega_h.h"
#endif

//*******************************************************
int m3dc1_ghost_load(int* nlayers)
//*******************************************************
{ 
#ifdef M3DC1_OMEGA_H
  if (m3dc1_model::instance()->local_planeid == 0)
  {
    assert(*nlayers>0);
    m3dc1_ghost::instance()->mesh = apf::makeEmptyMdsMesh(m3dc1_model::instance()->model, 2, false);
    apf::disownMdsModel(m3dc1_ghost::instance()->mesh);
    m3dc1_ghost::instance()->nlayers = *nlayers;
    
    // Set up ghosted mesh via omega_h
    osh_t osh_mesh = osh::fromAPF(m3dc1_mesh::instance()->mesh);
    osh_ghost(osh_mesh, m3dc1_ghost::instance()->nlayers);
    osh::toAPF(osh_mesh, m3dc1_ghost::instance()->mesh);
    osh_free(osh_mesh);   
  } else {
    assert(0);
  }
  m3dc1_ghost::instance()->initialize();
  // unlike nodes, global id for elements is not transferred so it should be generated explicitly
  generate_elem_global_numbering(m3dc1_ghost::instance()->mesh);

  // build orig_node_flag for "isghost" query
  m3dc1_ghost::instance()->org_node_flag = new vector<bool>;
  m3dc1_ghost::instance()->org_node_flag->resize(m3dc1_mesh::instance()->num_global_ent[0]);
  for (std::vector<bool>::iterator vit=m3dc1_ghost::instance()->org_node_flag->begin();
        vit!=m3dc1_ghost::instance()->org_node_flag->end(); ++vit)
    *vit = false;

  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))

    m3dc1_ghost::instance()->org_node_flag->at(get_ent_globalid(m3dc1_mesh::instance()->mesh, e, 0)) = true;
  m3dc1_mesh::instance()->mesh->end(it);

  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 INFO] "<<__func__<<" (nlayers="<<*nlayers<<")\n";
  return M3DC1_SUCCESS;
#else
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<": not supported\n";
  return M3DC1_FAILURE;
#endif
}

//*******************************************************
int m3dc1_ghost_delete()
//*******************************************************
{
#ifdef M3DC1_OMEGA_H
  // Destroy global numbering and associated field
  // Delete existing numbering
  destroy_node_global_numbering(m3dc1_ghost::instance()->mesh);
  destroy_elem_global_numbering(m3dc1_ghost::instance()->mesh);
  m3dc1_ghost::instance()->org_node_flag->clear();
  delete m3dc1_ghost::instance()->org_node_flag;
  m3dc1_ghost::destroy();
  return M3DC1_SUCCESS;
#else
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<": not supported\n";
  return M3DC1_FAILURE;
#endif
}

// Adjacency search

//*******************************************************
int m3dc1_mesh_search(int* initial_simplex,
		      double* final_position,
		      int* final_simplex)
//*******************************************************
{
  bool located = false;
  apf::MeshEntity* e = NULL;
  apf::MeshEntity* simplex = NULL;
  apf::Mesh2* m = NULL;
  if (!m3dc1_ghost::instance()->num_local_ent[0]) 
    m = m3dc1_mesh::instance()->mesh;
  else
    m = m3dc1_ghost::instance()->mesh;
  apf::Adjacent adjacent;
  int edge_curr_index, edge_prev_index;  
  int simplex_dim = m->getDimension();
  int vertex_dim = 0, edge_dim = 1;
  int bridge_dim = simplex_dim - 1;
  int count = 0, max_count = m->count(simplex_dim);
  double tol = 1e-15;
  
  simplex = apf::getMdsEntity(m, simplex_dim, *initial_simplex);
  while (not located) {
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
    bool b_negative[3] = {false, false, false};
    int bneg_index[2] = {0, 0};
    int bneg_count = 0;
    for (int j = 0; j < 3; ++j)
      if (b_coords[j] < 0) {
	b_negative[j] = true;
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
      for (size_t j = 0; j < ne; ++j) {
	int edge_tmp_index = apf::getMdsIndex(m, edges[j]);
	if ((edge_tmp_index != edge_curr_index) &&
	    (edge_tmp_index != edge_prev_index)) {
	  e = edges[j];
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
      return M3DC1_FAILURE;

    }
  }
  return M3DC1_SUCCESS;
}





