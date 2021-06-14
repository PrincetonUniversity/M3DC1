/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_mesh.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <PCU.h>
#include "apfMDS.h"
#include "ReducedQuinticImplicit.h"
#include "m3dc1_slnTransfer.h"
#include "apfShape.h" // getLagrange
#include "pumi.h"
#include "ma.h"
#include "gmi_base.h"

using namespace apf;

int get_local_planeid(int p)
{
  return p/m3dc1_model::instance()->group_size;
}

void update_remotes(apf::Mesh2* mesh, MeshEntity* e)
{
  int myrank=PCU_Comm_Self();

  apf::Copies remotes;
  apf::Copies new_remotes;
  apf::Parts parts;
  mesh->getRemotes(e,remotes);

  APF_ITERATE(apf::Copies,remotes,rit)
  {
    int p=rit->first;
    if (get_local_planeid(p)==m3dc1_model::instance()->local_planeid) // the same local plane
    {
      new_remotes[p]=rit->second; 
      parts.insert(p);
    }
  }
  parts.insert(myrank);
  mesh->clearRemotes(e);
  mesh->setRemotes(e,new_remotes);
  mesh->setResidence(e, parts);  // set pclassification
}

void m3dc1_mesh::remove_wedges()
{       
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshEntity* e;
  int myrank = PCU_Comm_Self();
  int num_2d_face=mesh->count(3);
 
  // delete regions
  MeshIterator* ent_it = mesh->begin(3);
  while ((e = mesh->iterate(ent_it)))
    mesh->destroy(e);
  mesh->end(ent_it);
  m3dc1_mesh::instance()->num_own_ent[3]=m3dc1_mesh::instance()->num_local_ent[3]=0;

  // delete faces
  int cnt=0;
  ent_it = mesh->begin(2);
  while ((e = mesh->iterate(ent_it)))
  {
    // original entity for 2D plane
    if (cnt<num_2d_face)
      update_remotes(mesh, e);
    else
       mesh->destroy(e);
    ++cnt;
  }
  mesh->end(ent_it);
  m3dc1_mesh::instance()->num_local_ent[2] = mesh->count(2);

  // remote edges and vertices
  for (int dim=1; dim>=0; --dim)
  {
    ent_it = mesh->begin(dim);
    while ((e = mesh->iterate(ent_it)))
    {
      apf::Adjacent adjacent;
      mesh->getAdjacent(e,dim+1,adjacent);
      if (adjacent.getSize()) // original entity for 2D plane
        update_remotes(mesh, e);
      else
        mesh->destroy(e);
    }
    mesh->end(ent_it);
    m3dc1_mesh::instance()->num_local_ent[dim]=mesh->count(dim);
  }

  mesh->acceptChanges(); // update partition model
  changeMdsDimension(mesh, 2);
  set_mcount();

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Wedges between 2D planes removed ***\n";
}

void m3dc1_mesh::create_wedges()
{

}
void compute_size_and_frame_fields(apf::Mesh2* m, double* size_1, double* size_2, 
     double* angle, apf::Field* sizefield, apf::Field* framefield)
{
  MeshEntity* e;
  MeshIterator* ent_it = m->begin(0);
  int i=0;
  while ((e = m->iterate(ent_it)))
  {
    double h1 = size_1[i];
    double h2 = size_2[i];

    double angle_1[3];
    angle_1[0] = angle[(i*3)];
    angle_1[1] = angle[(i*3)+1];
    angle_1[2] = angle[(i*3)+2];

 // Calculate the second unit vector
/*  double a, b;
    double frac_1, frac_2;
    frac_1 = (angle_1[0])*(angle_1[0]);
    frac_2 = (angle_1[0])*(angle_1[0]) + (angle_1[1])*(angle_1[1]);


    b = sqrt (frac_1/frac_2);
    a = -(angle_1[1]*b)/angle_1[0];

    double mag = sqrt (a*a + b*b);
    double dir_2[3];
    dir_2[0] = a /mag;
    dir_2[1] = b /mag;
    dir_2[2] = 0.0;
*/
    double dir_2[3];
    dir_2[0] = -angle_1[1];
    dir_2[1] =  angle_1[0];
    dir_2[2] =  0.0;

    ma::Vector h(h1, h2, h2);

    ma::Matrix r;
    r[0][0]=angle_1[0];
    r[0][1]=angle_1[1];
    r[0][2]=0.0;

    r[1][0]= dir_2[0];
    r[1][1]= dir_2[1];
    r[1][2]=0.0;

    r[2][0]=0;
    r[2][1]=0;
    r[2][2]=1.;

    apf::setVector(sizefield, e, 0, h);
    apf::setMatrix(framefield, e, 0, r);
//  apf::setMatrix(framefield, e, 0, apf::transpose(r));	// For Shock Test Case
  }
  m->end(ent_it);

  // sync the fields to make sure verts on part boundaries end up with the same size and frame
  apf::synchronize(sizefield);
  apf::synchronize(framefield);
}

// *********************************************************
void copy_field (apf::Mesh2* mesh, apf::Field* f, double* master_data_array)
// *********************************************************
{
  int num_2d_vtx=mesh->count(0);

  // copy the existing dof data
  int total_ndof = apf::countComponents(f);

  double* dof_val = new double[total_ndof*num_2d_vtx];

  freeze(f); // switch dof data from tag to array

  PCU_Comm_Begin();

  int proc=PCU_Comm_Self()%m3dc1_model::instance()->group_size;
  if (m3dc1_model::instance()->local_planeid)
  {
    memcpy(&(dof_val[0]), apf::getArrayData(f), total_ndof*num_2d_vtx*sizeof(double));
    PCU_Comm_Pack(proc, &num_2d_vtx, sizeof(int));
    PCU_Comm_Pack(proc, &(dof_val[0]), total_ndof*num_2d_vtx*sizeof(double));
  }
  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_val;

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      int index=total_ndof*recv_num_ent*from;
      PCU_Comm_Unpack(&(master_data_array[index]), total_ndof*recv_num_ent*sizeof(double));
    }
  }
  delete [] dof_val;
}


// *********************************************************
void copy_field (apf::Mesh2* mesh, apf::Field* f, 
                 map<int, apf::Field*>& new_fields, int field_id)
// *********************************************************
{
  int num_2d_vtx=mesh->count(0);
  int num_rank = pumi_size();
  // copy the existing dof data
  int total_ndof = apf::countComponents(f);

  double* dof_data = new double[total_ndof*num_2d_vtx];

  freeze(f); // switch dof data from tag to array

  PCU_Comm_Begin();

  int proc=PCU_Comm_Self()%m3dc1_model::instance()->group_size;
  if (m3dc1_model::instance()->local_planeid)
  {
    memcpy(&(dof_data[0]), apf::getArrayData(f), total_ndof*num_2d_vtx*sizeof(double));
    PCU_Comm_Pack(proc, &num_2d_vtx, sizeof(int));
    PCU_Comm_Pack(proc, &(dof_data[0]), total_ndof*num_2d_vtx*sizeof(double));
  }
  PCU_Comm_Send();
  int recv_num_ent;
  double* recv_dof_data = new double[total_ndof*num_2d_vtx];

  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&recv_num_ent, sizeof(int));
      if (!new_fields[field_id*num_rank+from])
        std::cout<<"p"<<PCU_Comm_Self()<<"] cannot find new_fields["<<field_id*num_rank+from<<"\n";
      dof_data = apf::getArrayData(new_fields[field_id*num_rank+from]);
      PCU_Comm_Unpack(&(dof_data[0]), total_ndof*recv_num_ent*sizeof(double));
#ifdef DEBUG
      std::cout<<"[p"<<PCU_Comm_Self()<<"] saving field "
               <<field_id<<"@p"<<from<<" into new_fields["<<field_id*num_rank+from<<"] \n";
#endif
    }
  }
}

// run mesh adaptation in multiple planes in 3D
void adapt_mesh (int field_id_h1, int field_id_h2, double* dir,
    int shouldSnap, int shouldRunPreZoltan ,int shouldRunPostZoltan,
    int shouldRefineLayer, int maximumIterations, double goodQuality)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;

  apf::Field* f_h1 = (*m3dc1_mesh::instance()->field_container)[field_id_h1]->get_field();
  synchronize_field(f_h1);
  int num_dof = countComponents(f_h1);
  if (!isFrozen(f_h1)) freeze(f_h1);
  double* data_h1;

  if (num_dof==1)
    data_h1 = apf::getArrayData(f_h1);
  else 
  {
    data_h1 = new double[mesh->count(0)];
    double* temp_data = apf::getArrayData(f_h1);
    for (int i=0; i<mesh->count(0); ++i)
      data_h1[i] = temp_data[i*num_dof];
  }

  apf::Field* f_h2 = (*m3dc1_mesh::instance()->field_container)[field_id_h2]->get_field();
  synchronize_field(f_h2);
  if (!isFrozen(f_h2)) freeze(f_h2);

  double* data_h2;
  if (num_dof==1)
    data_h2 = apf::getArrayData(f_h2);
  else
  {
    data_h2= new double[mesh->count(0)];
    double* temp_data = apf::getArrayData(f_h2);
    for (int i=0; i<mesh->count(0); ++i)
      data_h2[i] = temp_data[i*num_dof];
  }

  apf::Field* size_field = apf::createField(mesh, "size_field", apf::VECTOR, apf::getLagrange(1));
  apf::Field* frame_field = apf::createField(mesh, "frame_field", apf::MATRIX, apf::getLagrange(1));

  compute_size_and_frame_fields(mesh, data_h1, data_h2, dir, size_field, frame_field);
  	 
   m3dc1_field_delete (&field_id_h1);
   m3dc1_field_delete (&field_id_h2);

  if (num_dof>1)
  {
    delete [] data_h1;
    delete [] data_h2;
  }

  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*>::iterator it=m3dc1_mesh::instance()->field_container->begin();

  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 1);
    fields.push_back(field);
    it++;
  }
  
  apf::unfreezeFields(mesh); // turning field data from array to tag

  int num_field = fields.size();
  int* field_id = NULL;
  int* num_dofs_per_value = NULL;

  if (fields.size()>0)
  {
    num_field = m3dc1_mesh::instance()->field_container->size();
    field_id = new int [num_field];
    num_dofs_per_value = new int [num_field];
    int i=0;
    it=m3dc1_mesh::instance()->field_container->begin();
    for(;it!=m3dc1_mesh::instance()->field_container->end();++it)
    {
      field_id[i] = it->second->get_id();
      num_dofs_per_value[i] = it->second->get_dof_per_value();
    }
  }

  if (!PCU_Comm_Self())
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": "<<fields.size() <<" fields will be transfered\n";

  // delete all the matrix
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    apf::destroyNumbering(n);
  }
	
  ReducedQuinticImplicit shape;
  ReducedQuinticTransfer slnTransfer(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh, size_field, frame_field, &slnTransfer);
	
  in->shouldSnap = 0; // FIXME: crash if *shouldSnap==1;
  in->shouldTransferParametric = 0;
  in->shouldRunPreZoltan = shouldRunPreZoltan;
  in->shouldRunPostZoltan = shouldRunPostZoltan;
  in->shouldRunMidParma = 0;
  in->shouldRunPostParma = 0;
  in->shouldRefineLayer = shouldRefineLayer;
  in->maximumIterations=maximumIterations;
  in->goodQuality = goodQuality;

  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": snap "<<shouldSnap
  	  <<", runPreZoltan "<<shouldRunPreZoltan<<", runPostZoltan "<<shouldRunPostZoltan<<"\n";

#ifdef DEBUG
  apf::writeVtkFiles("before-adapt", mesh);
  mesh->writeNative("mesh.smb");
#endif

  apf::MeshEntity* e;

  gmi_model* model = m3dc1_model::instance()->model;
  int** ge_tag = new int*[3];
  gmi_iter* g_it;
  gmi_ent* ge;
  int cnt;

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    // save the current model tag and revert it back to the original
    for (int i = 0; i < 3; ++i)
      ge_tag[i] = new int[model->n[i]];

    for (int dim=0; dim<=2; ++dim)
    {
      g_it = gmi_begin(model, dim);
      cnt=0;
      while ((ge = gmi_next(model, g_it)))
      {
        ge_tag[dim][cnt++]= gmi_tag(model,ge);
        gmi_base_set_tag (model, ge, cnt);
      }
      gmi_end(model, g_it);
    }

    // remove 3D entities
    m3dc1_mesh::instance()->remove_wedges();
    apf::printStats(mesh);

    apf::writeVtkFiles("after-wedge-removal", mesh);
    mesh->writeNative("after-wedge-removal.smb");
  }

  ma::adapt(in);

  mesh->removeField(size_field);
  mesh->removeField(frame_field);
  apf::destroyField(size_field);
  apf::destroyField(frame_field);

  m3dc1_mesh::instance()->set_mcount();

  std::cout<<"After adaptation: ";
  m3dc1_mesh::instance()->print(__LINE__);

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    // re-construct wedges
    m3dc1_mesh::instance()->create_wedges();

    for(int i = 0; i < 3; ++i)
      delete [] ge_tag[i];
    delete [] ge_tag;
  } // 3d 

  if (num_field>0)
  {
    delete [] field_id;
    delete [] num_dofs_per_value;
  }

  reorderMdsMesh(mesh);

  // FIXME: crash in 3D 
  apf::writeVtkFiles("after-adapt", mesh);

  m3dc1_mesh::instance()->set_mcount();
  compute_globalid(mesh, 0);
  compute_globalid(mesh, mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    if (!isFrozen(field)) freeze(field);
    synchronize_field(field);
    it++;
  }
}
