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
      dof_data = apf::getArrayData(new_fields[field_id*num_rank+from]);
      PCU_Comm_Unpack(&(dof_data[0]), total_ndof*recv_num_ent*sizeof(double));
      std::cout<<"get dof from p "<<from<<" into new_fields["<<field_id*num_rank+from<<"] \n";
    }
  }
}

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
  	 
  // m3dc1_field_delete (&field_id_h1);
  // m3dc1_field_delete (&field_id_h2);
  destroyField(f_h1);
  destroyField(f_h2);

  // remove f from field container
  delete (*m3dc1_mesh::instance()->field_container)[field_id_h1];
  delete (*m3dc1_mesh::instance()->field_container)[field_id_h2];
  m3dc1_mesh::instance()->field_container->erase(field_id_h1);
  m3dc1_mesh::instance()->field_container->erase(field_id_h2);


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
    if (isFrozen(field)) unfreeze(field);
    fields.push_back(field);
    it++;
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

  map<int, apf::Field*> new_fields; 

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

    int mypid = pumi_rank(), field_id=0, num_rank=pumi_size(), ndof_per_plane;
    int num_field=m3dc1_mesh::instance()->field_container->size();

    it=m3dc1_mesh::instance()->field_container->begin();
    while(it!=m3dc1_mesh::instance()->field_container->end())
    {
      apf::Field* field = it->second->get_field();

      // master plane
      apf::Field* new_field=NULL;
      double* data_array = NULL;
      int num_dof=countComponents(field);
      if (!m3dc1_model::instance()->local_planeid)
      {
        std::string f_name(getName(field));
        f_name += "_p";
        for (int i=1; i<m3dc1_model::instance()->num_plane; ++i)
        { 
          int pid=field_id*num_rank+mypid+i*m3dc1_model::instance()->group_size;
          stringstream ss;
          ss << pid;
          std::string f_name_2=f_name + ss.str();
          new_fields[pid] = createPackedField(m3dc1_mesh::instance()->mesh, f_name_2.c_str(), num_dof);
          freeze(new_fields[pid]);
          std::cout<<"field "<<field_id<<": creating new field ["
                   <<pid<<"] ID "<<f_name_2.c_str()<<"\n";
        }  
      }

      copy_field(mesh, field, new_fields, field_id);

      // let's keep the old one for now
      if (isFrozen(field)) unfreeze(field);
      it++;
      ++field_id;
    }

    // remove entities on non-master plane
    if (m3dc1_model::instance()->local_planeid)
    {
      for (int dim=2; dim>=0; --dim)
      {
        MeshIterator* ent_it = mesh->begin(dim);
        while ((e = mesh->iterate(ent_it)))
          mesh->destroy(e);
        mesh->end(ent_it);
      }
    }
  }
  
  ma::adapt(in);

  mesh->removeField(size_field);
  mesh->removeField(frame_field);
  apf::destroyField(size_field);
  apf::destroyField(frame_field);

  m3dc1_mesh::instance()->initialize(false);

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    // re-construct wedges
    // compute the fields to copy from master plane
    int num_field = 0;
    int* field_id = NULL;
    int* num_dofs_per_value = NULL;

    if (!m3dc1_model::instance()->local_planeid && m3dc1_mesh::instance()->field_container)
    {
      num_field = m3dc1_mesh::instance()->field_container->size();
      field_id = new int [num_field];
      num_dofs_per_value = new int [num_field];
      int i=0;
      it=m3dc1_mesh::instance()->field_container->begin();
      while(it!=m3dc1_mesh::instance()->field_container->end())
      {
        field_id[i] = it->second->get_id();
        num_dofs_per_value[i] = it->second->get_dof_per_value();
        ++i;
      }
    }

    m3dc1_mesh::instance()->build3d(num_field, field_id, num_dofs_per_value, 
                                    ge_tag, false, &new_fields);

    for(int i = 0; i < 3; ++i)
      delete [] ge_tag[i];
    delete [] ge_tag;

    if (num_field>0)
    {
      delete [] field_id;
      delete [] num_dofs_per_value;
    }
  } // 3d 

  reorderMdsMesh(mesh);

  // FIXME: crash in 3D 
  if (m3dc1_model::instance()->num_plane==1)
    apf::writeVtkFiles("after-adapt", mesh);
 
  m3dc1_mesh::instance()->initialize(false);

  compute_globalid(mesh, 0);
  compute_globalid(mesh, mesh->getDimension());

  if (m3dc1_mesh::instance()->field_container && 
      m3dc1_mesh::instance()->field_container->size())
  {
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
}
