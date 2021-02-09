/******************************************************************************

  (c) 2005-2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "m3dc1_mesh.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include <PCU.h>
#include "apfMDS.h"
#include "ReducedQuinticImplicit.h"
#include "m3dc1_slnTransfer.h"
#include "apfShape.h" // getLagrange
#include "pumi.h"
#include "ma.h"

using namespace apf;

void compute_size_and_frame_fields(apf::Mesh2* m, double* size_1, double* size_2, 
     double* angle, apf::Field* sizefield, apf::Field* framefield)
{
  for (int i = 0; i<m->count(0); ++i)
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

    apf::MeshEntity* vert = getMdsEntity(m, 0, i);
    apf::setVector(sizefield, vert, 0, h);
    apf::setMatrix(framefield, vert, 0, r);
//  apf::setMatrix(framefield, vert, 0, apf::transpose(r));	// For Shock Test Case
  }
  // sync the fields to make sure verts on part boundaries end up with the same size and frame
  apf::synchronize(sizefield);
  apf::synchronize(framefield);
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

  // delete all the matrix
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
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

  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
//    if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": numbering "<<getName(n)<<" deleted\n";
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

  MPI_Comm groupComm;
  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    m3dc1_mesh::instance()->remove_wedges();
    pumi_mesh_print(mesh);
    int group_size = PCU_Comm_Peers()/m3dc1_model::instance()->num_plane;
    int groupRank = PCU_Comm_Self()%group_size; // modulo

    MPI_Comm_split(m3dc1_model::instance()->oldComm, m3dc1_model::instance()->local_planeid, groupRank, &groupComm);
    PCU_Switch_Comm(groupComm);
  }

  if (m3dc1_model::instance()->local_planeid == 0) // master plane
    ma::adapt(in);

  mesh->removeField(size_field);
  mesh->removeField(frame_field);
  apf::destroyField(size_field);
  apf::destroyField(frame_field);

  // remove non-master plane
  if (m3dc1_model::instance()->local_planeid)
  {
    apf::MeshEntity* e;
    for (int d=2; d>=0; --d)
    {
      apf::MeshIterator* ent_it = mesh->begin(d);
      while ((e = mesh->iterate(ent_it)))
        mesh->destroy(e);
      mesh->end(ent_it);
    }
  }

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
    MPI_Comm_free(&groupComm);

    // re-construct wedges
   std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<": "<<__LINE__<<"\n";

    // compute the fields to copy from master plane
    int num_field = 0;
    int* field_id = NULL;
    int* num_dofs_per_value = NULL;

    if (!m3dc1_model::instance()->local_planeid && m3dc1_mesh::instance()->field_container)
    {
      num_field = m3dc1_mesh::instance()->field_container->size();
      field_id = new int [num_field];
      num_dofs_per_value = new int [num_field];
      for (int i=0; i<num_field; ++i)
      {
        m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
        field_id[i] = i;
        num_dofs_per_value[i] = mf->get_dof_per_value();
      }
    }

    m3dc1_mesh::instance()->build3d(num_field, field_id, num_dofs_per_value, false);
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

  m3dc1_mesh::instance()->initialize();
  compute_globalid(mesh, 0);
  compute_globalid(mesh, mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    if (!isFrozen(field)) freeze(field);

#ifdef DEBUG
    int isnan;
    int fieldId= it->first;
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    synchronize_field(field);

#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
}
