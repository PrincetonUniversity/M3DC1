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
#include <pcu_util.h>
#include "apfMDS.h"
#include "apfField.h"
#include "ReducedQuinticImplicit.h"
#include "m3dc1_slnTransfer.h"
#include "apfShape.h" // getLagrange
#include "pumi.h"
#include "ma.h"
#include "gmi_base.h"
#include "apfMesh.h"
#include <assert.h>

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

// this will remove 3D entities and non-master planes
void m3dc1_mesh::remove3D()
{       
  if (!(m3dc1_model::instance()->num_plane)) // 2D
    return;

  m3dc1_model::instance()->save_gtag();

  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  MeshEntity* e;
  int cnt = 0;
  int myrank = PCU_Comm_Self();
  int num_2d_face=mesh->count(3);
 
  // delete regions
  MeshIterator* ent_it = mesh->begin(3);
  while ((e = mesh->iterate(ent_it)))
    mesh->destroy(e);
  mesh->end(ent_it);
  m3dc1_mesh::instance()->num_own_ent[3]=m3dc1_mesh::instance()->num_local_ent[3]=0;

  // delete faces
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

  if (m3dc1_model::instance()->local_planeid)
    for (int dim = 2; dim >= 0; dim--) 
    {
      ent_it = mesh->begin(dim);
      while ((e = mesh->iterate(ent_it)))
        mesh->destroy(e);
      mesh->end(ent_it);
    }

  mesh->acceptChanges();
  changeMdsDimension(mesh, 2);

  set_mcount();

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Wedges and non-master 2D planes removed ***\n";
}

// this static function only used in m3dc1_mesh::rebuildPointersOnNonMasterPlane
static int get_id_in_container(
    std::map<FieldID, m3dc1_field*>& fcontainer,
    apf::Field* f)
{
  typedef std::map<FieldID, m3dc1_field*> fct;

  int id = -1;
  for (fct::iterator it = fcontainer.begin(); it != fcontainer.end(); ++it)
  {
    m3dc1_field* mf = it->second;
    if (f == mf->get_field())
    {
      id = it->first;
      break;
    }
  }
  return id;
}
// this will rebuild the mesh data-structure on non-master plane after remove3D
//
// Notes
// (A) this is done to ensure "restore3D" behaves the same as "build3d".
// More specifically, the order of entity creation and iteration remains the same
// on all non-master planes.
//
// (B) the field pointers passed in the arrays will be updated for non-master planes
void m3dc1_mesh::rebuildPointersOnNonMasterPlane(
    std::vector<std::vector<apf::Field*> >& pFields, // multi-plane fields
    std::vector<apf::Field*>& zFields)             // fields that need to be zero-d
{
  if (!(m3dc1_model::instance()->num_plane)) // if 2D do nothing
    return;

  // store the names so they can be updated later on
  std::vector<std::vector<std::string> > pFieldNames; pFieldNames.clear();
  std::vector<std::string>              zFieldNames; zFieldNames.clear();
  std::vector<std::string>              tempNames;

  for (int i = 0; i < (int)pFields.size(); i++)
  {
    tempNames.clear();
    for (int j = 0; j < pFields[i].size(); j++)
    {
      tempNames.push_back(std::string(apf::getName(pFields[i][j])));
    }
    PCU_ALWAYS_ASSERT(pFields[i].size() == tempNames.size());
    pFieldNames.push_back(tempNames);
  }
  PCU_ALWAYS_ASSERT(pFieldNames.size() == pFields.size());

  for (int i = 0; i < (int)zFields.size(); i++)
    zFieldNames.push_back(std::string(apf::getName(zFields[i])));
  PCU_ALWAYS_ASSERT(zFieldNames.size() == zFields.size());


  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  std::vector<int> fncs; fncs.clear();
  std::vector<int> ids; ids.clear();  // if the field is in the m3dc1_mesh::instance()->field_container hold the id here
  std::vector<std::string> fnames; fnames.clear();
  std::vector<apf::FieldShape*> fshapes; fshapes.clear();

  for (int i = 0; i < mesh->countFields(); i++) {
    apf::Field* f = mesh->getField(i);
    fncs.push_back(f->countComponents());
    fnames.push_back(std::string(apf::getName(f)));
    fshapes.push_back(apf::getShape(f));
    int id = get_id_in_container(*m3dc1_mesh::instance()->field_container, f);
    ids.push_back(id);
  }

  PCU_Barrier();

  // now remove everything on non-master planes, and recreate them
  // Removal Phase
  // =============
  // a-fields
  // b-numberings
  // c-tags
  // d-internal mesh data-structure by calling destroyNative
  // =============
  // Recreate Phase
  // a-create empty mds meshes on non-master planes
  // b-set local_entid_tag, own_partid_tag, num_global_adj_node_tag, num_own_adj_node_tag to NULL
  if (m3dc1_model::instance()->local_planeid != 0)
  {
    // manually delete all the fields/numberings and associated tags
    while ( mesh->countFields() )
    {
      apf::Field* f = mesh->getField(0);
      mesh->removeField(f);
      apf::destroyField(f);
    }

    while ( mesh->countNumberings() )
    {
      apf::Numbering* n = mesh->getNumbering(0);
      mesh->removeNumbering(n);
      apf::destroyNumbering(n);
    }

    apf::DynamicArray<apf::MeshTag*> tags;
    mesh->getTags(tags);
    for (int i=0; i<tags.getSize(); i++)
    {
      if (mesh->findTag("norm_curv")==tags[i]) continue;
      for (int idim=0; idim<4; idim++)
        apf::removeTagFromDimension(mesh, tags[i], idim);
      mesh->destroyTag(tags[i]);
    }

    // destroy the native mesh
    apf::disownMdsModel(mesh);
    mesh->destroyNative();
    apf::destroyMesh(mesh);
    m3dc1_mesh::instance()->mesh = pumi_mesh_create(pumi::instance()->model, 2, false);
    mesh = m3dc1_mesh::instance()->mesh;

    local_entid_tag = NULL;
    own_partid_tag = NULL;
    num_global_adj_node_tag = NULL;
    num_own_adj_node_tag = NULL;
  }

  PCU_Barrier();


  // update the field pointers in m3dc1_mesh::field_container
  if (m3dc1_model::instance()->local_planeid != 0)
  {
    for (int i = 0; i < (int)fncs.size(); i++)
      {
	apf::Field* newf = apf::createPackedField(mesh, fnames[i].c_str(), fncs[i], fshapes[i]);
	apf::zeroField(newf);
	if (ids[i] > -1)
	{
	  m3dc1_field* mf = (*m3dc1_mesh::instance()->field_container)[ids[i]];
	  mf->set_field(newf);
	}
      }
  }

  // update the field pointers in pFields
  if (m3dc1_model::instance()->local_planeid != 0)
  {
    for (int i = 0; i < (int)pFieldNames.size(); i++)
    {
      for (int j = 0; j < pFieldNames[i].size(); j++)
      {
	apf::Field* f = mesh->findField(pFieldNames[i][j].c_str());
	PCU_ALWAYS_ASSERT(f);
	pFields[i][j] = f;
      }
    }
  }


  // update the field pointers in zFields
  if (m3dc1_model::instance()->local_planeid != 0)
  {
    for (int i = 0; i < (int)zFieldNames.size(); i++)
    {
      apf::Field* f = mesh->findField(zFieldNames[i].c_str());
      PCU_ALWAYS_ASSERT(f);
      zFields[i] = f;
    }
  }


  // since we have removed the Linear numbering for non-master planes, we add it beck here
  apf::Numbering* linnumbering = mesh->findNumbering("Linear");
  if (!linnumbering)
    apf::createNumbering(mesh, "Linear", mesh->getShape(), 1);

  mesh->acceptChanges();
  set_mcount();

  if (!PCU_Comm_Self())
    std::cout<<"\n*** Data structures have been re-initialized on non-master planes ***\n";
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
    double a, b;
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
    i++;    
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
void adapt_mesh (int field_id_h1, int field_id_h2, double* dir)
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
#ifdef OLDMA
  ma::Input* in = ma::configure(mesh, size_field, frame_field, &slnTransfer);
#else
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh, size_field, frame_field, &slnTransfer));
#endif
	
  in->shouldSnap = 0;
  in->shouldTransferParametric = 0;
  in->shouldRunMidZoltan = 1;
  in->shouldRunPreZoltan = 1;
  in->shouldRunPostZoltan = 1;
  in->maximumIterations = 5;
  in->userDefinedLayerTagName = "doNotAdapt";		// Works only when a tag "doNotAdapt" on elements is set
#ifdef _DEBUG
  if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] "<<__func__<<": runMidZoltan "<< in->shouldRunMidZoltan
  	  <<", runPreZoltan "<<in->shouldRunPreZoltan<<", runPostZoltan "<<in->shouldRunPostZoltan<<"\n"
          <<"should snap: " << in->shouldSnap
          <<", shouldTransferParametric: " << in->shouldTransferParametric << "\n";;

  static int ts=1;
  char fname[64];
  sprintf(fname, "before-adapt-ts%d", ts);
  apf::writeVtkFiles(fname, mesh);
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
    //m3dc1_mesh::instance()->remove_wedges(NULL);
    apf::printStats(mesh);

#ifdef DEBUG
    sprintf(fname, "after-wedge-removal-ts%d", ts);
    apf::writeVtkFiles(fname, mesh);
    sprintf(fname, "after-wedge-removal-ts%d.smb", ts);
    mesh->writeNative(fname);
#endif
  }

  mesh->acceptChanges();

  // do adaptation here
  ma::adaptVerbose(in);

  mesh->removeField(size_field);
  mesh->removeField(frame_field);
  apf::destroyField(size_field);
  apf::destroyField(frame_field);

  if (m3dc1_model::instance()->num_plane>1) // 3d
  {
    m3dc1_mesh::instance()->set_mcount();
    // re-construct wedges
    //m3dc1_mesh::instance()->create_wedges();

    for(int i = 0; i < 3; ++i)
      delete [] ge_tag[i];
    delete [] ge_tag;
  } // 3d 

  // FIXME: crash in 2D if no pre/post zoltan 
  reorderMdsMesh(mesh);
  m3dc1_mesh::instance()->initialize();

  compute_globalid(mesh, 0);
  compute_globalid(mesh, mesh->getDimension());

  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if (complexType) group_complex_dof(field, 0);
    synchronize_field(field);
    it++;
  }

  apf::freezeFields(mesh); // turning field data from tag to array
 
#ifdef DEBUG
  if (!PCU_Comm_Self()) std::cout<<"After adaptation: ";

  // FIXME: crash in 3D 
  sprintf(fname, "after-adapt-ts%d", ts++);
  apf::writeVtkFiles(fname, mesh);
/*
  MeshIterator* ent_it = mesh->begin(2);
  while ((e = mesh->iterate(ent_it)))
  {
    std::cout<<"Face "<<cnt++<<" ID "<<getMdsIndex(mesh, e)<<", local ID "<<get_ent_localid(mesh, e)<<"\n";
  }
  mesh->end(ent_it); 
*/
#endif

  if (!PCU_Comm_Self())
    cout<<"#global_ent: V "<<m3dc1_mesh::instance()->num_global_ent[0]
        <<", E "<<m3dc1_mesh::instance()->num_global_ent[1]
        <<", F "<<m3dc1_mesh::instance()->num_global_ent[2]
        <<", R "<<m3dc1_mesh::instance()->num_global_ent[3]<<"\n";
}

// NEW implementation of build3d for adapted mesh
// *********************************************************
void m3dc1_receiveVertices(Mesh2* mesh, std::vector<MeshEntity*>* vertices,
MeshTag* partbdry_id_tag, 
                           std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  int proc_grp_rank = PCU_Comm_Self()/m3dc1_model::instance()->group_size;
  int proc_grp_size = m3dc1_model::instance()->group_size;
  int myrank = PCU_Comm_Self();
  int num_ent, num_remote, own_partid;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      double* v_coords = new double[num_ent*3];
      double* v_params = new double[num_ent*3];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];

      PCU_Comm_Unpack(&(v_coords[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(v_params[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid;
      if (num_remote) 
      {
        e_rmt_pid = new int[num_remote];
        PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));
      }
      // create vertex
      int rinfo_pos=0;
      apf::Vector3 coord;
      apf::Vector3 param;     
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
        for (int k=0; k<3; ++k)
        {
          coord[k] = v_coords[index*3+k];
          param[k] = v_params[index*3+k];
        }
        new_ent = mesh->createVertex((ModelEntity*)geom_ent, coord, param);
        vertices->push_back(new_ent);

        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        own_partid=proc_grp_rank*proc_grp_size+e_own_partid[index];
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &own_partid);

        if (own_partid==myrank)
          ++(m3dc1_mesh::instance()->num_own_ent[0]);

        if (e_global_id[index]!=-1)
        {
          partbdry_entities[0][e_global_id[index]] = new_ent;
          mesh->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }

        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          mesh->addRemote(new_ent, proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index
      delete [] v_coords;
      delete [] v_params;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      if (num_remote)
        delete [] e_rmt_pid;
      m3dc1_mesh::instance()->num_local_ent[0] = mesh->count(0);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void m3dc1_receiveEdges(Mesh2* mesh, std::vector<MeshEntity*>* vertices,
        std::vector<MeshEntity*>* edges,
        MeshTag* partbdry_id_tag, std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  int myrank=PCU_Comm_Self();
  int proc_grp_rank = myrank/m3dc1_model::instance()->group_size;
  int proc_grp_size = m3dc1_model::instance()->group_size;

  int num_ent, num_remote, own_partid;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  Downward down_ent; 
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      int* e_down_lid = new int[num_ent*2];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];
      PCU_Comm_Unpack(&(e_down_lid[0]), num_ent*2*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid = new int[num_remote];
      if (num_remote) PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));

      // create edge
      int rinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
       //std::cout<<"#vertices "<<mesh->count(0)<<" vertices.size "<<vertices.size()
       //          <<", e_down_lid["<<index*2<<"] "<< e_down_lid[index*2]<<", e_down_lid["<<index*2+1<<"] "<<e_down_lid[index*2+1]<<"\n";
        down_ent[0] =  vertices->at(e_down_lid[index*2]);
        down_ent[1] =  vertices->at(e_down_lid[index*2+1]);
        new_ent = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)geom_ent, down_ent);
        edges->push_back(new_ent);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        own_partid=proc_grp_rank*proc_grp_size+e_own_partid[index];
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &own_partid);

        if (own_partid==myrank)
          ++(m3dc1_mesh::instance()->num_own_ent[1]);

        if (e_global_id[index]!=-1)
        {
          partbdry_entities[1][e_global_id[index]] = new_ent;
          mesh->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }
        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          //ent_bps[new_ent].insert(proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos]);
          mesh->addRemote(new_ent, proc_grp_rank*proc_grp_size+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index      
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      delete [] e_rmt_pid;
      m3dc1_mesh::instance()->num_local_ent[1] = mesh->count(1);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())

}

// *********************************************************
void m3dc1_receiveFaces(Mesh2* mesh, std::vector<MeshEntity*>* edges, std::vector<MeshEntity*>* faces)
// *********************************************************
{
  int num_ent, num_down;
  MeshEntity* new_ent;
  gmi_ent* geom_ent;
  Downward down_ent;
  int myrank = PCU_Comm_Self();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      PCU_Comm_Unpack(&num_down, sizeof(int));

      int* f_down_num = new int[num_ent];
      int* e_down_lid = new int[num_down];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      PCU_Comm_Unpack(&(f_down_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_down_lid[0]), num_down*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      // create face
      int dinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent =gmi_find(m3dc1_model::instance()->model, e_geom_type[index], e_geom_tag[index]);
        num_down = f_down_num[index];
        for (int i=0; i<num_down; ++i)
        {
          down_ent[i] = edges->at(e_down_lid[dinfo_pos]);
          ++dinfo_pos;
        }
        if (num_down==3)
          new_ent= mesh->createEntity(apf::Mesh::TRIANGLE, (ModelEntity*)geom_ent, down_ent);
        else
          new_ent= mesh->createEntity(apf::Mesh::QUAD, (ModelEntity*)geom_ent, down_ent);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->local_entid_tag, &index);
        mesh->setIntTag(new_ent, m3dc1_mesh::instance()->own_partid_tag, &myrank);
        faces->push_back(new_ent);
      } // for index
      delete [] f_down_num;
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      m3dc1_mesh::instance()->num_local_ent[2] = mesh->count(2);
      m3dc1_mesh::instance()->num_own_ent[2] = mesh->count(2);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

/*
void set_remote(Mesh2* m, MeshEntity* e, int p, MeshEntity* r)
{
  apf::Copies remotes;
  m->getRemotes(e,remotes);
  bool found=false;
  APF_ITERATE(apf::Copies, remotes, it)
    if (it->first==p) 
    {
      found=true;
      break;
    }
  if (found)
  {
    remotes[p] = r;
    m->setRemotes(e,remotes);
  }
  else
    m->addRemote(e, p, r);
}
// *********************************************************
void m3dc1_stitchLink(Mesh2* mesh, MeshTag* partbdry_id_tag, 
                      std::map<int, MeshEntity*>* partbdry_entities)
// *********************************************************
{
  PCU_Comm_Begin();
  MeshEntity* e;
  
  int ent_info[3];
  int global_id; 
  for (int dim=0; dim<2; ++dim)
  {
    for (std::map<int, MeshEntity*>::iterator ent_it = partbdry_entities[dim].begin();
         ent_it!= partbdry_entities[dim].end(); ++ent_it)
    {    
      e = ent_it->second;
      apf::Copies remotes;
      mesh->getRemotes(e, remotes);
      mesh->getIntTag(e, partbdry_id_tag, &global_id);
      ent_info[0] = dim;
      ent_info[1] = global_id; 
      ent_info[2] = PCU_Comm_Self(); 
      APF_ITERATE(apf::Copies, remotes, it)
      {  
        PCU_Comm_Pack(it->first, &(ent_info[0]),3*sizeof(int));
        PCU_COMM_PACK(it->first, e);
      }
    }
  } //  for (int dim=0; dim<2; ++dim)

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      int sender_ent_info[3];
      MeshEntity* sender;
      PCU_Comm_Unpack(&(sender_ent_info[0]), 3*sizeof(int));
      PCU_COMM_UNPACK(sender);
      e = partbdry_entities[sender_ent_info[0]][sender_ent_info[1]];
      set_remote(mesh, e, sender_ent_info[2], sender);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}
*/

// *********************************************************
void sendEntities(Mesh2* mesh, int dim, 
std::vector<MeshEntity*>* vertices,
std::vector<MeshEntity*>* edges,
std::vector<MeshEntity*>* faces,
MeshTag* partbdry_id_tag)
// *********************************************************
{
  int own_partid, num_ent = m3dc1_mesh::instance()->num_local_ent[dim];
  if (dim>2 || !num_ent) return;

  double* v_coords;
  double* v_params;
  int* e_geom_type = new int[num_ent];
  int* e_geom_tag = new int[num_ent];
  int* f_down_num; // for faces
  int* e_own_partid; // for vertices and edges
  int* e_global_id; // for vertices and edges
  int* e_rmt_num;  // for vertices and edges
  
  MeshEntity* e;
  gmi_ent* geom_ent;

  if (dim==0)
  {
    v_coords = new double[num_ent*3];  // for vertices
    v_params = new double[num_ent*3];
  }
  if (dim==2)
    f_down_num = new int[num_ent]; // for faces
  else // for vertices and edges
  {
    e_own_partid = new int[num_ent]; 
    e_global_id = new int[num_ent]; 
    e_rmt_num = new int[num_ent];  
  }

  int num_down=0, num_remote=0;
  MeshIterator* it = mesh->begin(dim);
  while ((e = mesh->iterate(it)))
  {
    switch (mesh->getType(e))
    {
      case apf::Mesh::TRIANGLE: num_down+=3; break;
      case apf::Mesh::QUAD: num_down+=4; break;
      case apf::Mesh::EDGE: num_down+=2; break;
      default: break;
    }

    if (mesh->isShared(e))
    {
      apf::Copies remotes;
      mesh->getRemotes(e, remotes);
      num_remote+=remotes.size();
    }
  }
  mesh->end(it);

  std::vector<int> e_down_lid;
  int* e_rmt_pid = new int[num_remote];

  int global_id, rinfo_pos=0, index=0, num_down_ent;

  it = mesh->begin(dim);
  apf::Vector3 coord;
  apf::Vector3 param;
  Downward down_ent;
  while ((e = mesh->iterate(it)))
  {
    switch (getDimension(mesh, e))
    {
      case 0: mesh->getPoint(e, 0, coord);
              mesh->getParam(e, param);
              for (int i=0; i<3; ++i)
              {
                v_coords[index*3+i] = coord[i];
                v_params[index*3+i] = param[i];
              }
              vertices->push_back(e);
              break;
       default: { 
                  num_down_ent =  mesh->getDownward(e, dim-1, down_ent); 
                  if (dim==2) 
                  {
                    faces->push_back(e);
                    f_down_num[index] = num_down_ent;
                  }
                  else
                    edges->push_back(e);
                  for (int i=0; i<num_down_ent; ++i)
                    e_down_lid.push_back(get_ent_localid(mesh, down_ent[i]));
                }
    } // switch

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    e_geom_type[index] = gmi_dim(m3dc1_model::instance()->model, geom_ent);
    e_geom_tag[index] = gmi_tag(m3dc1_model::instance()->model, geom_ent);

    if (dim!=2) // for vertices and edges
    {
      own_partid=get_ent_ownpartid(mesh,e);
      e_own_partid[index] = own_partid;
      apf::Copies remotes;
      mesh->getRemotes(e, remotes);
      if (mesh->hasTag(e,partbdry_id_tag)) // part-bdry entity
      {
        mesh->getIntTag(e, partbdry_id_tag, &global_id);
        e_global_id[index] = global_id; 
        e_rmt_num[index] = remotes.size();
        APF_ITERATE(apf::Copies, remotes, it)
        {
          e_rmt_pid[rinfo_pos]=it->first;
          ++rinfo_pos;
        }
      }
      else
      {
        e_global_id[index] = -1; 
        e_rmt_num[index] = 0;
      }
    }
    ++index;
  }
  mesh->end(it);

  int proc=PCU_Comm_Self()+m3dc1_model::instance()->group_size;
  while (proc<PCU_Comm_Peers())
  {  
    PCU_Comm_Pack(proc, &num_ent, sizeof(int));
    switch (dim)
    { 
      case 0: PCU_Comm_Pack(proc, &(v_coords[0]),num_ent*3*sizeof(double));
              PCU_Comm_Pack(proc, &(v_params[0]),num_ent*3*sizeof(double));
              break;
      default: if (dim==2)
               {
                 PCU_Comm_Pack(proc, &num_down, sizeof(int));
                 PCU_Comm_Pack(proc, &(f_down_num[0]),num_ent*sizeof(int));
               }
               PCU_Comm_Pack(proc, &e_down_lid.at(0), num_down*sizeof(int));
    }
    PCU_Comm_Pack(proc, &(e_geom_type[0]),num_ent*sizeof(int));
    PCU_Comm_Pack(proc, &(e_geom_tag[0]),num_ent*sizeof(int));
    if (dim<2)  
    {      
      PCU_Comm_Pack(proc, &(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_global_id[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_rmt_num[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &num_remote, sizeof(int));
      if (num_remote)
        PCU_Comm_Pack(proc, &(e_rmt_pid[0]), num_remote*sizeof(int));
    }
    proc+=m3dc1_model::instance()->group_size;
  }

  if (dim==0)
  {
    delete [] v_coords;
    delete [] v_params;
  }
  delete [] e_geom_type;
  delete [] e_geom_tag;

  if (dim==2)
    delete [] f_down_num;
  else 
  {    
    delete [] e_own_partid;
    delete [] e_global_id;
    delete [] e_rmt_num;
    if (num_remote) delete [] e_rmt_pid;
  }
}

// create wedges between planes
void create_wedges(apf::Mesh2* mesh, MeshTag* local_entid_tag, 
     std::vector<MeshEntity*>* vertices, std::vector<MeshEntity*>* edges,
     std::vector<MeshEntity*>* faces,
     MeshEntity** remote_vertices, MeshEntity** remote_edges, MeshEntity** remote_faces, 
     std::vector<MeshEntity*>& btw_plane_edges, std::vector<MeshEntity*>& btw_plane_faces, 
     std::vector<MeshEntity*>& btw_plane_regions)
{
  int local_partid=PCU_Comm_Self();
  MeshEntity* e;
  MeshEntity* new_ent;
  gmi_ent* geom_ent; 

  int new_id, local_id, local_id_0, local_id_1, local_id_2;
  double local_plane_phi = m3dc1_model::instance()->get_phi(local_partid);
  double next_plane_phi = m3dc1_model::instance()->get_phi(m3dc1_model::instance()->next_plane_partid);
  bool flip_wedge=false;
  if (local_plane_phi>next_plane_phi) // last plane
    flip_wedge=true;

  // create remote copy of vertices on next plane
  gmi_ent* new_geom_ent = NULL;
  int num_local_vtx=mesh->count(0);
  apf::Vector3 cur_coord;
  apf::Vector3 new_coord;

  for (int index=0; index<vertices->size(); ++index)
  {
    e = vertices->at(index);
    mesh->getPoint(e, 0, cur_coord);
    new_coord[0]=cur_coord[0];
    new_coord[1]=cur_coord[1];
    new_coord[2]=next_plane_phi;
    apf::Vector3 param(0,0,0);     
    mesh->getParam(e,param);
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    assert(geom_ent);
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);
    if (!new_geom_ent)
      std::cout<<"["<<PCU_Comm_Self()<<"] new_geom_ent for dim "<<gmi_dim(m3dc1_model::instance()->model, geom_ent)
               <<", tag "<<gmi_tag(m3dc1_model::instance()->model, geom_ent)<<" not found\n";
    assert(new_geom_ent);

    new_ent = mesh->createVertex((ModelEntity*)new_geom_ent, new_coord, param);
    new_id = num_local_vtx+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_vertices[index]=new_ent;
    // update the z-coord of the current coordinate based on the prev plane partid
    cur_coord[2] = local_plane_phi;
    mesh->setPoint(e, 0, cur_coord);
    if (index==num_local_vtx) break;
  }

  // create remote copy of edges on next plane
  int num_local_edge=mesh->count(1);
  Downward down_vtx;
  Downward new_down_vtx;
 
  for (int index=0; index<edges->size(); ++index)
  {  
    e = edges->at(index);
    mesh->getDownward(e, 0, down_vtx);
    local_id_0 = get_ent_localid(mesh, down_vtx[0]);
    local_id_1 = get_ent_localid(mesh, down_vtx[1]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);
    if (!new_geom_ent)
      std::cout<<"["<<PCU_Comm_Self()<<"] new_geom_ent for dim "<<gmi_dim(m3dc1_model::instance()->model, geom_ent)
               <<", tag "<<gmi_tag(m3dc1_model::instance()->model, geom_ent)<<" not found\n";
    assert(new_geom_ent);

    new_down_vtx[0] = remote_vertices[local_id_0];
    new_down_vtx[1] = remote_vertices[local_id_1];

    new_ent = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);
    new_id = num_local_edge+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_edges[index]=new_ent;
  }

  // create remote copy of faces on next plane
  int num_local_face=mesh->count(2);
  Downward down_edge;
  Downward new_down_edge;

  for (int index=0; index<faces->size(); ++index)
  {  
    e = faces->at(index);
    mesh->getDownward(e, 1, down_edge);
    local_id_0 = get_ent_localid(mesh, down_edge[0]);
    local_id_1 = get_ent_localid(mesh, down_edge[1]);
    local_id_2 = get_ent_localid(mesh, down_edge[2]);

    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntNextPlane(geom_ent);

    new_down_edge[0] = remote_edges[local_id_0];
    new_down_edge[1] = remote_edges[local_id_1];
    new_down_edge[2] = remote_edges[local_id_2];

    new_ent = mesh->createEntity(apf::Mesh::TRIANGLE, (ModelEntity*)new_geom_ent, new_down_edge);
    new_id = num_local_face+index;
    mesh->setIntTag(new_ent, local_entid_tag, &new_id);
    remote_faces[index]=new_ent;
  }

  // create edges, faces and regions between planes. 
  // if edges are classified on geom edge, a new face is classified on geom_surface_face
  // otherwise, a new face is classified on geom_region

  Downward quad_faces;
  Downward wedge_faces;

  int edge_counter=0, face_counter=0, rgn_counter=0;
  std::vector<MeshEntity*> vertex_edges;
  std::vector<MeshEntity*> edge_faces;

  Downward edgesNextPlane;
  Downward edgesBtwPlane;
  int num_upward;

  for (int index=0; index<faces->size(); ++index)
  {
    e = faces->at(index);
    mesh->getDownward(e, 0, down_vtx);
    mesh->getDownward(e, 1, down_edge);

    for (int pos=0; pos<3; ++pos)
    {
      local_id = get_ent_localid(mesh, down_edge[pos]);
      edgesNextPlane[pos] = remote_edges[local_id];
    }

    /**create edges between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the new created face btw planes*/
      edgesBtwPlane[pos]=NULL;
      vertex_edges.clear();
      num_upward = mesh->countUpward(down_vtx[pos]);
      for (int i=0; i<num_upward; ++i)
        vertex_edges.push_back(mesh->getUpward(down_vtx[pos], i));

      for (unsigned int i=0; i<vertex_edges.size();++i)
      {
        // get the local id
        local_id = get_ent_localid(mesh, vertex_edges[i]);
        if (local_id>=num_local_edge) // edge is between planes
        {
          edgesBtwPlane[pos]=vertex_edges[i];
          break; // get out of for loop
        }
      }

      if (edgesBtwPlane[pos]!=NULL) continue;

      // create new edges between vertices[pos] and its remote vertex if not found
      local_id = get_ent_localid(mesh, down_vtx[pos]);

      geom_ent = (gmi_ent*)(mesh->toModel(down_vtx[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
      new_down_vtx[0] = down_vtx[pos];
      new_down_vtx[1] = remote_vertices[local_id],

      edgesBtwPlane[pos] = mesh->createEntity(apf::Mesh::EDGE, (ModelEntity*)new_geom_ent, new_down_vtx);

      new_id =num_local_edge*2+edge_counter;
      mesh->setIntTag(edgesBtwPlane[pos], local_entid_tag, &new_id);
      btw_plane_edges.push_back(edgesBtwPlane[pos]);
      ++edge_counter;
    }// for (int pos=0; pos<3; ++pos)

    /**create quads between planes*/
    for (int pos=0; pos<3; ++pos)
    {
      /** seek all the faces of the edge, find the newly created face btw planes*/
      quad_faces[pos]=NULL;
      edge_faces.clear();

      num_upward = mesh->countUpward(down_edge[pos]);

      for (int i=0; i<num_upward; ++i)
      {
        edge_faces.push_back(mesh->getUpward(down_edge[pos],i));
      }
      for (unsigned int i=0; i<edge_faces.size(); ++i)
      {
        local_id = get_ent_localid(mesh, edge_faces[i]);
        if (local_id>=num_local_face) // face is between planes
        {
          quad_faces[pos]=edge_faces[i];
          break; // get out of for loop
        }
      }
      if (quad_faces[pos]!=NULL) continue;

      /**create new quad between two planes if found*/
      Downward quad_edges;
      quad_edges[0]=down_edge[pos];
      Downward down_edge_vtx;
      mesh->getDownward(down_edge[pos], 0, down_edge_vtx);
      MeshEntity* vtx_1=down_edge_vtx[1];
      for (int i=0; i<3; ++i)
      { 
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_1 || edgeBtw_down[1]==vtx_1)
        {
          quad_edges[1]=edgesBtwPlane[i];
          break; // get out of for loop
        }
      }

      quad_edges[2]=edgesNextPlane[pos];
      MeshEntity* vtx_2 = down_edge_vtx[0];

      for (int i=0; i<3; ++i)
      {
        Downward edgeBtw_down;
        mesh->getDownward(edgesBtwPlane[i], 0, edgeBtw_down);
        if (edgeBtw_down[0]==vtx_2 || edgeBtw_down[1]==vtx_2)
        {
          quad_edges[3]=edgesBtwPlane[i];
          break;// get out of for loop
        }
      }
      geom_ent = (gmi_ent*)(mesh->toModel(down_edge[pos]));
      new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);

      quad_faces[pos]= mesh->createEntity(apf::Mesh::QUAD, (ModelEntity*)new_geom_ent, quad_edges);
      btw_plane_faces.push_back(quad_faces[pos]);
      new_id = num_local_face*2+face_counter;
      mesh->setIntTag(quad_faces[pos], local_entid_tag, &new_id);
      ++face_counter;
    } //     for (int pos=0;pos<3; ++pos)

    // create regions per face on local plane
    wedge_faces[0]=e;
    wedge_faces[1]=quad_faces[0];
    wedge_faces[2]=quad_faces[1];
    wedge_faces[3]=quad_faces[2];
    local_id = get_ent_localid(mesh, e);
    wedge_faces[4]=remote_faces[local_id];

    if (flip_wedge) // flip top & bottom of wedge to avoid negative volume in the last plane
    {
      wedge_faces[4]=e;
      wedge_faces[0]=remote_faces[local_id];
    }

    /**create new region between two planes*/
    geom_ent = (gmi_ent*)(mesh->toModel(e));
    new_geom_ent = m3dc1_model::instance()->geomEntBtwPlane(geom_ent);
//    std::cout<<"[M3D-C1 INFO] (p"<<PCU_Comm_Self()<<") create PRISM with face "<<get_ent_localid(mesh, wedge_faces[0])<<"(t="<< mesh->getType(wedge_faces[0])<<"), "<<get_ent_localid(mesh, wedge_faces[1])<<"(t="<< mesh->getType(wedge_faces[1])<<"), "<<get_ent_localid(mesh, wedge_faces[2])<<"(t="<< mesh->getType(wedge_faces[2])<<"), "<<get_ent_localid(mesh, wedge_faces[3])<<"(t="<< mesh->getType(wedge_faces[3])<<"), "<<get_ent_localid(mesh, wedge_faces[4])<<"(t="<< mesh->getType(wedge_faces[4])<<") (#face="<<mesh->count(2)<<")"<<std::endl;
    new_ent = mesh->createEntity(apf::Mesh::PRISM, (ModelEntity*)new_geom_ent, wedge_faces);
    btw_plane_regions.push_back(new_ent);
    mesh->setIntTag(new_ent, local_entid_tag, &rgn_counter);
    ++rgn_counter;
  }
}

void create_localid(apf::Mesh2* mesh, int dim)
{
  MeshIterator* ent_it = mesh->begin(dim);
  apf::MeshEntity* e;

  // assign sequential local ID starting with 0
  int index=0;
  while ((e = mesh->iterate(ent_it)))
  {
    mesh->setIntTag(e, m3dc1_mesh::instance()->local_entid_tag, &index);
    ++index;
  }
  mesh->end(ent_it);
}

// *********************************************************
void m3dc1_mesh::restore3D()
// *********************************************************
{
    /* // compute the fields to copy from master plane */
    int num_field = 0;
    int* field_id = NULL;
    int* num_dofs_per_value = NULL;

    /* if (!m3dc1_model::instance()->local_planeid && m3dc1_mesh::instance()->field_container) */
    /* { */
    /*   num_field = m3dc1_mesh::instance()->field_container->size(); */
    /*   cout<<__func__<<": #fields to copy from master to non-master plane "<<num_field<<"\n"; */
    /*   field_id = new int [num_field]; */
    /*   num_dofs_per_value = new int [num_field]; */
    /*   for (int i=0; i<num_field; ++i) */
    /*   { */
    /*     m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id]; */
    /*     field_id[i] = i; */
    /*     num_dofs_per_value[i] = mf->get_dof_per_value(); */
    /*   } */
    /* } */

  int local_partid=PCU_Comm_Self();

  MeshIterator* ent_it;
  apf::MeshEntity* e;
  int index; 

  set_mcount();  // reset mesh entity count

  // assign uniq id to part bdry entities
  MeshTag* partbdry_id_tag = mesh->findTag("m3dc1_pbdry_globid");
  if (!partbdry_id_tag)
   partbdry_id_tag = mesh->createIntTag("m3dc1_pbdry_globid", 1);

  for (int dim=0; dim<2; ++dim)
    assign_uniq_partbdry_id(mesh, dim, partbdry_id_tag);

  // copy 2D mesh in process group 0 to other process groups
  std::map<int, MeshEntity*> partbdry_entities[2];

  std::vector<apf::MeshEntity*>* vertices = new std::vector<apf::MeshEntity*>;
  std::vector<apf::MeshEntity*>* edges = new std::vector<apf::MeshEntity*>;
  std::vector<apf::MeshEntity*>* faces = new std::vector<apf::MeshEntity*>;

  for (int dim=0; dim<3; ++dim)
  {
    if (mesh->count(dim))  create_localid (mesh, dim);
    PCU_Comm_Begin();
    sendEntities(mesh, dim, vertices, edges, faces, partbdry_id_tag);
    PCU_Comm_Send();
    switch (dim)
    {
      case 0: m3dc1_receiveVertices(mesh, vertices, partbdry_id_tag, partbdry_entities); break;
      case 1: m3dc1_receiveEdges(mesh, vertices, edges, partbdry_id_tag, partbdry_entities); break;
      case 2: m3dc1_receiveFaces(mesh, edges, faces); break;
      default: break;
    }   
  }

  m3dc1_stitchLink(mesh, partbdry_id_tag, partbdry_entities);

  for (int dim=0; dim<2; ++dim)
  {
    ent_it = mesh->begin(dim);
    while ((e = mesh->iterate(ent_it)))
    {
      if (!mesh->isShared(e)) continue;
      apf::Copies remotes;
      apf::Parts parts;
      mesh->getRemotes(e, remotes);
      APF_ITERATE(apf::Copies, remotes, it)
        parts.insert(it->first);
      parts.insert(local_partid);
      mesh->setResidence (e, parts); // set pclassification
      mesh->removeTag(e, partbdry_id_tag);
    }
    mesh->end(ent_it);
  }
  mesh->destroyTag(partbdry_id_tag);

  // update global ent counter
  MPI_Allreduce(num_own_ent, num_global_ent, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // construct 3D model
  changeMdsDimension(mesh, 3);
  

  m3dc1_model::instance()->restore_gtag();

  int num_local_vtx = num_local_ent[0];
  int num_local_edge = num_local_ent[1];
  int num_local_face = num_local_ent[2]; 

  MeshEntity** remote_vertices=new MeshEntity*[mesh->count(0)];
  MeshEntity** remote_edges=new MeshEntity*[mesh->count(1)];
  MeshEntity** remote_faces=new MeshEntity*[mesh->count(2)];
  std::vector<MeshEntity*> btw_plane_edges;
  std::vector<MeshEntity*> btw_plane_faces;
  std::vector<MeshEntity*> btw_plane_regions;

  create_wedges(mesh, local_entid_tag, vertices, edges, faces, 
      remote_vertices, remote_edges, remote_faces, 
      btw_plane_edges, btw_plane_faces, btw_plane_regions);


  // exchange remote copies to set remote copy links
  PCU_Comm_Begin();
  int num_entities = mesh->count(0)+mesh->count(1)+mesh->count(2);
  MeshEntity** entities = new MeshEntity*[num_entities];
  int pos=0;
  for (int index=0; index<num_local_vtx; ++index)
    entities[pos++]=remote_vertices[index];
  for (int index=0; index<num_local_edge; ++index)
    entities[pos++]=remote_edges[index];

  for (int index=0; index<num_local_face; ++index)
    entities[pos++]=remote_faces[index];

  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &num_entities, sizeof(int));
  PCU_Comm_Pack(m3dc1_model::instance()->next_plane_partid, &(entities[0]),num_entities*sizeof(MeshEntity*));
  PCU_Comm_Send();
  delete [] entities;

  // receive phase begins
  std::vector<MeshEntity*> ent_vec;
  std::map<MeshEntity*, MeshEntity*> partbdry_ent_map;
  while (PCU_Comm_Listen())
  {
    int index=0;
    while (!PCU_Comm_Unpacked())
    {
      int num_ent;
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      MeshEntity** s_ent = new MeshEntity*[num_ent];
      PCU_Comm_Unpack(&(s_ent[0]), num_ent*sizeof(MeshEntity*));
      pos=0;
      //index=0;
      //ent_it = mesh->begin(0);
      //while ((e = mesh->iterate(ent_it)))
      for (int index=0; index<vertices->size(); ++index)
      {
        e = vertices->at(index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }
      
      
      for (int index=0; index<edges->size(); ++index)
      {
        e = edges->at(index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        if (mesh->isShared(e))
          partbdry_ent_map[e]=s_ent[pos];
        ent_vec.push_back(e);
        ++pos;
      }

      for (int index=0; index<faces->size(); ++index)
      {
        e = faces->at(index);
        mesh->addRemote(e, m3dc1_model::instance()->prev_plane_partid, s_ent[pos]);
        ent_vec.push_back(e);
        ++pos;
      }
      delete [] s_ent;
    }
  }

  push_new_entities(mesh, partbdry_ent_map);

  bounce_orig_entities(mesh, ent_vec, m3dc1_model::instance()->prev_plane_partid,
                       remote_vertices,remote_edges,remote_faces);

  // update partition classification
  for (int index=0; index<vertices->size(); ++index)
  {
    e = vertices->at(index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_vtx; ++index)
  {
    e = remote_vertices[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<edges->size(); ++index)
  {
    e = edges->at(index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first); 
    parts.insert(local_partid);  
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_edge; ++index)
  {
    e = remote_edges[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<faces->size(); ++index)
  {
    e = faces->at(index);
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  for (int index=0; index<num_local_face; ++index)
  {
    e = remote_faces[index];
    apf::Copies remotes;
    apf::Parts parts;
    mesh->getRemotes(e, remotes);
    APF_ITERATE(apf::Copies, remotes, it)
      parts.insert(it->first);   
    parts.insert(local_partid);
    mesh->setResidence(e, parts);
  }

  // clear entity containers with creation order preserved
  vertices->clear();
  edges->clear();
  faces->clear();
  delete vertices;
  delete edges;
  delete faces;

  // update partition model
  mesh->acceptChanges();

  // update m3dc1_mesh internal data (# local, owned and global entities)
  update_partbdry(remote_vertices, remote_edges, remote_faces, btw_plane_edges, btw_plane_faces, btw_plane_regions);

  // delete existing local numbering
  apf::Numbering* local_n = mesh->findNumbering(mesh->getShape()->getName());
  if (local_n) destroyNumbering(local_n);
  
  /* // FIXME: re-create the field and copy field data on master process group to non-master */
  /* for (int i=0; i<num_field; ++i) */
  /*   update_field(field_id[i], num_dofs_per_value[i], num_local_vtx, remote_vertices); */

  // clear temp memory
  delete [] remote_vertices;
  delete [] remote_edges;
  delete [] remote_faces;

  if (num_field>0)
  {
    delete [] field_id;
    delete [] num_dofs_per_value;
  }

  set_node_adj_tag();

  compute_globalid(mesh, 0);
  compute_globalid(mesh, 3);

  if (!PCU_Comm_Self())
    std::cout<<"\n*** MESH SWITCHED TO 3D ***\n\n";
}

