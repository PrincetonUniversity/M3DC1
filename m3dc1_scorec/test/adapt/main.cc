/****************************************************************************** 

  (c) 2005-2021 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "m3dc1_model.h"
#include "apf.h"
#include "apfMesh.h"
#include <iostream>
#include "ma.h"
#include "pumi.h"
#include <mpi.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "petscksp.h"
#include "PCU.h"

static void print_local_and_owned_ents(int dim)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  int cnt = 0;
  int cnt_owned = 0;

  while ( (e = m->iterate(it)) )
  {
    if (m->isOwned(e)) cnt_owned++;
    cnt++;
  }
  m->end(it);

  if (dim == 0)
    printf("self %d: vert cnt  all/owned %d/%d\n", PCU_Comm_Self(), cnt, cnt_owned);
  else if (dim == 1)
    printf("self %d: edge cnt  all/owned %d/%d\n", PCU_Comm_Self(), cnt, cnt_owned);
  else if (dim == 2)
    printf("self %d: face cnt  all/owned %d/%d\n", PCU_Comm_Self(), cnt, cnt_owned);
  else
    printf("self %d: regi cnt  all/owned %d/%d\n", PCU_Comm_Self(), cnt, cnt_owned);

  if(!PCU_Comm_Self())
    printf("\n\n");

  PCU_Barrier();
}

static MPI_Comm change_comm(int numplane)
{
  MPI_Comm newComm;
  int self = PCU_Comm_Self();
  int peers = PCU_Comm_Peers();
  int groupsize = peers/numplane;
  int localpid = self/groupsize;
  int grouprank = self%groupsize;

  MPI_Comm_split(m3dc1_model::instance()->oldComm, localpid, grouprank, &newComm);
  return newComm;
}


int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init(); // this calls PCU_Comm_Init()
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  lion_set_verbosity(1);

  if (argc != 4)
    if ( !PCU_Comm_Self() )
    {
      printf("Usage: %s <model> <mesh> <#planes>\n", argv[0]);
      PetscFinalize();
      m3dc1_scorec_finalize();
      MPI_Finalize();
    }

  int zero=0, num_plane;
  num_plane = atoi(argv[3]);

  m3dc1_model_load(argv[1]);


  if (num_plane>1 && PCU_Comm_Peers()%num_plane==0)
    m3dc1_model_setnumplane (&num_plane); // this switches the Comm


  // -- load the 2D mesh
  // -- build the 3D mesh
  // -- write to vtk for debugging
  m3dc1_mesh_load(argv[2]); // load the 2d mesh
  m3dc1_mesh_build3d(&zero, &zero, &zero); // make the 3d mesh
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::writeVtkFiles("01_mesh_initial", m);

  // add a made up size field
  apf::Field* sfld = apf::createFieldOn(m, "size", apf::SCALAR);
  apf::MeshEntity* ent;
  apf::MeshIterator* it = m->begin(0);
  while ( (ent = m->iterate(it)) )
  {
    apf::Vector3 p;
    m->getPoint(ent, 0, p);
    double x = p[0];
    double s = x>2. ? 0.24 : 0.55;
    apf::setScalar(sfld, ent, 0, s);
  }
  m->end(it);

  // remove 3D and non-master plane and write to vtk for debugging
  m3dc1_mesh::instance()->remove3D();
  apf::printStats(m);

  PCU_Barrier();

  // change the communications so that there are two groups (corresponding to the original 2 parts in the 2D mesh)
  // Note this is exactly similar to what happens inside m3dc1_model_setnumplane ()
  MPI_Comm groupComm = change_comm(m3dc1_model::instance()->num_plane);
  PCU_Switch_Comm(groupComm);

  // do an adapt on the master plane
  if (m3dc1_model::instance()->local_planeid == 0)
  {
    // write the mesh before adapt
    // for this one will only have 2 .vtk files
    apf::writeVtkFiles("02_mesh_2d_before_adapt", m);
#ifdef OLDMA
    ma::Input* in = ma::configure(m, sfld);
#else
    ma::Input* in = ma::makeAdvanced(ma::configure(m, sfld));
#endif

    in->shouldSnap = false;
    in->shouldTransferParametric = false;
    ma::adapt(in);
    // write the mesh after adapt
    // for this one will only have 2 .vtk files
    apf::writeVtkFiles("03_mesh_2d_after_adapt", m);

    // clean up fields, numberings, and tags
    int numFields = m->countFields();
    int numNumberings = m->countNumberings();
    int numTags;
    apf::DynamicArray<apf::MeshTag*> tags;
    m->getTags(tags);
    numTags = tags.getSize();
    if (!PCU_Comm_Self())
      printf("there are %d/%d/%d fields/numberings/tags on the mesh\n", numFields, numNumberings, numTags);

    while (m->countFields()) {
      apf::Field* f = m->getField(0);
      if (!PCU_Comm_Self())
      	printf("removing field with name %s\n", apf::getName(f));
      m->removeField(f);
      apf::destroyField(f);
    }
    while (m->countNumberings()) {
      apf::Numbering* n = m->getNumbering(0);
      if (!PCU_Comm_Self())
      	printf("removing Numbering with name %s\n", apf::getName(n));
      m->removeNumbering(n);
      apf::destroyNumbering(n);
    }
    numFields = m->countFields();
    numNumberings = m->countNumberings();
    m->getTags(tags);
    numTags = tags.getSize();
    if (!PCU_Comm_Self())
      printf("after removing fields and numberings, there are %d/%d/%d fields/numberings/tags on the mesh\n", numFields, numNumberings, numTags);
    if (numTags > 0)
    {
      for (int i = 0; i < numTags; i++) {
	if (!PCU_Comm_Self())
	  printf("tag %d's name is %s and will be removed\n", i, m->getTagName(tags[i]));
	for (int d = 0; d < 4; d++) {
	  apf::removeTagFromDimension(m, tags[i], d);
	}
	m->destroyTag(tags[i]);
      }
      // Note this is note done in m3dc1_scorec/test/adapt_3d/main.cc but seems to be necessary
      // otherwise failure happens earlier
      m3dc1_mesh::instance()->local_entid_tag = NULL;
      m3dc1_mesh::instance()->own_partid_tag = NULL;
      m3dc1_mesh::instance()->num_global_adj_node_tag = NULL;
      m3dc1_mesh::instance()->num_own_adj_node_tag = NULL;
    }

    apf::reorderMdsMesh(m);
  }

  PCU_Barrier();

  // switch to global comm
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);

  // some stats
  print_local_and_owned_ents(0);
  print_local_and_owned_ents(1);
  print_local_and_owned_ents(2);
  print_local_and_owned_ents(3);

  m3dc1_mesh::instance()->restore3D();
  apf::printStats(m);

  apf::writeVtkFiles("04_mesh_3d_after_adapt", m);
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
