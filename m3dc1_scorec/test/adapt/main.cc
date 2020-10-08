/****************************************************************************** 

  (c) 2005-2020 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include "apfMesh.h"
#include <iostream>
#include "ma.h"
#include "pumi.h"
#include <mpi.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "petscksp.h"
#include <math.h>
#include <stdlib.h>

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

  if (argc<3 & !pumi_rank())
  {
    cout<<"Usage: ./main  model mesh #planes"<<endl;
    return M3DC1_FAILURE;
  }
 
  int num_plane=1;
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && pumi_size()%num_plane==0)  
    	m3dc1_model_setnumplane (&num_plane);

  }
  if (m3dc1_model_load(argv[1])) // model loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }

  // m3dc1_model_print();

  if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
  {
    PetscFinalize();
    m3dc1_scorec_finalize();
    MPI_Finalize();
    return 0;
  }
  
  
  // Set input for m3dc1_mesh_build3d()
  int num_field = 0;
  int field_id = 0;
  int dof_per_value = 0;
 
  if (num_plane>1)
    m3dc1_mesh_build3d(&num_field, &field_id, &dof_per_value);    // Build 3d

  pumi_mesh_print(m3dc1_mesh::instance()->mesh);
 
  double center[3]={1.75, 0, 0.};  
  apf::MeshTag* psiTag[6];
  apf::Mesh2*  mesh =m3dc1_mesh::instance()->mesh;
  int tagFound=1;

  for(int i=0; i<6; i++)
  {
  	char buff[256];
    	sprintf(buff, "psi%d", i);
    	psiTag[i]=mesh->findTag(buff);
    	if(!psiTag[i]) tagFound=0;
  }
  if(tagFound)
  {
    apf::Field* psiField=createPackedField(mesh, "psi", 6);
    apf::MeshIterator* it = mesh->begin(0);
    it = mesh->begin(0);
    while (apf::MeshEntity * e = mesh->iterate(it))
    {
      double psiValue[6];
      for(int i=0; i<6; i++)
      {
        mesh->getDoubleTag(e, psiTag[i], psiValue+i);
      }
      apf::setComponents(psiField, e,0,psiValue);
    }
    mesh->end(it);
    ReducedQuinticImplicit shape;
    vector<apf::Field*> fields;
    fields.push_back(psiField);
    //ReducedQuinticTransfer slnTrans(mesh, fields, &shape);
    double psi0 = 0.48095979306833486;
    double psil = 0.20875939867733129;
    double param[13]={0.98, 2, 1, .05, .5, .05, .5, .05, .01, 100., 100., .07, .16};
    SizeFieldPsi sf (psiField, psi0, psil, param, 2);
    std::cout << psi0 << "\n";
  }
  
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
