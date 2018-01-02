/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <apf.h>
#include <apfPUMI.h>
#include <pumi.h>
#include <pumi_mesh.h>
#include <apfNumberingAlgorithms.h>
#include <iostream>
using namespace std;
pMeshMdl globalMesh;
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (!PCU_Comm_Self())
  {
    cout<<"Usage: ./main  inputMesh partioneded(0 for no, 1 for yes) outputMesh(.smb)"<<endl;
  }  
  assert(argc>=4);
  const char* meshFile = argv[1];
  PUMI_Init(MPI_COMM_WORLD);
  pMeshMdl mesh;
  PUMI_Mesh_Create(NULL,mesh);
  PUMI_Mesh_LoadFromFile(mesh,meshFile,atoi(argv[2]));
  FMDB_Mesh_WriteToFile (mesh, "part.sms", 1);
  apf::Mesh2* m = apf::createMesh(mesh);
  m->writeNative(argv[3]);
  if(argc==5) apf::writeVtkFiles("mesh",m);
  apf::destroyMesh(m);
  PUMI_Mesh_Del(mesh);
  PUMI_Finalize();
}
