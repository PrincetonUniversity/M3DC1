/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "ma.h"
#include "PCU.h"
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  if (argc!=3 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh "<<endl;
    return M3DC1_FAILURE;
  } 
  const char* meshFile = argv[1];
  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]); 
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
  double center[3]={1.75, 0, 0.};
  //double param[13]={0.98, 2, 1, .05, .5, .05, .5, .1, .01, 100., 100., .07, .16};
  double param[13]={0.98, 2, 1, .05, .5, .05, .5, .05, .01, 100., 100., .07, .16};

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
    ReducedQuinticTransfer slnTrans(mesh, fields, &shape);
    double psi0 = 0.48095979306833486;
    double psil = 0.20875939867733129;
    //double param[13]={0.98, 2, 1, .05, .5, .05, .5, .1, .01, 100., 100., .07, .16};
    double param[13]={0.98, 2, 1, .05, .5, .05, .5, .05, .01, 100., 100., .07, .16};
    SizeFieldPsi sf (psiField, psi0, psil, param, center);
    ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
  in->shouldRunPreZoltan = true;
    in->maximumIterations = 9;
    ma::adapt(in);
  }
  else
  {
    double mmax[2], mmin[2];
    m3dc1_model_getmaxcoord(mmax,mmax+1);
    m3dc1_model_getmincoord(mmin,mmin+1);
    Vortex sfv(mesh, center, mmax[0]-mmin[0]);
    ma::Input* in = ma::configure(mesh,&sfv,0);
    in->maximumIterations = 9;
  in->shouldRunPreZoltan = true;
    ma::adapt(in);
  }
//  in->shouldFixShape=0;
  //in->goodQuality = 0.5;
//check(m);
  //mesh->verify();
  apf::writeVtkFiles("after",mesh);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  MPI_Finalize();
  return 0;
}
