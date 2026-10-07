/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <ma.h>
#include <apf.h>
#include <gmi_analytic.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>

double  a_param = 0;
double  b_param = 0;
double  c_param = 0;
double  d_param = 0;
double  e_param = 0;
double meshSize=0.1;
void getCenter(double x[3])
{
   x[0]=a_param;
   x[1]=d_param;
   x[2]=0.;
}
void edgeFunction(double const p[2], double x[3], void*)
{
  double phi = p[0];
  x[0] = a_param + b_param*(cos(phi + c_param*sin(phi)));
  x[1] = d_param + e_param*sin(phi);
  x[2] = 0;
}

void snap2Model(apf::Mesh2* mesh, gmi_model* model)
{
  for(int i=0; i<5; i++)
  {
    apf::MeshEntity* vt = getMdsEntity(mesh, 0, i);
    assert(vt);
    gmi_ent* gent= (gmi_ent*)(mesh->toModel(vt));
    apf::Vector3 param(0,0,0);
    int gType = gmi_dim(model,gent);
    double xyz[3];
    if(gType==1)
    {
      mesh->getParam(vt,param);
      double pa[]={param[0],0,0};
      edgeFunction(pa, xyz, 0);
    }
    else
    {
      assert(gType==2);
      getCenter(xyz);
    }
    mesh->setPoint(vt, 0, xyz);
  }
}
void faceFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}

gmi_model* makeModel(char* modelfile)
{
  FILE* fp=fopen(modelfile,"r");
  if(!fp)
  {
    std::cout<<"SCOREC ERROR: fail to open the file \""<<modelfile<<"\"\n";
    return NULL;
  }
  fscanf(fp, "%lf %lf %lf %lf %lf\n", &a_param, &b_param, &c_param, &d_param, &e_param);
  fclose(fp);
  gmi_model* model = gmi_make_analytic();
  int edgePeriodic = 1;
  double edgeRange[2] = {0, 2 * apf::pi};
  gmi_add_analytic(model, 1, 1, edgeFunction, &edgePeriodic, &edgeRange, 0);
  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,0},{0,0}};
  gmi_add_analytic(model, 2, 1, faceFunction, facePeriodic, faceRanges, 0);
  return model;
}

class Size : public ma:: IsotropicFunction
{
  public:
    Size(ma::Mesh* m){}
    virtual double getValue(ma::Entity* vert) {return meshSize;};
};

static void testIndexing(apf::Mesh2* m)
{
  for (int d = 0; d <= m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    int i = 0;
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      assert( apf::getMdsIndex(m, e) == i );
      assert( apf::getMdsEntity(m, d, i) == e );
      ++i;
    }
    m->end(it);
  }
}

static void fusionAdapt(apf::Mesh2* m)
{
  Size sf(m);
  for(int i=0; i<10; i++)
  {
    std::cout<<" iter "<<i<<std::endl;
    #ifdef OLDMA
  	ma::Input* in = ma::configure(m, &sf);
    #else
	ma::Input* in = ma::makeAdvanced(ma::configure(m,&sf));
    #endif
    in->maximumIterations = 1;
    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    ma::adapt(in);
  }
}

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  if (argc!=3)
  {
    std::cout<<"SCOREC ERROR: wrong input arguments. \nUsage: ./create_smb model_file relative_mesh_size"
             <<"\n - model_file shall contain 5 double: a_param, b_param, c_param, d_param, e_param \n"
             <<" - given relative_mesh_size, mesh edge length is about relative_mesh_size * the longest edge of the bounding box of the model\n";
    PCU_Comm_Free();
    MPI_Finalize();    
    return 1;
  }

  char* modelfile = argv[1];
  meshSize =atof(argv[2]);

  // create model
  gmi_model* model = makeModel(modelfile);
  if (!model)
  {
    PCU_Comm_Free();
    MPI_Finalize();    
    return 1;
  }
  // generate and refine mesh
  apf::Mesh2* mesh = apf::loadMdsMesh(model, "seed.smb");
  snap2Model(mesh, model);
  mesh->verify();
  testIndexing(mesh);
  fusionAdapt(mesh);
  // export to files
  char meshfile[256];
  sprintf(meshfile,"%s.smb",modelfile); 
  mesh->writeNative(meshfile);
  apf::writeVtkFiles("Vtk-Mesh-File", mesh);
  std::cout<<"\n";
  mesh->verify();
  std::cout<<"\n";
  // clean and finalize
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
