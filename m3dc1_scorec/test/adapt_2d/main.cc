#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include "apfMesh.h"
#include <gmi_mesh.h>
#include <iostream>
#include "ma.h"
#include "pumi.h"
#include <mpi.h>
#include <lionPrint.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "petscksp.h"
#include <math.h>
#include <stdlib.h>


// Function Prototype
int get_field (double aver,double* boundingbox, double*  pos, double &size_h1,double &size_h2, double* dir_1);

// Function to get Field (Dummy For now)
int get_field (double aver,double* boundingbox, double*  pos, double &size_h1,double &size_h2, double* dir_1)    
{
        double average = aver;
        double lower = boundingbox[0];
        double upper = boundingbox[1];
        double x = (pos[0] - lower)/(upper - lower);
        double sizeFactor = 2;
        if (x < 0.5)
                sizeFactor = 5;
	if (x >= 0.5 && x < 0.8)
		sizeFactor = 0.5;
        size_h1 = average;
        size_h2 = average/sizeFactor;
        dir_1[0] = 1.0;
        dir_1[1] = 0.0;
        dir_1[2] = 0.0;

//	std::cout << " Size " << size_h1 << " , " << size_h2 << "\n"; 
        return M3DC1_SUCCESS;
}

// Class to define the sizeField
class SetSizeField : public ma::AnisotropicFunction
{
  public:
    SetSizeField(ma::Mesh* m)
    {
        mesh = m;
        average = ma::getAverageEdgeLength(m);
        ma::getBoundingBox(m,lower,upper);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
        ma::Vector p = ma::getPosition(mesh,v);
        double pos[3] = {p[0],p[1],0.0};
        double box[2] = {lower[0],upper[0]};
        double dir_1[3];
	
	std::cout << " Lower & Upper: " << lower[0] << " , " << upper[0] << "\n";
	std::cout << "Position Coordinates : " << pos[0] << " , " << pos[1] << "\n";
	std::cout << "Average = " << average <<"\n";

        double aver = average;
        double size_h1;     
	double size_h2;     
        int field_success = get_field (aver ,box, pos, size_h1,size_h2, dir_1);

        // Calculate the second unit vector
        double a, b;
        double frac_1, frac_2;
        frac_1 = (dir_1[0])*(dir_1[0]);
        frac_2 = (dir_1[0])*(dir_1[0]) + (dir_1[1])*(dir_1[1]);


        b = sqrt (frac_1/frac_2);
        a = -(dir_1[1]*b)/dir_1[0];

        double mag = sqrt (a*a + b*b);
        double dir_2[3];
        dir_2[0] = a /mag;
        dir_2[1] = b /mag;
        dir_2[2] = 0.0;

	ma::Vector h(size_h1, size_h2,size_h2);

	R[0][0]=dir_1[0];
  	R[0][1]=dir_1[1];
  	R[0][2]=0.0;
 
  	R[1][0]= dir_2[0];
  	R[1][1]= dir_2[1];
  	R[1][2]=0.0;
	
  	R[2][0]=0;
  	R[2][1]=0;
  	R[2][2]=1.;
		
	H = h;

    }
  private:
        ma::Mesh* mesh;
        double average;
        ma::Vector lower;
        ma::Vector upper;

};


int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  m3dc1_scorec_init();
  m3dc1_model_load(argv[1]);
  m3dc1_mesh_load(argv[2]);
  lion_set_verbosity(1);
  gmi_register_mesh();
  bool logInterpolation = true;
  
  apf::Mesh2*  mesh = m3dc1_mesh::instance()->mesh;


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
  apf::writeVtkFiles("adapt_before",mesh);


  SetSizeField sf(mesh);
  ma::Input* in = ma::configure(mesh, &sf,0,logInterpolation);
 // ma::Input* in = ma::configureUniformRefine (mesh, 1);

  std::cout << "cansnap() : " << mesh->canSnap() << "\n";
  in->shouldSnap = false;
  in->shouldTransferParametric = false;  
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;  
  in->goodQuality = 0.2;

  ma::adapt(in);
  mesh->verify();
  apf::writeVtkFiles("adapt_after",mesh);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return 0;
}



