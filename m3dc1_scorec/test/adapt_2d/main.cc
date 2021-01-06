#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "apf.h"
#include <apfMDS.h>
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


int find_isofield (ma::Mesh* m, std::vector <double>& size_h1,std::vector <double>& size_h2,std::vector <double>& dir_1);
int find_anisofield (ma::Mesh* m, std::vector <double>& size_h1,std::vector <double>& size_h2,std::vector <double>& dir_1);


// Function to create Isotropic Size Field
int find_isofield (ma::Mesh* m, std::vector <double>& size_h1,std::vector <double>& size_h2,std::vector <double>& dir_1)
{
	
        ma::Mesh* mesh;
        ma::Vector low;
        ma::Vector up;

	double average = ma::getAverageEdgeLength(m);	
        ma::getBoundingBox(m,low,up);
        
        double lower = low[0];
        double upper = up[0];

 	int numVert=m3dc1_mesh::instance()->mesh->count(0);
	std::cout << " Number of vertices are: " << numVert << "\n";
	for (int i=0; i<numVert; i++)
	{
	//	apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
		double coord[3];
		m3dc1_node_getcoord(&i, coord);
		double pos = coord[0];	
        	double sizeFactor = 2;

	//	std::cout << average/sizeFactor << "\n";
        	size_h1.push_back(average/sizeFactor);
        	size_h2.push_back(average/sizeFactor);
        	dir_1.push_back(1.0);
		dir_1.push_back(0.0);
		dir_1.push_back(0.0);
	}
        return M3DC1_SUCCESS;
}

//// Function to create AnIsotropic Size Field
int find_anisofield (ma::Mesh* m, std::vector <double>& size_h1,std::vector <double>& size_h2,std::vector <double>& dir_1)
{

        ma::Mesh* mesh;
        ma::Vector low;
        ma::Vector up;

        double average = ma::getAverageEdgeLength(m);
        ma::getBoundingBox(m,low,up);

        double lower = low[0];
        double upper = up[0];

        int numVert=m3dc1_mesh::instance()->mesh->count(0);
        std::cout << " Number of vertices are: " << numVert << "\n";
        for (int i=0; i<numVert; i++)
        {
                double coord[3];
                m3dc1_node_getcoord(&i, coord);
                double pos = coord[0];

                double x = (pos - lower)/(upper - lower);
                double sizeFactor = 2;
		
		if (x < 0.5)
                	sizeFactor = 5;
        	if (x >= 0.5 && x < 0.8)
                	sizeFactor = 0.5;
			
                size_h1.push_back(average);
                size_h2.push_back(average/sizeFactor);
                dir_1.push_back(1.0);
                dir_1.push_back(0.0);
                dir_1.push_back(0.0);
        }
        return M3DC1_SUCCESS;
}


class SetIsoSizeField : public ma::AnisotropicFunction
{
  public:
    SetIsoSizeField(ma::Mesh* m,std::vector <double>& size_h1,std::vector <double>& size_h2,std::vector <double>& dir_1 )
    {
	mesh = m;
	int arr_size = size_h1.size();
	double size_1[arr_size], size_2[arr_size];
        double angle[3*arr_size];
	for (int i = 0; i<arr_size; ++i)
	{
		size_1[i] = size_h1[i];
                size_2[i] = size_h2[i];
                angle[(i*3)] = dir_1[(i*3)];
		angle[(i*3)+1] = dir_1[(i*3)+1];
		angle[(i*3)+2] = dir_1[(i*3)+2];
	}	

    }

    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {

        double h1, h2;
        double angle_1[3];

        for (int i = 0; i<arr_size; ++i)
        {
                h1 = size_1[i];
                h2 = size_2[i];
                
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

		ma::Vector h(h1, h2,h2);

		R[0][0]=angle_1[0];
  		R[0][1]=angle_1[1];
  		R[0][2]=0.0;

  		R[1][0]= dir_2[0];
  		R[1][1]= dir_2[1];
  		R[1][2]=0.0;

  		R[2][0]=0;
  		R[2][1]=0;
  		R[2][2]=1.;

		H = h;
	}	
    }
  private:
        ma::Mesh* mesh;
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
  

  std::vector <double> size_h1;
  std::vector <double> size_h2;
  std::vector <double> dir_1;

  int field_success = find_isofield (mesh, size_h1,size_h2, dir_1);
  
  std::cout << "Debug-A" << "\n";
  SetIsoSizeField sf(mesh, size_h1, size_h2, dir_1);
  ma::Input* in = ma::configure(mesh, &sf,0,logInterpolation);

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



