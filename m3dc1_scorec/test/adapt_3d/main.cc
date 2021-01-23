/****************************************************************************** 

  (c) 2005-2021 Scientific Computation Research Center, 
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
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "petscksp.h"
#include "PCU.h"


// Isotropic  Mesh Size Field
void find_isofield (ma::Mesh* m, int* field_id1, int* field_id2, double* dir)
{
  double average = ma::getAverageEdgeLength(m);

  double data = average/2;

  int size_data=1;
  for (int i=0; i<m->count(0); ++i)
  {
    m3dc1_node_setfield(&i, field_id1, &data, &size_data);
    m3dc1_node_setfield(&i, field_id2, &data, &size_data);

    dir[i*3+0] = 1.0;
    dir[i*3+1] = 0.0;
    dir[i*3+2] = 0.0;
  }
  m3dc1_field_sync(field_id1);
  m3dc1_field_sync(field_id2);
}

// Anisotropic Mesh Size Field
int find_anisofield (ma::Mesh* m, int* field_id1, int* field_id2, double* dir)
{
	ma::Vector low;
        ma::Vector up;

        double average = ma::getAverageEdgeLength(m);
	ma::getBoundingBox(m,low,up);
        double lower = low[0];
        double upper = up[0];
 
	
	double data_1, data_2;
	int size_data=1;
    
  
	for (int i=0; i<m->count(0); ++i)
        {
                double coord[3];
                m3dc1_node_getcoord(&i, coord);
                double pos = coord[0];

                double x = (pos - lower)/(upper - lower);
                double sizeFactor = average*(4*x+2)/3; 


                if (x < 0.5)
                        sizeFactor = 3;
                if (x >= 0.5 && x < 0.8)
                        sizeFactor = 2;
		
		
		data_1 = sizeFactor;		
		data_2 = average/sizeFactor;
		
		m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
    		m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);		

    		dir[i*3+0] = 1.0;
    		dir[i*3+1] = 0.0;
    		dir[i*3+2] = 0.0;
        }
  	m3dc1_field_sync(field_id1);
    	m3dc1_field_sync(field_id2);
	
}

// Anisotropic Mesh Size Field- Shock Function
int find_anisofield_shock (ma::Mesh* m, int* field_id1, int* field_id2, double* dir)
{
        ma::Vector low;
        ma::Vector up;

        double average = ma::getAverageEdgeLength(m);
        ma::getBoundingBox(m,low,up);
        double lower = low[0];
        double upper = up[0];

	double shock_cx = (low[0]+up[0])/2;             // Center of Shock x component
        double shock_cy = (low[1]+up[1])/2;             // Center of shock y component
        double shock_radius = 0.5;
        double shock_thickness = 0.1;

        double min_size = 0.001;
        double max_size = 0.1;
	
        double data_1, data_2;
        int size_data=1;
	

        for (int i=0; i<m->count(0); ++i)
        {
		// Find the coordinates of the vertices
                double coord[3];
                m3dc1_node_getcoord(&i, coord);

		
                double pos[2];
		pos[0] = coord[0];
		pos[1] = coord[1];

                double x = pos[0] - shock_cx;
                double y = pos[1] - shock_cy;
 		double r2 = x*x + y*y;
		double r = sqrt(r2);
		
		if (abs(r - shock_radius) > shock_thickness) 			// outside shock
		{
			data_1 = max_size;
			data_2 = max_size;
                	dir[i*3+0] = 1.0;
                	dir[i*3+1] = 0.0;
                	dir[i*3+2] = 0.0;
		}		
		else								// Inside Shock
		{
			double size = min_size + (max_size - min_size) * abs(r - shock_radius)/shock_thickness;
			data_1 = size;
			data_2 = max_size;
			dir[i*3+0] = x/r;
                        dir[i*3+1] = y/r;
                        dir[i*3+2] = 0.0;			
		}
	

                m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
                m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);

        }
          m3dc1_field_sync(field_id1);
          m3dc1_field_sync(field_id2);

}

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  lion_set_verbosity(1);  

  if (argc<3 & !pumi_rank())
  {
    cout<<"Usage: ./main  model mesh #planes"<<endl;
    return M3DC1_FAILURE;
  }
  
  int num_plane=pumi_size();
  gmi_model* g;
  apf::Mesh2* m;
  
  if (!atoi(argv[3]))
  {
    gmi_register_mesh();
    g = gmi_load(argv[1]);
    m3dc1_mesh::instance()->mesh = apf::loadMdsMesh(g, argv[2]);
    m3dc1_mesh::instance()->initialize();
  }
  else
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && pumi_size()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);

    if (m3dc1_model_load(argv[1])) // model loading failed
    {
      PetscFinalize();
      m3dc1_scorec_finalize();
      MPI_Finalize();
      return 0;
    }

    if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
    {
      PetscFinalize();
      m3dc1_scorec_finalize();
      MPI_Finalize();
      return 0;
    }
    int zero=0;

    if (num_plane>1)
      m3dc1_mesh_build3d(&zero, &zero, &zero);
  }

  m = m3dc1_mesh::instance()->mesh; 
  int num_adapt_iter = 1; 
  if (argc>5) 
 	num_adapt_iter=atoi(argv[4]);
  if (!PCU_Comm_Self()) std::cout << "size of adapt loops = " << num_adapt_iter << "\n";
  int shouldSnap=0;
  if (argc>6) 
    shouldSnap=atoi(argv[5]);
  else
    if (!PCU_Comm_Self()) std::cout <<"Missing input args: argv[5] shouldSnap, argv[6] shouldRunPreZoltan, "
         <<"argv[7] shouldRunPostZoltan, argv[8] shouldefineLayer\n";

  int shouldRunPreZoltan=1;
  int shouldRunPostZoltan=1;
  int shouldRefineLayer=1;
  int maximumIterations=5;
  double goodQuality =0.4;

  if (argc>7) shouldRunPreZoltan=atoi(argv[6]);
  if (argc>8) shouldRunPostZoltan=atoi(argv[7]);
  if (argc>9) shouldRefineLayer=atoi(argv[8]);
  if (argc>10) maximumIterations=atoi(argv[9]);
  if (argc>11) goodQuality =atof(argv[10]);
    
 // apf::writeVtkFiles("before-adapt",m3dc1_mesh::instance()->mesh);

  for (int n = 0; n<num_adapt_iter; ++n)
  {
  	int fid_size1=1;
  	int fid_size2=2;
  	int scalar_type=0;
  	int num_value=1;

  	m3dc1_field_create (&fid_size1, "size 1", &num_value, &scalar_type, &num_value);
  	m3dc1_field_create (&fid_size2, "size 2", &num_value, &scalar_type, &num_value);
   
   	double* dir = new double[m->count(0)*3];

  //    Select one size field from the below provided three options
	
  //	find_isofield (m, &fid_size1, &fid_size2, dir);
  //	find_anisofield (m, &fid_size1, &fid_size2, dir);   
   	find_anisofield_shock (m, &fid_size1, &fid_size2, dir);

  	if (!PCU_Comm_Self()) std::cout << "start adaptation with # max iterations "
                                        <<maximumIterations<<", goodQuality "<<goodQuality<<"\n";
  
  	m3dc1_mesh_adapt(&fid_size1, &fid_size2, dir, &shouldSnap, &shouldRunPreZoltan,
		&shouldRunPostZoltan, &shouldRefineLayer, &maximumIterations, &goodQuality);
  

  	pumi_mesh_print(m);

  	if (!PCU_Comm_Self()) std::cout << "adaptation completed\n";
  
  	delete [] dir;
  }
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
