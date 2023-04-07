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


void setSizeFieldOnVertex(int indx, double& xSize, double& ySize, double* dirVector, double* meshParameters);
void doubleSizeField( int refIndx, int indx, double& xSize, double& ySize, double* dirVector, int* field_id1, int* field_id2, double factor);


// Anisotropic Mesh Size Field
int find_anisofield (ma::Mesh* m, gmi_model* model, int* field_id1, int* field_id2, double* dir)
{
	ma::Vector low;
        ma::Vector up;

        double average = ma::getAverageEdgeLength(m);
	ma::getBoundingBox(m,low,up);
        double lower = low[0];
        double upper = up[0];
 	
	double data_1, data_2;
	int size_data=1;

	apf::MeshTag* sizeFieldSet = m->createIntTag("sfDefined", 1);
	double meshParameters[3] = {average, lower, upper};
	std::vector <int> adaptModelFaces;
	adaptModelFaces.push_back(1);		// Insert the Ids of model faces to adapt
	for (int nFace = 0; nFace < adaptModelFaces.size(); ++nFace)
	{
		int modelFaceIndx = adaptModelFaces[nFace];
		int iterCount = 0;
		for (int i=0; i<m->count(0); ++i)
		{
			apf::MeshEntity* mV = getMdsEntity(m, 0, i);
			gmi_ent* gent= (gmi_ent*)(m->toModel(mV));
			int gTag = gmi_tag(m3dc1_model::instance()->model, gent);
			int gType = gmi_dim(m3dc1_model::instance()->model,gent);
			int vTag = 0;
			if (gType == 2 && gTag == modelFaceIndx)
			{
				setSizeFieldOnVertex(i,data_1, data_2, dir, meshParameters);
				m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
				m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);
				vTag = 1;
				m->setIntTag(mV, sizeFieldSet, &vTag);				
			}
			else
				m->setIntTag(mV, sizeFieldSet, &vTag);
		}
	
		// In Second Step look at the boundaries (model edges bounding the face). 
		// For now we dont have full model adjacency so hard code it.
		std::vector <int> adaptModelEdges;
		adaptModelEdges.push_back(1);			// Model Adjacencies will provide this information.
		for (int nEdge = 0; nEdge < adaptModelEdges.size(); ++nEdge)
		{
			int modelEdgeIndx =  adaptModelEdges[nEdge];
			std::vector <apf::MeshEntity*> meshV, meshVTemp;
			int modelNodeIndx = 1;
			for (int i=0; i<m->count(0); ++i)
			{
				apf::MeshEntity* mV = getMdsEntity(m, 0, i);
				gmi_ent* gent= (gmi_ent*)(m->toModel(mV));
				int gTag = gmi_tag(m3dc1_model::instance()->model, gent);
				int gType = gmi_dim(m3dc1_model::instance()->model,gent);
				if ((gType == 1 && gTag == modelEdgeIndx) || (gType == 0 && gTag == modelNodeIndx))
				{
					setSizeFieldOnVertex(i,data_1, data_2, dir, meshParameters);
					m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
					m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);
					m->removeTag(mV, sizeFieldSet);
					int vTag = 1;
					m->setIntTag(mV, sizeFieldSet, &vTag);
					meshV.push_back(mV);
				}
			}	
		
			while (meshV.size() > 0)
			{
				apf::MeshEntity* mV = meshV[0];
				double d_1, d_2;
				int vId= getMdsIndex(m, mV);
				m3dc1_node_getfield(&vId, field_id1, &d_1, &size_data);
				m3dc1_node_getfield(&vId, field_id2, &d_2, &size_data);
				apf::Adjacent edges;
				m->getAdjacent(mV,1, edges);
				for (int i = 0; i < edges.getSize(); ++i)
				{
					apf::MeshEntity* edge = edges[i];
					gmi_ent* gent= (gmi_ent*)(m->toModel(edge));
                                        int gTag = gmi_tag(m3dc1_model::instance()->model, gent);
                                        int gType = gmi_dim(m3dc1_model::instance()->model,gent);
                                        if ((gType == 2 && gTag == modelFaceIndx) || (gType == 1 && gTag == modelEdgeIndx))
                                                continue;
					else
					{
						apf::MeshEntity* otherV = getEdgeVertOppositeVert(m, edge, mV);
						int vIndx = getMdsIndex(m, otherV);
						int vTagSet = -1;
						m->getIntTag(otherV, sizeFieldSet, &vTagSet);
						if (vTagSet == 1)
							continue;
						double factor = 2;
						doubleSizeField(vId, vIndx, data_1, data_2, dir, field_id1, field_id2,factor);
						m3dc1_node_setfield(&vIndx, field_id1, &data_1, &size_data);
						m3dc1_node_setfield(&vIndx, field_id2, &data_2, &size_data);
						m->removeTag(otherV, sizeFieldSet);	
						int vTag = 1;
						m->setIntTag(otherV, sizeFieldSet, &vTag);
						meshVTemp.push_back(otherV);				
					}
				}
				meshV.erase(meshV.begin());
				if (meshV.size() == 0 && iterCount < 2)
				{
					for (int i = 0; i < meshVTemp.size(); ++i)
						meshV.push_back(meshVTemp[i]);
					meshVTemp.clear();
					iterCount++;
				}		

			} 
		}
	
	}

	for (int i=0; i<m->count(0); ++i)
        {
		apf::MeshEntity* mV = getMdsEntity(m, 0, i);
		int vTagSet = -1;
		m->getIntTag(mV, sizeFieldSet, &vTagSet);

		if (vTagSet == 1)
			continue;
				
		setSizeFieldOnVertex(i,data_1, data_2, dir, meshParameters);
		m3dc1_node_setfield(&i, field_id1, &data_1, &size_data);
    		m3dc1_node_setfield(&i, field_id2, &data_2, &size_data);		
        }

  	m3dc1_field_sync(field_id1);
    	m3dc1_field_sync(field_id2);
}

// A function to set field given the index of vertex
void setSizeFieldOnVertex(int indx, double& xSize, double& ySize, double* dirVector, double* meshParameters)
{
	double average = meshParameters[0]; 
	double lower = meshParameters[1];
	double upper = meshParameters[2];
	double coord[3];
	m3dc1_node_getcoord(&indx, coord);
	double pos = coord[0];
        double x = (pos - lower)/(upper - lower);
        double sizeFactor = average*(4*x+2)/3;

	if (x < 0.5)
		sizeFactor = 2;
	if (x >= 0.5 && x < 0.8)
		sizeFactor = 3;

	xSize = (average/sizeFactor)*0.5;
	ySize = (average/sizeFactor)*0.5;

	dirVector[indx*3+0] = 1.0;
	dirVector[indx*3+1] = 0.0;
	dirVector[indx*3+2] = 0.0;
}

void doubleSizeField( int refIndx, int indx, double& xSize, double& ySize, double* dirVector, int* field_id1, int* field_id2, double factor)
{
	double data1, data2;
	int size_data=1;
	m3dc1_node_getfield(&refIndx, field_id1, &data1, &size_data);
	m3dc1_node_getfield(&refIndx, field_id2, &data2, &size_data);
	
	xSize = data1*factor;
	ySize = data2*factor;

        dirVector[indx*3+0] = 1.0;
        dirVector[indx*3+1] = 0.0;
        dirVector[indx*3+2] = 0.0;
	
}


int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  lion_set_verbosity(1);  

  if (argc<3 & !pumi_rank())
  {
    cout<<"Usage: ./main  model mesh #planes model/mesh_options (0-2) \n"
        <<"\tOption 0 (default): load .txt model and 2D mesh (3D mesh is constructed if #planes>1)\n"
        <<"\tOption 1: load .txt model and distributed 3D mesh (no 3D construction needed)\n"
        <<"\tOption 2: load .dmg model and 2D mesh\n";
    return M3DC1_FAILURE;
  }
  
  int zero=0, num_plane;
  gmi_model* g;
  apf::Mesh2* m;
  
  if (argc==3 || atoi(argv[3])==-1 || (argc>4 && atoi(argv[4])==2))
  {
    gmi_register_mesh();
    g = gmi_load(argv[1]);
    m3dc1_mesh::instance()->mesh = apf::loadMdsMesh(g, argv[2]);
    m3dc1_mesh::instance()->initialize();
  }
  else
  {
    if (m3dc1_model_load(argv[1])) // model loading failed
    {
      PetscFinalize();
      m3dc1_scorec_finalize();
      MPI_Finalize();
      return 0;
    }
    num_plane = atoi(argv[3]);
    if (num_plane>1 && pumi_size()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);
  
    // loading m3dc1 model and mesh directly -- no 3d mesh buildup 
    if (argc>4 && atoi(argv[4])==1)
      m3dc1_mesh_load_3d(argv[2], &num_plane);
    else
    {
      if (m3dc1_mesh_load(argv[2]))  // mesh loading failed
      {
        PetscFinalize();
        m3dc1_scorec_finalize();
        MPI_Finalize();
        return 0;
      }

      if (num_plane>1)
        m3dc1_mesh_build3d(&zero, &zero, &zero);
    }
    m3dc1_field_importall();
  }

  // print model face adjacency;
  if (!PCU_Comm_Self())
  {
    gmi_iter* it = gmi_begin(m3dc1_model::instance()->model, 2);
    gmi_ent* gent;
    int gTag, gDim;
    vector<int> adj_ids;

    while ((gent = gmi_next(m3dc1_model::instance()->model, it))) 
    {
      std::cout<<"Model face "<<gTag<<" adjacency info\n";
      gTag = gmi_tag(m3dc1_model::instance()->model, gent);
      gDim = gmi_dim(m3dc1_model::instance()->model, gent);
      adj_ids.clear();
      get_gent_adj(2, gTag, 1, adj_ids); 
      std::cout<<"\t adj edge: ";
      for (vector<int>::iterator vit=adj_ids.begin(); vit!=adj_ids.end(); ++vit)
        std::cout<<*vit<<" ";
      std::cout<<"\n";
      adj_ids.clear();
      std::cout<<"\t adj vertex: ";
       get_gent_adj(2, gTag, 0, adj_ids);  
      for (vector<int>::iterator vit=adj_ids.begin(); vit!=adj_ids.end(); ++vit)
        std::cout<<*vit<<" ";
      std::cout<<"\n";
    }
    gmi_end(m3dc1_model::instance()->model, it);
  }

  m = m3dc1_mesh::instance()->mesh; 
  int num_adapt_iter = 1; 
  if (argc>5 && atoi(argv[4])>0) 
    num_adapt_iter=atoi(argv[4]);

  int fid_size1, fid_size2;

  for (int n = 0; n<num_adapt_iter; ++n)
  {
    m3dc1_field_getnewid (&fid_size1);
    fid_size2 = fid_size1+1;
    int scalar_type=0;
    int num_value=1;

    m3dc1_field_create (&fid_size1, "size 1", &num_value, &scalar_type, &num_value);
    m3dc1_field_create (&fid_size2, "size 2", &num_value, &scalar_type, &num_value);
   
    double* dir = new double[m->count(0)*3];

    // Set Size Field	
    find_anisofield (m, m3dc1_model::instance()->model, &fid_size1, &fid_size2, dir);   
  
    // Identify the elements where no adapt is required
    apf::MeshTag* doNotAdaptFace = m->createIntTag("doNotAdapt", 1);
    for (int i=0; i<m->count(2); ++i)
    {
       apf::MeshEntity* ele = getMdsEntity(m, 2, i);
       gmi_ent* gent= (gmi_ent*)(m->toModel(ele));
       int gTag = gmi_tag(m3dc1_model::instance()->model, gent);
       int gType = gmi_dim(m3dc1_model::instance()->model,gent);
       int elemTag = -1;
       apf::MeshEntity*  vertices[3];
       m->getDownward(ele, 0, vertices);
       int count = 0;
       apf::MeshTag* vTag = m->findTag("sfDefined");
       for (int j =0; j < 3; ++j)
       {
               int tagNode = -1;
               m->getIntTag(vertices[j], vTag, &tagNode);
               if (tagNode == 1)
                       count++;
       }

        if (count == 3)
                elemTag = 0;
        else
                elemTag = 1;

        m->setIntTag(ele, doNotAdaptFace, &elemTag);
    }
    apf::writeVtkFiles("before-adapt", m); 
    m3dc1_mesh_adapt(&fid_size1, &fid_size2, dir);
    apf::writeVtkFiles("adapt_after",m);
    if (!PCU_Comm_Self()) std::cout << "adaptation completed\n";
  
    delete [] dir;
  } // adaptive loop

  apf::printStats(m3dc1_mesh::instance()->mesh); 
  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();

  return 0;
}
