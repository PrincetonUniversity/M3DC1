/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "slntransferUtil.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <m3dc1_scorec.h>
#include "PCU.h"
// return origin para
void reorder( double vtxCoord[3][2], int order[3]);
void evalGeoPara(double vtxCoord[3][2], double origin[2], double para[5],int order[3] )
{
  // first get the order
  // set the longest edge as xi axis
  reorder( vtxCoord, order); 
  double abcCoord[3][2];
  for(int i=0;i<3;i++)
  {
    int idx=order[i];
    abcCoord[i][0]=vtxCoord[idx][0];
    abcCoord[i][1]=vtxCoord[idx][1];
  }
  double ab[2]={abcCoord[1][0]-abcCoord[0][0],abcCoord[1][1]-abcCoord[0][1]};
  double ablensq=ab[0]*ab[0]+ab[1]*ab[1];
  double ac[2]={abcCoord[2][0]-abcCoord[0][0],abcCoord[2][1]-abcCoord[0][1]};
  double buffer=(ab[0]*ac[0]+ab[1]*ac[1])/ablensq;
  double ac2ab[2]={ab[0]*buffer,ab[1]*buffer};
  origin[0]=ac2ab[0]+abcCoord[0][0];
  origin[1]=ac2ab[1]+abcCoord[0][1];
  double a[2]={origin[0]-abcCoord[0][0],origin[1]-abcCoord[0][1]};
  double b[2]={origin[0]-abcCoord[1][0],origin[1]-abcCoord[1][1]}; 
  double c[2]={origin[0]-abcCoord[2][0],origin[1]-abcCoord[2][1]};
  para[1]=sqrt(a[0]*a[0]+a[1]*a[1]);
  para[0]=sqrt(b[0]*b[0]+b[1]*b[1]);
  para[2]=sqrt(c[0]*c[0]+c[1]*c[1]);
  assert(fabs(para[0]+para[1]-sqrt(ablensq))<tol*ablensq);
  para[3]=ab[1]/sqrt(ablensq);
  para[4]=ab[0]/sqrt(ablensq);
}

void reorder( double vtxCoord[3][2], int order[3])
{
  double v1v2[2]={vtxCoord[1][0]-vtxCoord[0][0],vtxCoord[1][1]-vtxCoord[0][1]};
  double v1v3[2]={vtxCoord[2][0]-vtxCoord[0][0],vtxCoord[2][1]-vtxCoord[0][1]};
  double v2v3[2]={vtxCoord[2][0]-vtxCoord[1][0],vtxCoord[2][1]-vtxCoord[1][1]};
  int counterclockwise=1;
  if(v1v2[0]*v1v3[1]-v1v2[1]*v1v3[0]<0) counterclockwise=0;
  double v1v2lensq=v1v2[0]*v1v2[0]+v1v2[1]*v1v2[1];
  double v1v3lensq=v1v3[0]*v1v3[0]+v1v3[1]*v1v3[1];
  double v2v3lensq=v2v3[0]*v2v3[0]+v2v3[1]*v2v3[1];
  if(v1v2lensq>v1v3lensq&&v1v2lensq>v2v3lensq)
  {
    if(counterclockwise)
    {
      order[0]=0;
      order[1]=1;
      order[2]=2;
    }
    else
    {
      order[0]=1;
      order[1]=0;
      order[2]=2;
    }
  }
  else
  {
    if(v1v3lensq>v2v3lensq)
    {
      if(counterclockwise)
      {
        order[0]=2;
        order[1]=0;
        order[2]=1;
      }
      else
      {
        order[0]=0;
        order[1]=2;
        order[2]=1;
      }
    }
    else
    {
      if(counterclockwise)
      {
        order[0]=1;
        order[1]=2;
        order[2]=0;
      }
      else
      {
        order[0]=2;
        order[1]=1;
        order[2]=0;
      }
    }
  }
}

void printmatrix(double* A, int M, int N, char* filename)
{
  FILE* fp= fopen (filename, "w");
  for ( int i=0; i<M; i++)
  {
    for( int j=0; j<N; j++)
     fprintf(fp, "%.16e ",A[i*M+j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void getTriangle(double v1[2], double v2[2], double newv[2], double & a, double origin[2], double theta[2])
{
  origin[0]=0.5*(v2[0]+v1[0]);
  origin[1]=0.5*(v2[1]+v1[1]);
  double AO[2]={origin[0]-v1[0],origin[1]-v1[1]};
  a=sqrt(AO[0]*AO[0]+AO[1]*AO[1]);
  // rotate AO by pi/2 counterclock wise
  double OC[2]={-AO[1],AO[0]};
  newv[0]=origin[0]+OC[0];
  newv[1]=origin[1]+OC[1];
  theta[0]=AO[1]/a;
  theta[1]=AO[0]/a;
}

inline int checkMagnitude(double a);
int containPoint(double vtxcoord[3][2], double pt[2])
{
  double x1=vtxcoord[0][0];
  double y1=vtxcoord[0][1];
  double x2=vtxcoord[1][0];
  double y2=vtxcoord[1][1];
  double x3=vtxcoord[2][0];
  double y3=vtxcoord[2][1];
  double x=pt[0];
  double y=pt[1];
  double detT=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
  double lamda_1=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
  double lamda_2=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
  double lamda_3=1.0-lamda_1-lamda_2;
  return checkMagnitude(lamda_1)&&checkMagnitude(lamda_2)&&checkMagnitude(lamda_3);
}
//check if a is [0,1]
inline int checkMagnitude (double a)
{
  return -tol<a&& a<1+tol;
}

void output_face_data (int * size, double * data, char * vtkfile)
{
  FILE *fp;
	
  char filename_buff[256];
  strcpy(filename_buff,vtkfile);
#ifdef USEVTK
  sprintf(filename_buff, "%s_%d.vtk",filename_buff,PCU_Comm_Self());
  /** write out the mesh file first */
  fp = fopen(filename_buff, "w");
  int zero=0,numnodes;
  m3dc1_mesh_getnument(&zero, &numnodes);
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "%d\n", 0);
  fprintf(fp, "ASCII\n\n");
  /** write out the points */
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %d float\n", numnodes);

  for(int i=0;i<numnodes;i++)
  {
    double v_xyz[3];
    m3dc1_node_getcoord(&i,v_xyz);
    fprintf(fp, "%.16e %.16e %.16e\n", v_xyz[0], v_xyz[1], v_xyz[2]);
  }
	
  fprintf(fp, "\n");
  /** write out the cells */
  int two=2,numfaces;
  m3dc1_mesh_getnument(&two, &numfaces);
  int numcells=numfaces;
  int numdata=numfaces*3 + numfaces;
  fprintf(fp, "CELLS %d %d\n",numcells , numdata);

  for(int i=0;i<numfaces;i++)
  {	
    fprintf(fp, "%d ", 3);
    int nodes[3];
    int num_get=3, num_get_t;
    m3dc1_ent_getadj (&two, &i, &zero, nodes, &num_get, &num_get_t);
    for(int i=0; i<3; i++)
    {
      fprintf(fp, "%d ", nodes[i]);
    }
 	
   fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  // write out the cell format
  fprintf(fp, "CELL_TYPES %d\n", numcells);
  for(int i=1;i<=numfaces;i++)
  {
    fprintf(fp, "%d\n", 5);
  }
  fprintf(fp, "\n");

  fprintf(fp, "CELL_DATA %d\n", numfaces);
  fprintf(fp, "SCALARS error float %d\n", *size);
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int face=0;face<numfaces;face++)
  {
    for(int i=0; i<*size; i++)
      fprintf(fp, "%e ", data[i*numfaces+face]);
    fprintf(fp, "\n");
  }
  fclose(fp);
#else
  if(!PCU_Comm_Self())
  {
    sprintf(filename_buff, "%s.pvtu",vtkfile);
    fp=fopen(filename_buff, "w");
    fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\">\n");
    fprintf(fp,"<PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf(fp,"<PPoints>\n");
    fprintf(fp,"<PDataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\"/>\n");
    fprintf(fp,"</PPoints>\n");
    fprintf(fp,"<PCellData>\n");
    fprintf(fp,"<PDataArray type=\"Float64\" Name=\"error\" NumberOfComponents=\"%d\" format=\"ascii\"/> \n", *size);
    fprintf(fp,"</PCellData>\n");
    for(int i=0; i<PCU_Comm_Peers(); i++)
    {
      fprintf(fp, "<Piece Source=\"%s%d.vtu\"/>\n",vtkfile, i);
    }
    fprintf(fp, "</PUnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
  }
  sprintf(filename_buff, "%s%d.vtu",vtkfile,PCU_Comm_Self());
  fp=fopen(filename_buff, "w");
  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  int zero=0,numnodes;
  m3dc1_mesh_getnument(&zero, &numnodes); 
  int two=2,numfaces;
  m3dc1_mesh_getnument(&two, &numfaces);
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numnodes, numfaces);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numnodes;i++)
  {
    double v_xyz[3];
    m3dc1_node_getcoord(&i,v_xyz);
    fprintf(fp, "%lf %lf %lf\n", v_xyz[0], v_xyz[1], v_xyz[2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<numfaces;i++)
  {
    int nodes[3];
    int num_get=3, num_get_t;
    m3dc1_ent_getadj (&two, &i, &zero, nodes, &num_get, &num_get_t);
    for(int i=0; i<3; i++)
    {
      fprintf(fp, "%d ", nodes[i]);
    }
   fprintf(fp, "\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for(int i=0;i<numfaces;i++)
  {
    fprintf(fp, "%d\n", 3*(i+1));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<numfaces;i++)
  {
    fprintf(fp, "%d\n",5);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"error\" NumberOfComponents=\"%d\" format=\"ascii\">\n", *size);
  for (int face=0;face<numfaces;face++)
  {
    for(int i=0; i<*size; i++)
      fprintf(fp, "%e ", data[i*numfaces+face]);
    fprintf(fp, "\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
#endif
}
