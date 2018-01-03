/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef NOM3DC1
#define NOM3DC1
#endif
#include "m3dc1_scorec.h"
#include <mpi.h>
#include "PCU.h"
#include "m3dc1_mesh.h"
#include <iostream>
#include <assert.h>

using namespace std;

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff);

int main( int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  if (argc<4 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes"<<endl;
    return M3DC1_FAILURE;
  }

  m3dc1_model_load(argv[1]);
  int  num_plane =  atoi(argv[3]);
  m3dc1_model_setnumplane (&num_plane);
  if (!PCU_Comm_Self()) std::cout<<" num plane "<<num_plane<<endl;
  m3dc1_mesh_load(argv[2]);
  // set/get field dof values
  int num_vertex, vertex_dim=0;

  int zero=0;
  m3dc1_mesh_build3d(&zero, &zero, &zero);
  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);
  int field=1;
  int num_values=1, num_dofs=2, value_type[]={0};
  int num_dofs_node = num_values * num_dofs;
  m3dc1_field_create (&field, "test_field", &num_values, value_type, &num_dofs);
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &field, &num_dofs_node, xyz);
  }
  m3dc1_field_sum_plane(&field);
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    double dofs[2];
    int num_dofs_t;
    m3dc1_ent_getdofdata(&vertex_dim, &inode, &field, &num_dofs_t, dofs);
    for(int i=0; i<2; i++)
      assert(AlmostEqualDoubles(dofs[i], xyz[i]*num_plane, 1e-12, 1e-12));
  }

  m3dc1_field_delete(&field);
  m3dc1_scorec_finalize();
  MPI_Finalize();
}

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff)
{
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

