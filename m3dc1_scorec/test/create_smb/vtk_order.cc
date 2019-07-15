/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <deque>
#include <set>
#include <PCU.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <apfVector.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi.h>
#include <gmi_null.h>

// seol -- reorder input mesh before conversion
//         starting vtx: a vtx with min Y
struct Queue {
  bool has(apf::MeshEntity* e) { return h.count(e); }
  void push(apf::MeshEntity* e)
  {
    q.push_back(e);
    h.insert(e);
  }
  void pushVector(std::vector<apf::MeshEntity*> const& l)
  {
    for (size_t i = 0; i < l.size(); ++i)
      push(l[i]);
  }
  apf::MeshEntity* pop()
  {
    apf::MeshEntity* e;
    e = q.front();
    q.pop_front();
    h.erase(e);
    return e;
  }
  bool empty() { return q.empty(); }
  std::deque<apf::MeshEntity*> q;
  std::set<apf::MeshEntity*> h;
};


apf::MeshEntity* findFirst(apf::Mesh* m)
{
  apf::MeshEntity* v;
  apf::MeshEntity* best;
  apf::MeshIterator* it = m->begin(0);
  best = m->iterate(it);
  apf::Vector3 coord;
  m->getPoint(best, 0, coord);
  double min_Y=coord[1];
  while ((v = m->iterate(it)))
  {
    m->getPoint(v, 0, coord);  
    if (min_Y > coord[1])
    {
      best = v;
      min_Y = coord[1];
    }
  }
  m->end(it);
  return best;
}

bool visited(Queue& q, apf::Numbering* nn, apf::MeshEntity* e)
{
  return apf::isNumbered(nn, e, 0, 0) || q.has(e);
}

bool hasNode(apf::Mesh* m, apf::MeshEntity* e)
{
  if (m->getShape()->countNodesOn(m->getType(e))>0) 
    return true;
  return false;
}

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  if (argc!=2)
  {
    std::cout<<"SCOREC ERROR: wrong input arguments. \n"
             <<"  Usage: ./reorder inmesh.smb\n";
    PCU_Comm_Free();
    MPI_Finalize();    
    return 1;
  }
  gmi_register_null();
  gmi_model* nullModel = gmi_load(".null");

  apf::Mesh2* m = apf::loadMdsMesh(nullModel,argv[1]);
  int mesh_dim = m->getDimension();

  apf::Numbering* n_org=apf::numberOwnedDimension(m, "elem", mesh_dim);
  writeVtkFiles("org-mesh.vtk", m);
  std::cout << "[INFO] vtk files for original ordering saved in \"org-mesh.vtk\"\n";
  destroyNumbering(n_org);

  // reorder mesh
  apf::Numbering* nn = apf::createNumbering(m, "node", apf::getConstant(0), 1);
  apf::Numbering* en = apf::createNumbering(m, "elem", apf::getConstant(mesh_dim), 1);

  Queue q;
  q.push(findFirst(m));

  // node and element number starts from 0
  int labelnode = 0;
  int labelelem = 0;

  std::vector<apf::MeshEntity*> node_arr;
  std::vector<apf::MeshEntity*> elem_arr;

  node_arr.resize(m->count(0)+1);
  elem_arr.resize(m->count(mesh_dim)+1);

  apf::MeshEntity* otherVtx;
  apf::MeshEntity* edge;
  apf::MeshEntity* elem;

  while (!q.empty()) 
  {
    apf::MeshEntity* vtx = q.pop();
    if (!apf::isNumbered(nn, vtx, 0, 0))
    {
      node_arr[labelnode] = vtx;
      apf::number(nn, vtx, 0, 0, labelnode);

      ++labelnode;
    }

    std::vector<apf::MeshEntity*> entities;
    apf::Adjacent edges;
    m->getAdjacent(vtx,1, edges);
    for (size_t i = 0; i < edges.getSize(); ++i) 
    {
      edge = edges[i];
      apf::Adjacent adjacent;
      m->getAdjacent(edge, mesh_dim, adjacent);      
      for (size_t j = 0; j < adjacent.getSize(); ++j) 
      {
        elem = adjacent[j];
        if (!apf::isNumbered(en, elem, 0, 0))
        {
          elem_arr[labelelem] = elem;
          apf::number(en, elem, 0, 0, labelelem);
          ++labelelem;
        }
      }
      otherVtx = apf::getEdgeVertOppositeVert(m, edge, vtx);
      if (!visited(q, nn, otherVtx))
        entities.push_back(otherVtx);
    }
    q.pushVector(entities);
  } // while
  destroyNumbering(nn);

  writeVtkFiles("reordered-mesh.vtk", m);
  std::cout << "[INFO] vtk files for adjacency-based ordering saved in \"reordered-mesh.vtk\"\n";
  destroyNumbering(en);

  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
