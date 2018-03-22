/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_SLNTRANSFER_H
#define M3DC1_SLNTRANSFER_H
#include "maSolutionTransfer.h"
#include <apf.h>
#include <apfVector.h>
#include <apfNumbering.h>
#include <iostream>
#include <vector>
#include <ReducedQuintic.h>
using namespace std;

class ReducedQuinticTransfer : public ma::SolutionTransfer
{
  public:
    ReducedQuinticTransfer(apf::Mesh* m, vector<apf::Field*>& fds, ReducedQuintic* shape)
      : fields(fds), mesh(m), thecase(shape)
    {
      int maxComp=0;
      for(size_t ii=0; ii < fds.size(); ii++)
        maxComp=max(maxComp,apf::countComponents(fds.at(ii)));
      value.allocate(6*maxComp);
    }
    virtual void onVertex(
        apf:: MeshElement* parent,
        ma:: Vector const& xi,
        ma:: Entity* vert);
    virtual bool hasNodesOn(int dimension)
    {
      return dimension==0;
    }
  protected:
    vector<apf::Field*> fields;
    apf::Mesh* mesh;
    apf::NewArray<double> value;
    ReducedQuintic* thecase;
    static int dofNode;
};
#endif
