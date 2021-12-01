/****************************************************************************** 

  (c) 2005-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_slnTransfer.h"
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "m3dc1_field.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include <assert.h>

#include <pcu_util.h>
#include <PCU.h>

int ReducedQuinticTransfer::dofNode = C1TRIDOFNODE;
///////////////////////////////////////////////////////////////////////////////////////////////////////
void  ReducedQuinticTransfer::onVertex(apf::MeshElement* parent, ma::Vector const& xi, ma::Entity* vert)
{
  apf::MeshEntity* oldEdge= apf::getMeshEntity(parent);

#ifdef DEBUG
  assert(apf::getDimension(m3dc1_mesh::instance()->mesh,oldEdge)==1);
#endif

  apf::MeshEntity* vertices[3];
  int num_face=-1;
  apf:: Up facEdg;
  mesh->getUp(oldEdge,facEdg);
  num_face=facEdg.n;

#ifdef DEBUG
  assert(num_face<=2);
#endif

  for(int fi=0; fi<fields.size(); fi++)
  {
    apf::Field* field = fields.at(fi);
    int numComp = apf::countComponents(field);
    vector<vector<double> > dofsVertex(2);
    apf::MeshEntity*  vertices[2];
    mesh->getDownward(oldEdge, 0, vertices);
    apf::Vector3 xyz;
    m3dc1_mesh::instance()->mesh->getPoint(vert, 0, xyz);
    apf::Vector3 xyz2[2];
    for (int i = 0; i < 2; i++) {
      m3dc1_mesh::instance()->mesh->getPoint(vertices[i], 0, xyz2[i]);
    }
    for( int i=0; i<2; i++)
    {
      dofsVertex.at(i).resize(numComp);
      // TODO most likely this can be done by directly inquiring the field at vertex
      // and there is no need for create element
      apf::getComponents(field, vertices[i], 0, &(dofsVertex[i][0]));
      /* apf::Element* vertex = apf::createElement(field,vertices[i]); */
      /* apf::getComponents(vertex,xi,&(dofsVertex[i][0])); */
      /* m3dc1_mesh::instance()->mesh->getPoint(vertices[i], 0, xyz2[i]); */
      /* apf::destroyElement(vertex); */
    }
    double len1= sqrt((xyz2[0][0]-xyz2[1][0])*(xyz2[0][0]-xyz2[1][0])+(xyz2[0][1]-xyz2[1][1])*(xyz2[0][1]-xyz2[1][1]));
    double len2= sqrt((xyz2[0][0]-xyz[0])*(xyz2[0][0]-xyz[0])+(xyz2[0][1]-xyz[1])*(xyz2[0][1]-xyz[1]));
    double wt = len2/len1;
    assert(wt>=0&&wt<=1);
    if(numComp==1) // mesh size field
    {
      double value = wt*dofsVertex[0][0]+(1-wt)*dofsVertex[1][0];
      assert(value>0);
      apf::setComponents(field,vert,0,&value);
      continue;
    }
    vector<double> newdofs(numComp);
    int face_skip=0;
    for( int iface=0; iface<num_face; iface++)
    {
      apf::MeshEntity* element = facEdg.e[iface];
      apf::MeshEntity*  vertices[3];
      mesh->getDownward(element, 0, vertices);
      double coords[3][2];
      for( int i=0; i<3; i++)
      {
        apf::Vector3 xyz;
        m3dc1_mesh::instance()->mesh->getPoint(vertices[i], 0, xyz);
        coords[i][0]=xyz[0];
        coords[i][1]=xyz[1];
        //debug
        //cout<<" vetex "<<i<<" "<<vertices[i]<<" coord" <<cd[0]<<" "<<cd[1]<<endl;
      }
      // set up reduced quintic shape fns in the face
      thecase->setCoord(coords);
      //get dofs of the three vertex 
      int miss_flag = 0;
      for(int i=0;i<3;i++)
      {
        if(!hasEntity(field,vertices[i]))
        {
          miss_flag=1;
          break;
        }
	// TODO most likely this can be done by directly inquiring the field at vertex
	// and there is no need for create element
        apf::getComponents(field, vertices[i], 0, &(value[numComp*i]));
        /* apf::Element* vertex = apf::createElement(field,vertices[i]); */
        /* apf::getComponents(vertex,xi,&(value[numComp*i])); */
	/* apf::destroyElement(vertex); */
      }
      assert(!miss_flag);
      if(miss_flag)
      {
        face_skip++;
        continue;
      }
        // interpolate dofs on the new vertex
      double coord_newv[3];
      apf::Vector3 xyz;
      m3dc1_mesh::instance()->mesh->getPoint(vert, 0, xyz);
      for(int i=0; i<3; i++)
        coord_newv[i]=xyz[i];
      //cout<<" new vtx "<<ent<<" coord "<<coord_newv[0]<<" "<<coord_newv[1]<<endl;
      int numField=numComp/dofNode;
      for( int ifield=0; ifield<numField; ifield++)
      {
        vector<double> dofs_one_field(3*dofNode);
        vector<double> dofs_cacu(dofNode);
        for( int j=0; j<dofNode; j++)
        {
          dofs_one_field.at(j)=value[ifield*dofNode+j];
          dofs_one_field.at(dofNode+j)=value[numComp+ifield*dofNode+j];
          // use trivial dofs for third vertex
          dofs_one_field.at(2*dofNode+j)=value[numComp*2+ifield*dofNode+j];
        }
        thecase->setDofs(&(dofs_one_field[0]));
        thecase->eval_g(coord_newv,&(dofs_cacu[0]));
        int idx_start=ifield*dofNode;
        for(int i=0; i<dofNode; i++)
        {
          newdofs.at(idx_start+i)+=dofs_cacu[i]; 
        }
      }
    }
    if(face_skip<num_face)
    {
      for( int i=0; i<newdofs.size(); i++)
      {
        newdofs.at(i)/=(double)(num_face-face_skip);
      }
    }
    else
    {
      apf::MeshEntity*  vertices[2];
      mesh->getDownward(oldEdge, 0, vertices);
      int numVtx=0;
      vector<double> dofsBuff(numComp,0);  
      for(int i=0; i<2; i++)
      {
        if(hasEntity(field,vertices[i]))
        {
          getComponents(field, vertices[i], 0, &dofsBuff[0]);
          for(int j=0; j<numComp; j++)
            newdofs[j]+=dofsBuff[j];
          numVtx++;
        }
      }
      assert(numVtx);
      for(int i=0; i<numComp; i++) newdofs[i]/=(double)numVtx;
    }
    //for(int i=0; i<6; i++)
    // cout<<newdofs[i]<<" ";
    //cout<<endl;
    // for second order derivative use nodal average
    assert(wt>=0&&wt<=1);
    for(int i=0; i<numComp/dofNode; i++)
    {
      int idx = i*dofNode+3;
      for(int j=idx; j<3; j++)
         newdofs.at(j)=wt*dofsVertex[0].at(j)+(1-wt)*dofsVertex[1].at(j);
    }
    apf::setComponents(field,vert,0,&(newdofs[0]));
  }
}


/* void  ReducedQuinticTransfer::onVertex(apf::MeshElement* parent, ma::Vector const& xi, ma::Entity* vert) */
/* { */
/*   apf::MeshEntity* oldEdge= apf::getMeshEntity(parent); */

/*   apf::MeshEntity* vertices[3]; */
/*   int num_face=-1; */
/*   apf:: Up faces; */
/*   mesh->getUp(oldEdge,faces); */
/*   num_face=faces.n; */
/*   PCU_ALWAYS_ASSERT_VERBOSE((num_face>0 && num_face<=2), "expecting 0 < num_face <= 2\n"); */


/*   // from the six dofs (phi, phi_x, phi_y, phi_xx, phi_xy, phi_yy) */
/*   // the first 3 are continuous at element boundaries. Therefore, we can */
/*   // simply use one one of the upward adjacent faces to oldEdge to interpolate */
/*   apf::MeshEntity* face = faces.e[0]; */
/*   apf::MeshEntity* dvs[3]; */
/*   mesh->getDownward(face, 0, dvs); */
/*   double coords[3][2]; */
/*   for (int i = 0; i < 3; i++) { */
/*     apf::Vector3 p; */
/*     mesh->getPoint(dvs[i], 0, p); */
/*     coords[i][0] = p[0]; */
/*     coords[i][1] = p[1]; */
/*   } */

/*   // coordinate of the new vertex */
/*   double newCoord[2]; */
/*   apf::Vector3 p; */
/*   mesh->getPoint(vert, 0, p); */
/*   newCoord[0] = p[0]; */
/*   newCoord[1] = p[1]; */


/*   for(int fi=0; fi<fields.size(); fi++) */
/*   { */
/*     apf::Field* field = fields.at(fi); */
/*     int numComp = apf::countComponents(field); */
/*     if (numComp == 1) // if size_field is added to the list take care of it now */
/*     { */
/*       apf::Element* el = apf::createElement(field, parent); */
/*       apf::setScalar(field, vert, 0, apf::getScalar(el, xi)); */
/*       apf::destroyElement(el); */
/*     } */
/*     else */
/*     { */

/*       // get all the dofs associated with "face"'s downward verts in "dvs" */

/*       // for setting the case, set the coordinates first */
/*       thecase->setCoord(coords); */

/*       // then set the dofs */
/*       apf::NewArray<double> dofs_all_fields(3*numComp); */
/*       for (int i = 0; i < 3; i++) { */
/* 	apf::getComponents(field, dvs[i], 0, &(dofs_all_fields[i*numComp])); */
/*       } */

/*       int numField=numComp/dofNode; */
/*       vector<double> newdofs(numComp); */
/*       for( int ifield=0; ifield<numField; ifield++) */
/*       { */
/* 	vector<double> dofs_one_field(3*dofNode); */
/* 	vector<double> dofs_cacu(dofNode); */
/* 	for( int i=0; i<dofNode; i++) */
/* 	{ */
/* 	  dofs_one_field.at(i)=dofs_all_fields[ifield*dofNode+i]; */
/* 	  dofs_one_field.at(dofNode+i)=dofs_all_fields[numComp+ifield*dofNode+i]; */
/* 	  dofs_one_field.at(2*dofNode+i)=dofs_all_fields[numComp*2+ifield*dofNode+i]; */
/* 	} */
/* 	thecase->setDofs(&(dofs_one_field[0])); */
/* 	thecase->eval_g(newCoord,&(dofs_cacu[0])); */
/* 	for(int i=0; i<dofNode; i++) */
/* 	  newdofs.at(ifield*dofNode+i)=dofs_cacu[i]; */
/*       } */
/*       apf::setComponents(field,vert,0,&(newdofs[0])); */
/*     } */
/*   } */
/* } */


