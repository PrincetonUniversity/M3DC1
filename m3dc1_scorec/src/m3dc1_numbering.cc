#include "m3dc1_numbering.h"
#include "m3dc1_mesh.h" // temp for ownership
#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
void aggregateNumbering(MPI_Comm cm, apf::Numbering * num, int nv, int ndfs)
{
  apf::Field * fld = apf::getField(num);
  apf::Mesh2 * msh = static_cast<apf::Mesh2*>(apf::getMesh(fld));
  apf::FieldShape * shp = apf::getShape(fld);
  //int cmps = apf::countComponents(fld);
  // currently this assertion will only work for non-complex fields, need to create a complex field type
  //  the function and numbering should work fine regardless
  //assert(cmps == nv * nd);
  int dim = msh->getDimension();
  int strd = 0;
  apf::MeshIterator * it = NULL;
  apf::MeshEntity * ent = NULL;
  // calculate the local stride
  for(int dd = 0; dd < dim; dd++)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if(is_ent_original(msh,ent))
        {
          int tp = msh->getType(ent);
          int nds = shp->countNodesOn(tp);
          strd += nds;
        }
      }
      msh->end(it);
    }
  }
  // sum the stride across the comm
  int lcl_strd = strd;
  MPI_Allreduce(&lcl_strd,&strd,1,MPI_INTEGER,MPI_SUM,cm);
  //int dof_head = 0;
  int nd_idx = 0;
  for(int dd = 0; dd < dim; ++dd)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if(is_ent_original(msh,ent))
        {
          int tp = msh->getType(ent);
          int nds = shp->countNodesOn(tp);
          for(int nd = 0; nd < nds; ++nd)
          {
            for(int var = 0; var < nv; ++var)
            {
              for(int dof = 0; dof < ndfs; ++dof)
              {
                int cmp = var * ndfs + dof;
                int nbr = var * strd + nd_idx * ndfs + dof;
                apf::number(num,ent,nd,cmp,nbr);
              }
            }
            nd_idx++;
          }
        }
      }
    }
  }
  // get intra-comm offsets
  int rnk = -1;
  MPI_Comm_rank(cm,&rnk);
  int intra_comm_offset = 0;
  int lcl_intra_offset = 0;
  MPI_Scan(&intra_comm_offset,&lcl_intra_offset,1,MPI_INTEGER,MPI_SUM,cm);
  // get inter-comm offsets
  // only rank zero in each comm cm will have a nonzero value
  int inter_comm_offset = (int)(!(bool)rnk) * (lcl_strd * nv * ndfs);
  int lcl_offset = 0;
  // todo : replace MPI_COMM_WORLD WITH MSI_COMM_WORLD when we pull this into msi
  MPI_Scan(&inter_comm_offset,&lcl_offset,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
  apf::synchronize(num);
}
