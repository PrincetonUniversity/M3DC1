#include "m3dc1_numbering.h"
#include "m3dc1_field.h"
#include "m3dc1_mesh.h"
#include "m3dc1_scorec.h"
#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "PCU.h"
#include <cassert>
extern MPI_Comm M3DC1_COMM_WORLD;
void aggregateNumbering(MPI_Comm agg_cm,
                        apf::Numbering * num,
                        int nv,
                        int ndfs,
                        MPI_Comm gbl_cm)
{
  MPI_Comm o_gbl_cm = PCU_Get_Comm();
  PCU_Switch_Comm(gbl_cm);
  int gbl_sz = 0;
  MPI_Comm_size(gbl_cm,&gbl_sz);
  bool lcl = (gbl_sz == 1);
  apf::Mesh2 * msh = static_cast<apf::Mesh2*>(apf::getMesh(num));
  apf::FieldShape * shp = apf::getShape(num);
  // currently this assertion will only work for non-complex fields, need to create a complex field type
  //  the function and numbering should work fine regardless
  //assert(cmps == nv * nd);
  int dim = msh->getDimension();
  int strd = 0;
  apf::MeshIterator * it = NULL;
  apf::MeshEntity * ent = NULL;
  // calculate the local stride
  for (int dd = 0; dd < dim; dd++)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if(is_ent_original(msh,ent) || lcl)
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
  MPI_Allreduce(&lcl_strd,&strd,1,MPI_INTEGER,MPI_SUM,agg_cm);
  //int dof_head = 0;
  // the first dof id is the exscan of the #of nodes times the number of dofs in the first block on the node on the agg_cm
  //int dof_id_offset = lcl_strd * ndfs;
  int nd_idx = 0;
  MPI_Exscan(&lcl_strd,&nd_idx,1,MPI_INTEGER,MPI_SUM,agg_cm);
  for(int dd = 0; dd < dim; ++dd)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if (is_ent_original(msh,ent) || lcl)
        {
          int tp = msh->getType(ent);
          int nds = shp->countNodesOn(tp);
          for(int nd = 0; nd < nds; ++nd)
          {
            for(int blk = 0; blk < nv; ++blk)
            {
              for(int dof = 0; dof < ndfs; ++dof)
              {
                int cmp = blk * ndfs + dof;
                int nbr = (blk * strd + nd_idx) * ndfs + dof;
                apf::number(num,ent,nd,cmp,nbr);
              }
            }
            nd_idx++;
          }
        }
      }
      msh->end(it);
    }
  }
  // get intra-comm offsets
  /* // won't work when cm != MPI_COMM_SELF and isn't needed when cm == MPI_COMM_SELF
  int intra_comm_offset = nd_idx; * apf::countComponents(fld); // would prefer to directly record this value
  int lcl_intra_offset = 0;
  MPI_Exscan(&intra_comm_offset,&lcl_intra_offset,1,MPI_INTEGER,MPI_SUM,cm);
  */
  // get inter-comm offsets
  // only rank zero in each comm cm will have a nonzero value
  int scp_rnk = -1;
  MPI_Comm_rank(agg_cm,&scp_rnk);
  int inter_comm_offset = (int)(!(bool)scp_rnk) * (strd * nv * ndfs);
  int lcl_offset = 0;
  // todo : replace with MSI_COMM_WORLD when we pull this into msi
  MPI_Scan(&inter_comm_offset,&lcl_offset,1,MPI_INTEGER,MPI_SUM,gbl_cm);
  // subtract local part from the scan (can't just use exscan cause the comms are arbitrary sizes
  //  and we don't want to bother making new comms)
  lcl_offset -= (strd * nv * ndfs);
  apf::setNumberingOffset(num,lcl_offset,m3dc1_mesh::instance()->get_ownership());
  if(!lcl)
    apf::synchronize(num,m3dc1_mesh::instance()->get_ownership(),false);
  PCU_Switch_Comm(o_gbl_cm);
}

