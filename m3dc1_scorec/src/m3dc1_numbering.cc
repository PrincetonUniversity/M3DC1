#include "m3dc1_numbering.h"
#include "m3dc1_mesh.h" // temp for ownership
#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "PCU.h"
#include <cassert>
extern MPI_Comm M3DC1_COMM_WORLD;

//*******************************************************
void synchronize_numbering(apf::Mesh2* msh, apf::Numbering * num, 
                           apf::FieldShape * shp, int nv, int ndfs)
//*******************************************************
{
  apf::MeshEntity* ent;
  int nnbr=nv*ndfs;
  int* nbr = new int[nnbr];
  int dim = msh->getDimension();
  apf::MeshIterator * it = NULL;

  PCU_Comm_Begin();

  for(int dd = 0; dd < dim; ++dd)
  {
    if (shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        if (!is_ent_original(msh,ent)) continue;

        int tp = msh->getType(ent);
        int nds = shp->countNodesOn(tp);

        for (int nd = 0; nd < nds; ++nd)
        {
          for (int cmp = 0; cmp < nv*ndfs; ++cmp)
            nbr[cmp]=apf::getNumber(num,ent,nd,cmp);
        }
        if (msh->isShared(ent))
        {
          apf::Copies remotes;
          msh->getRemotes(ent,remotes);
          APF_ITERATE(apf::Copies,remotes,rit)
          {
            PCU_COMM_PACK(rit->first,rit->second);
            PCU_Comm_Pack(rit->first,&(nbr[0]), nnbr*sizeof(int));
          }
        }
        if (msh->isGhosted(ent))
        {
          apf::Copies ghosts;
          msh->getGhosts(ent,ghosts);
          APF_ITERATE(apf::Copies,ghosts,git)
          {
            PCU_COMM_PACK(git->first,git->second);
            PCU_Comm_Pack(git->first,&(nbr[0]),nnbr*sizeof(int));
          }
        }
      } // while
      msh->end(it);
    } // if
  } // for 
  delete [] nbr;

  PCU_Comm_Send();

  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      int* recv_nbr = new int[nnbr];
      PCU_Comm_Unpack(&(recv_nbr[0]),nnbr*sizeof(int));

      int tp = msh->getType(r);
      int nds = shp->countNodesOn(tp);
      for(int nd = 0; nd < nds; ++nd)
        for (int cmp = 0; cmp < nv*ndfs; ++cmp)
          apf::number(num, r, nd, cmp, recv_nbr[cmp]);
      delete [] recv_nbr;
    }
}

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
  // the first dof id is the exscan of the #of nodes times the number of dofs in the first block on the node on the cm
  int dof_id_offset = lcl_strd * ndfs;
  int nd_idx = 0;
  MPI_Exscan(&dof_id_offset,&nd_idx,1,MPI_INTEGER,MPI_SUM,cm);
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
            for(int blk = 0; blk < nv; ++blk)
            {
              for(int dof = 0; dof < ndfs; ++dof)
              {
                int cmp = blk * ndfs + dof;
                int nbr = blk * strd + nd_idx * ndfs + dof;
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
  int rnk = -1;
  MPI_Comm_rank(cm,&rnk);
  /* // won't work when cm != MPI_COMM_SELF and isn't needed when cm == MPI_COMM_SELF
  int intra_comm_offset = nd_idx; * apf::countComponents(fld); // would prefer to directly record this value
  int lcl_intra_offset = 0;
  MPI_Exscan(&intra_comm_offset,&lcl_intra_offset,1,MPI_INTEGER,MPI_SUM,cm);
  */
  // get inter-comm offsets
  // only rank zero in each comm cm will have a nonzero value
  int inter_comm_offset = (int)(!(bool)rnk) * (lcl_strd * nv * ndfs);
  int lcl_offset = 0;
  // todo : replace MPI_COMM_WORLD WITH MSI_COMM_WORLD when we pull this into msi
  MPI_Exscan(&inter_comm_offset,&lcl_offset,1,MPI_INTEGER,MPI_SUM,M3DC1_COMM_WORLD);
  apf::SetNumberingOffset(num,lcl_offset);
  // this does work with PUMI default ownershiop
  apf::synchronize(num);

  // synchronize numbering manually
  // synchronize_numbering(msh, num, shp, nv, ndfs);

#ifdef DEBUG
  for(int dd = 0; dd < dim; ++dd)
  {
    if(shp->hasNodesIn(dd))
    {
      it = msh->begin(dd);
      while((ent = msh->iterate(it)))
      {
        int tp = msh->getType(ent);
        int nds = shp->countNodesOn(tp);
        for(int nd = 0; nd < nds; ++nd)
          for(int dof = 0; dof < nv*ndfs; ++dof)
            assert(apf::isNumbered(num,ent,nd,dof));
      }
      msh->end(it);
    }
  }
#endif
}
