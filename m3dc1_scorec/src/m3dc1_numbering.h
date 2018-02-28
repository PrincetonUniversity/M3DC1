#ifndef M3DC1_NUMBERING_H_
#define M3DC1_NUMBERING_H_
#include "apfNumbering.h"
#include <mpi.h>
void aggregateNumbering(MPI_Comm cm, apf::Numbering * num, int nv, int nd);
#endif
