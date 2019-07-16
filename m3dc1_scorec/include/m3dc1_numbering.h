#ifndef M3DC1_NUMBERING_H_
#define M3DC1_NUMBERING_H_
#include "m3dc1_scorec.h"
#include <apfNumbering.h>
extern MPI_Comm M3DC1_COMM_WORLD;
void aggregateNumbering(MPI_Comm agg_cm,
                        apf::Numbering * num,
                        int nv,
                        int ndfs,
                        MPI_Comm gbl_cm = M3DC1_COMM_WORLD);
#endif
