#ifndef M3DC1_MATRIX_ALLOCATE_H_
#define M3DC1_MATRIX_ALLOCATE_H_
//#include <petscmat.h>
#include "m3dc1_mesh.h"
#include "m3dc1_matrix.h"
#include "m3dc1_field.h"
void allocateMatrix(Mat A,
                    m3dc1_matrix * mat,
                    m3dc1_mesh * msh,
                    m3dc1_field * fld);
#endif
