#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "mpi.h"
#include "stdio.h"

/* the matrix id's need to correspond to the same id's in
   the sparse module of M3Dmodules.f90 */
enum FortranMatrixID {
  s6matrix_sm = 1,
  s8matrix_sm = 2,
  s7matrix_sm = 3,
  s4matrix_sm = 4,
  s3matrix_sm = 5,
  s5matrix_sm = 6,
  s1matrix_sm = 7,
  s2matrix_sm = 8,
  d1matrix_sm = 9,
  d2matrix_sm = 10,
  d4matrix_sm = 11,
  d8matrix_sm = 12,
  q1matrix_sm = 13,
  r2matrix_sm = 14,
  r8matrix_sm = 15,
  q2matrix_sm = 16,
  q8matrix_sm = 17,
  gsmatrix_sm = 18,
  s9matrix_sm = 19,
  d9matrix_sm = 20,
  r9matrix_sm = 21,
  q9matrix_sm = 22,
  r14matrix_sm = 23,
  s10matrix_sm = 24,
  d10matrix_sm = 25,
  q10matrix_sm = 26,
  r10matrix_sm = 27
};

/* 
   below sets the PETSc matrix type for the matrix with id "matrixid"
   PETSc documentation is available at 
   http://www-unix.mcs.anl.gov/petsc/petsc-as/documentation/index.html
*/
int setPETScMat(int matrixid, Mat * A) {
  PetscErrorCode ierr;
  PetscInt ipetscmpiaij =1;
  PetscInt ipetscsuperlu=0;

  PetscOptionsGetInt(PETSC_NULL,"-ipetscsuperlu",&ipetscsuperlu,PETSC_NULL);
  if(ipetscsuperlu) ipetscmpiaij=0;

  if(ipetscmpiaij) {  /* create default distributed matrix type  */
    ierr = MatSetType(*A, MATMPIAIJ);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScMat %d to MATMPIAIJ\n", matrixid); 
  }
  else {  /* create specialized matrix for SuperLU/SuperLU_DIST */
    ierr = MatSetType(*A, MATSUPERLU_DIST);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScMat %d to MATSUPERLU_DIST\n", matrixid); 
  }
  return 0;
}

/* 
   below sets the PETSc options for the solver for matrix with id "matrixid"
   PETSc documentation is available at 
   http://www-unix.mcs.anl.gov/petsc/petsc-as/documentation/index.html
*/
int setPETScKSP(int matrixid, KSP * ksp, Mat * A) {
  PetscErrorCode ierr;
  ierr = KSPCreate(MPI_COMM_WORLD, ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(*ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetTolerances(*ksp, .000001, .000000001, 
			  PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
  return 0;
}
