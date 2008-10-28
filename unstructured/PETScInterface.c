//#include "/u/xluo/develop.newCompiler.constraint/mctk/Examples/PPPL/PPPL/Matrix.h"
#include "/u/xluo/develop.newCompiler.constraint/mctk/Examples/PPPL/PPPL/MatrixInterface.h"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "stdio.h"
#include <stdlib.h>

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
  ierr = KSPSetOperators(*ksp, *A, *A, SAME_PRECONDITIONER /*DIFFERENT_NONZERO_PATTERN*/);CHKERRQ(ierr);
  ierr = KSPSetTolerances(*ksp, .000001, .000000001, 
			  PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
  return 0;
}


/* 
   solve2 new strategy for preconditioner 
   src/ksp/ksp/examples/tutorials/ex6f.F
   preconditioner is update every MAX_SAME_PC_COUNT times
   only for matrix with Id=5,1,6
*/
#define MAX_SAME_PC_COUNT 20
#define MAX_LINEAR_SYSTEM 3
#define SOLVE2_DEBUG 0

typedef struct
{
  int matrixId;
  int same_pc_count; // count when a pc needs to be refreshed
  Mat A; // solver matrix and preconditioner
  KSP ksp;  // solver object
} KSP_ARRAY;

// global container
KSP_ARRAY ksp_array[MAX_LINEAR_SYSTEM];

int solve2_(int *matrixId, double * rhs_sol, int * ier)
{ // local variables
  static int start=0;
  int i, whichMatrix, flag, ldb, numglobaldofs;
  int rowSize, rowId;
  int colSize, *colId;
  int valType = 0;
  PetscScalar *values;
  PetscErrorCode ierr;
  PetscTruth     flg;
  PetscLogDouble  v1,v2; 

  // petsc structure
  PC     pc; 
  Vec    u,b;


  /* at the very begining allocate space */
  if(!start) {
     for(i=0; i<MAX_LINEAR_SYSTEM; i++) {
        ksp_array[i].matrixId= -1;
        ksp_array[i].same_pc_count=0 ;
     }
     start=1;
  } // end of if(!start)

  /* find the slot for each individual solve */
  for(i=0; i<MAX_LINEAR_SYSTEM; i++) {
     if(ksp_array[i].matrixId == *matrixId ) { // found it
        whichMatrix=i;
        if(SOLVE2_DEBUG)
        PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: found matrixId_%d in ksp_array[%d] \n",
                                       *matrixId, whichMatrix);
        break;
     }else{
        if(ksp_array[i].matrixId != -1 && (i+1)<MAX_LINEAR_SYSTEM) 
           continue; // there is enough space, but someone already take this one; go to avaialable slot
        else if(ksp_array[i].matrixId == -1 )  { // this is an empty slot; then take it
           whichMatrix=i;
           ksp_array[i].matrixId = *matrixId;
           PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: add matrixId_%d into ksp_array[%d]\n",
                                          *matrixId, whichMatrix);
           
           break;
        }else{
           PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: Need to increase the size of MAX_LINEAR_SYSTEM.\n");
           exit(1);
        }
     }
  } //i<MAX_LINEAR_SYSTEM

  /* Step 0: check the matrix status */
  checkMatrixStatus_(matrixId, &flag);
  if(flag==0) return 0; // the matrix does not need to besolved

  /* Step 1: number of local rows; number of global rows */
  getMatrixLocalDofNum_(matrixId, &ldb); 
  getMatrixGlobalDofs_(matrixId, &numglobaldofs);

  /* at ntime=0,1 of each individual solve : allocate solve space */
  if(!(ksp_array[whichMatrix].same_pc_count)) { 
     /* Step 2: construct matrix */
     int *d_nnz, *o_nnz;
     d_nnz = (int*)calloc(ldb, sizeof(int));
     o_nnz = (int*)calloc(ldb, sizeof(int)); 
     getMatrixPetscDnnzOnnz_(matrixId, &valType, d_nnz, o_nnz); 
    
     ierr = MatCreateMPIAIJ(MPI_COMM_WORLD, ldb, ldb, PETSC_DECIDE, PETSC_DECIDE,
                            PETSC_NULL, d_nnz, PETSC_NULL, o_nnz, &(ksp_array[whichMatrix].A));
     CHKERRQ(ierr);
     ierr=MatSetFromOptions((ksp_array[whichMatrix].A)/*, MATSUPERLU_DIST MATMPIAIJ*/);CHKERRQ(ierr); 
     //ierr=MatSetType((ksp_array[whichMatrix].A), MATSUPERLU_DIST /*MATMPIAIJ*/);CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: create A%d into ksp_array\n", *matrixId);
  
     free(d_nnz);
     free(o_nnz);

     /* Step 3: construct ksp */
     ierr=KSPCreate(PETSC_COMM_WORLD,&(ksp_array[whichMatrix].ksp));CHKERRQ(ierr);
     ierr=KSPAppendOptionsPrefix((ksp_array[whichMatrix].ksp), "solve2_");CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: create ksp_%d into ksp_array\n", *matrixId);
  } // end of if(!(ksp_array[whichMatrix].same_pc_count))

  /* Step 4 */
  ierr = KSPGetPC((ksp_array[whichMatrix].ksp),&pc); CHKERRQ(ierr);
  ierr = KSPSetFromOptions((ksp_array[whichMatrix].ksp)); CHKERRQ(ierr);
  
  /* Step 5 */
  ierr = PetscGetTime(&v1); CHKERRQ(ierr); 
  assemblePetscMatrixNNZs_(matrixId, &valType, &(ksp_array[whichMatrix].A)); 
  ierr = MatAssemblyBegin((ksp_array[whichMatrix].A),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd((ksp_array[whichMatrix].A),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscGetTime(&v2); CHKERRQ(ierr);

  /* Step 6: construct the rhs vec */
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,ldb,PETSC_DECIDE,rhs_sol, &b);CHKERRQ(ierr);
  ierr = VecDuplicate(b, &u); CHKERRQ(ierr); 

  /* Step 7: at the start of each new cycle, re-set A as the preconditioner 
             for next MAX_SAME_PRECONDITIONER_COUNT times 
  */
     if(! ((ksp_array[whichMatrix].same_pc_count) % (MAX_SAME_PC_COUNT)) ) { 
     ierr = KSPSetOperators((ksp_array[whichMatrix].ksp),
                         (ksp_array[whichMatrix].A),
                         (ksp_array[whichMatrix].A),SAME_NONZERO_PATTERN); CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: pc_SAME_NONZERO_PATTERN for %d at %d\n",
                               *matrixId, ksp_array[whichMatrix].same_pc_count); 
     }else{
     ierr = KSPSetOperators((ksp_array[whichMatrix].ksp),
                         (ksp_array[whichMatrix].A),
                         (ksp_array[whichMatrix].A),SAME_PRECONDITIONER); CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: pc_SAME_PRECONDITIONER for %d at %d\n",
                               *matrixId, ksp_array[whichMatrix].same_pc_count); 
     }
  ierr = KSPSetInitialGuessNonzero((ksp_array[whichMatrix].ksp),PETSC_FALSE /*PETSC_TRUE*/);CHKERRQ(ierr); 
  ierr = KSPSetTolerances((ksp_array[whichMatrix].ksp), 1.e-6, 1.e-9, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(ierr);

  /* Step 8 */
  ierr = KSPSolve((ksp_array[whichMatrix].ksp), b, u); CHKERRQ(ierr); 

  // Step 9: put solution back to rhs_sol
  PetscScalar *value;
  ierr = VecGetArray(u, &value); CHKERRQ(ierr);
  for(i=0;i<ldb;i++) rhs_sol[i] = value[i];
  ierr = VecRestoreArray(u, &value); CHKERRQ(ierr); 
  
  // Step 10: set back the matrix solution to mesh
  setMatrixSoln_(matrixId, &valType, rhs_sol);
  
  /* Step 11: clean the stored data */
  cleanMatrixValues_(matrixId); 
  ierr = VecDestroy(u); CHKERRQ(ierr);
  ierr = VecDestroy(b); CHKERRQ(ierr); 
  (ksp_array[whichMatrix].same_pc_count)++;

  return 0;
}


//cj added oct 6, 2008
int dump_matrix_(int *matrixId)
{
//int fileid=33;
//writeToFile(*matrixId, fileid);

  int i, ldb;
  int valType = 0;
  int *d_nnz, *o_nnz;
  int offset=6;

     getMatrixLocalDofNum_(matrixId, &ldb); 
     d_nnz = (int*)calloc(ldb, sizeof(int));
     o_nnz = (int*)calloc(ldb, sizeof(int));
     getMatrixPetscDnnzOnnz_(matrixId, &valType, d_nnz, o_nnz);
  
     for(i=0;i<ldb;i=i+offset) 
     PetscPrintf(PETSC_COMM_WORLD, "\tvertex_%d:  %d  %d  %d  %d  %d  %d\n",
     i/offset,
     d_nnz[i+0]+o_nnz[i+0],
     d_nnz[i+1]+o_nnz[i+1],
     d_nnz[i+2]+o_nnz[i+2],
     d_nnz[i+3]+o_nnz[i+3],
     d_nnz[i+4]+o_nnz[i+4],
     d_nnz[i+5]+o_nnz[i+5]); 

     free(d_nnz);
     free(o_nnz);

  return 0;
} 

