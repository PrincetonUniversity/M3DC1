//#include "/p/tsc/m3dc1/lib/develop.newCompiler.constraint/mctk/Examples/PPPL/PPPL/Matrix.h"
//#include "/p/tsc/m3dc1/lib/develop.newCompiler.constraint/mctk/Examples/PPPL/PPPL/MatrixInterface.h"
//#include "/p/tsc/m3dc1/lib/develop.newCompiler.petscDev/mctk/Examples/PPPL/PPPL/MatrixInterface.h"
#include <MatrixInterface.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscvec.h>

#include <stdio.h>
#include <stdlib.h>

#ifdef USEHYBRID
#ifdef USECOMPLEX
#include "zpdslin_solver.h"
#include "zpdslin_util.h"
#include "pdslin_solver.h"
#else
#include "dpdslin_solver.h"
#include "dpdslin_util.h"
#include "pdslin_solver.h"
#endif 

#else
#include "superlu_defs.h"
#endif

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
   below sets the Superlu_Dist Solver options May 02, 2011
*/
int setSuperluOptions(int matrixid, superlu_options_t * options) {
  //  (*options).ColPerm = NATURAL ;
//    PetscPrintf(PETSC_COMM_WORLD, "\tsetSuperluOptions to ConditionNumber = YES\n");
  return 0;
}

/* 
   below sets the PETSc matrix type for the matrix with id "matrixid"
   PETSc documentation is available at 
   http://www-unix.mcs.anl.gov/petsc/petsc-as/documentation/index.html
*/
int setPETScMat(int matrixid, Mat * A) {
  PetscErrorCode ierr;
//PetscInt ipetscmpiaij =1;
//PetscInt ipetscsuperlu=0;

//PetscOptionsGetInt(PETSC_NULL,"-ipetscsuperlu",&ipetscsuperlu,PETSC_NULL);
//if(ipetscsuperlu) ipetscmpiaij=0;

//if(ipetscmpiaij) {  /* create default distributed matrix type  */
    ierr = MatSetType(*A, MATMPIAIJ);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScMat %d to MATMPIAIJ\n", matrixid); 
//}
//else {  /* create specialized matrix for SuperLU/SuperLU_DIST */
//  ierr = MatSetType(*A, MATSUPERLU_DIST);CHKERRQ(ierr);
//  PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScMat %d to MATSUPERLU_DIST\n", matrixid); 
//}
    ierr = MatSetFromOptions(*A);CHKERRQ(ierr); 
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
  PetscPrintf(PETSC_COMM_WORLD, "\tsetPETScKSP for %d\n", matrixid); 

  return 0;
}


/*
   dec 1, 2010 cj
   get_counter fix newvar_solve matrix print bug
*/
typedef struct
{
  int matrixId;
  int counter; // count when a pc needs to be refreshed
} MAT_COUNTER;
// global container
MAT_COUNTER matcounter[40];
int get_counter_(int *matrixId, int *counter)
{ // local variables
  static int start=0;
  int i;
  /* at the very begining */
  if(!start) {
     for(i=0; i<40; i++) {
        matcounter[i].matrixId= -1;
        matcounter[i].counter=0 ;
     }
     start=1;
  } // end of if(!start)
  /* find the slot for each individual matrix */
  for(i=0; i<40; i++) {
     if(matcounter[i].matrixId == *matrixId ) { // found it
        matcounter[i].counter ++ ;
        break;
     }else{
        if(matcounter[i].matrixId != -1 && (i+1)<40)
           continue; // there is enough space, but someone already take this one; go to avaialable slot
        else if(matcounter[i].matrixId == -1 )  { // this is an empty slot; then take it
           matcounter[i].matrixId = *matrixId;
           break;
        }else{
           PetscPrintf(PETSC_COMM_WORLD, "\tget_counter: Need to increase the container size.\n");
           exit(1);
        }
     }
  } //i<40 
  *counter = matcounter[i].counter;
  PetscPrintf(PETSC_COMM_WORLD, "\tget_counter: found matrix %d %d times.\n", *matrixId, *counter);
  return 0;
}


/* 
   solve2 new strategy for preconditioner 
   src/ksp/ksp/examples/tutorials/ex6f.F
   preconditioner is update every MAX_SAME_PC_COUNT times
   only for matrix with Id=5,1,6
*/
#define MAX_SAME_PC_COUNT 1
#define MAX_LINEAR_SYSTEM 40
#define SOLVE2_DEBUG 1

typedef struct
{
  int matrixId;
  int same_pc_count; // count when a pc needs to be refreshed
  int ldb, numglobaldofs;
  Mat A; // solver matrix and preconditioner
  KSP ksp;  // solver object
} KSP_ARRAY;

// global container
KSP_ARRAY ksp_array[MAX_LINEAR_SYSTEM];

int solve2_(int *matrixId, double * rhs_sol, int * ier)
{ // local variables
  static int start=0;
  int i, whichMatrix, flag, its;
  int rowSize, rowId;
  int colSize, *colId;
  int valType = 0;
  PetscScalar *values;
  PetscReal rms, norm;
  PetscErrorCode ierr;
#ifdef PetscDEV
  PetscBool     flg;
#else
  PetscTruth     flg;
#endif
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

  /* at ntime=0,1 of each individual solve : allocate solve space */
  if(!(ksp_array[whichMatrix].same_pc_count)) { 
     /* Step 1: number of local rows; number of global rows */
     getMatrixLocalDofNum_(matrixId, &(ksp_array[whichMatrix].ldb)); 
     getMatrixGlobalDofs_(matrixId, &(ksp_array[whichMatrix].numglobaldofs)); 
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_%d: numglobaldofs_%d ldb_%d \n",
                 *matrixId, ksp_array[whichMatrix].numglobaldofs, ksp_array[whichMatrix].ldb);

     /* Step 2: construct matrix */
     int *d_nnz, *o_nnz;
     d_nnz = (int*)calloc(ksp_array[whichMatrix].ldb, sizeof(int));
     o_nnz = (int*)calloc(ksp_array[whichMatrix].ldb, sizeof(int)); 
     getMatrixPetscDnnzOnnz_(matrixId, &valType, d_nnz, o_nnz); 
     /*for(i=0; i<ksp_array[whichMatrix].ldb; i++)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_%d_%d: d_nnz_%d o_nnz_%d \n", 
                                   *matrixId, i+1, d_nnz[i], o_nnz[i]);*/
    
     ierr = MatCreateMPIAIJ(MPI_COMM_WORLD, 
                            ksp_array[whichMatrix].ldb, ksp_array[whichMatrix].ldb,
                            PETSC_DECIDE, PETSC_DECIDE,
                            PETSC_NULL, d_nnz, PETSC_NULL, o_nnz,
                            &(ksp_array[whichMatrix].A));
     CHKERRQ(ierr);
     ierr=MatSetFromOptions((ksp_array[whichMatrix].A)/*, MATSUPERLU_DIST MATMPIAIJ*/);CHKERRQ(ierr); 
     //ierr=MatSetType((ksp_array[whichMatrix].A), MATSUPERLU_DIST /*MATMPIAIJ*/);CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: create A%d into ksp_array\n", *matrixId);
  
     free(d_nnz);
     free(o_nnz);

     /* Step 3: construct ksp */
     ierr=KSPCreate(PETSC_COMM_WORLD,&(ksp_array[whichMatrix].ksp));CHKERRQ(ierr);
     /*
     if(*matrixId==5 || *matrixId==6)
     ierr=KSPAppendOptionsPrefix((ksp_array[whichMatrix].ksp), "solve2_");CHKERRQ(ierr); 
     if(*matrixId==4 || *matrixId==2)
     ierr=KSPAppendOptionsPrefix((ksp_array[whichMatrix].ksp), "solve1_");CHKERRQ(ierr); 
     */
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: create ksp_%d into ksp_array\n", *matrixId);

     /* Step 4 */
     ierr = KSPGetPC((ksp_array[whichMatrix].ksp),&pc); CHKERRQ(ierr);
     ierr = KSPSetFromOptions((ksp_array[whichMatrix].ksp)); CHKERRQ(ierr);
  } // end of if(!(ksp_array[whichMatrix].same_pc_count))
  
  if(flag==1) {
     /* Step 5 */
     ierr = PetscGetTime(&v1); CHKERRQ(ierr); 
     assemblePetscMatrixNNZs_(matrixId, &valType, &(ksp_array[whichMatrix].A)); 
     ierr = MatAssemblyBegin((ksp_array[whichMatrix].A),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd((ksp_array[whichMatrix].A),MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = PetscGetTime(&v2); CHKERRQ(ierr);
  } // only if flag=1 we update matrix

  /* Step 6: construct the rhs vec
     ierr = PetscMalloc(ldb*sizeof(PetscScalar),&tmp);CHKERRQ(ierr); 
     ierr = PetscMemcpy((void *)tmp,rhs_sol,ldb*sizeof(PetscScalar));CHKERRQ(ierr);
     ierr = PetscFree(tmp);CHKERRQ(ierr);
     const PetscScalar *tmp;
     tmp=&(rhs_sol[0]); 
     ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,ldb,PETSC_DECIDE,tmp, &b);CHKERRQ(ierr);
  */
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,ksp_array[whichMatrix].ldb,
                               PETSC_DECIDE,rhs_sol, &b);CHKERRQ(ierr);
  ierr = VecDuplicate(b, &u); CHKERRQ(ierr); 

  /* Step 7: at the start of each new cycle, re-set A as the preconditioner 
             for next MAX_SAME_PRECONDITIONER_COUNT times 
  */
     if(! ((ksp_array[whichMatrix].same_pc_count) % (MAX_SAME_PC_COUNT)) &&
        flag==1) { 
     ierr = KSPSetOperators((ksp_array[whichMatrix].ksp),
                         (ksp_array[whichMatrix].A),
                         (ksp_array[whichMatrix].A),SAME_NONZERO_PATTERN); CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: new preconditioner for %d at %d\n",
                               *matrixId, ksp_array[whichMatrix].same_pc_count); 
     }else{
     ierr = KSPSetOperators((ksp_array[whichMatrix].ksp),
                         (ksp_array[whichMatrix].A),
                         (ksp_array[whichMatrix].A),SAME_PRECONDITIONER); CHKERRQ(ierr);
     if(SOLVE2_DEBUG)
     PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: old preconditioner for %d at %d\n",
                               *matrixId, ksp_array[whichMatrix].same_pc_count); 
     }
  ierr = KSPSetInitialGuessNonzero((ksp_array[whichMatrix].ksp),PETSC_FALSE /*PETSC_TRUE*/);
  CHKERRQ(ierr); 
  //ierr = KSPSetTolerances((ksp_array[whichMatrix].ksp), 1.e-6, 1.e-9, PETSC_DEFAULT, PETSC_DEFAULT);
  ierr = KSPSetTolerances((ksp_array[whichMatrix].ksp), 
                          PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);

  /* Step 8 */
  ierr = KSPSolve((ksp_array[whichMatrix].ksp), b, u); CHKERRQ(ierr); 
  ierr = KSPGetIterationNumber(ksp_array[whichMatrix].ksp, &its); CHKERRQ(ierr);
  ierr = MatNorm(ksp_array[whichMatrix].A,NORM_1,&norm); CHKERRQ(ierr); 
  ierr = VecNorm(u,NORM_2,&rms); CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD, "\tsolve2_: %d its_%d rms = %e matnorm = %e\n",
              *matrixId, its, rms/(float)ksp_array[whichMatrix].numglobaldofs, norm);

  // Step 9: put solution back to rhs_sol
  PetscScalar *value;
  ierr = VecGetArray(u, &value); CHKERRQ(ierr);
  for(i=0;i<ksp_array[whichMatrix].ldb;i++) rhs_sol[i] = value[i];
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


#ifdef USEHYBRID
typedef struct
{
  int matrixId;
  int init; // count when a pc needs to be refreshed
  dPDSLinMatrix A;
  pdslin_param  input;
  pdslin_stat  stat;
  int  info;
  double *x_loc, *b_loc; 
  int lnnz;
} HYBRIDSOLVER_ARRAY; 
// global container
HYBRIDSOLVER_ARRAY hs_array[MAX_LINEAR_SYSTEM];
int hybridsolve_(int *matrixId, double *rhs_sol, int *ier)
{ // local variables
  static int start=0, myrank, size;
  int rowSize, rowId, colSize, *colId;
  int valType = 0;
  double *values;
  int i, j, flag, counter, whichMatrix;

  int n, m, mloc, frow, nnz, lnnz;
  int *lrowptr, *lcolind;
  double *lnzval;
  double *oneb_loc, *onex_loc, rms, grms, oneb, onex, goneb, gonex;


  /* at the very begining allocate space */
  if(!start) {
     for(i=0; i<MAX_LINEAR_SYSTEM; i++) {
        hs_array[i].matrixId= -1;
        hs_array[i].init=0 ;
     }
     start=1;
  } // end of if(!start)

  /* find the slot for each individual hybrid solve */
  for(i=0; i<MAX_LINEAR_SYSTEM; i++) {
     if(hs_array[i].matrixId == *matrixId ) { // found it
        whichMatrix=i;
        if(SOLVE2_DEBUG)
        PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: found matrixId_%d in hs_array[%d] \n",
                                       *matrixId, whichMatrix);
        break;
     }else{
        if(hs_array[i].matrixId != -1 && (i+1)<MAX_LINEAR_SYSTEM) 
           continue; // there is enough space, but someone already take this one; go to avaialable slot
        else if(hs_array[i].matrixId == -1 )  { // this is an empty slot; then take it
           whichMatrix=i;
           hs_array[i].matrixId = *matrixId;
           PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: add matrixId_%d into hs_array[%d]\n",
                                          *matrixId, whichMatrix);
           
           break;
        }else{
           PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: Need to increase the size of MAX_LINEAR_SYSTEM.\n");
           exit(1);
        }
     }
  } //i<MAX_LINEAR_SYSTEM

  /* Step 1: check the matrix status */
  checkMatrixStatus_(matrixId, &flag);
  if(flag==0) return 0; // the matrix does not need to besolved
  PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d status_%d\n", *matrixId, flag); 

  /* Step 8: hybrid solve */
      /* 8.0 hybrid solve: clean solver content 
         dpdslin_solver( double *b0, double *x, dPDSLinMatrix *matrix, pdslin_param *input, pdslin_stat *stat, int *info );
       */
        if(hs_array[whichMatrix].init >= 1 && flag==1) {
      hs_array[whichMatrix].input.job =  PDSLin_CLEAN;
      dpdslin_solver( hs_array[whichMatrix].b_loc, 
                      hs_array[whichMatrix].x_loc, 
                      &(hs_array[whichMatrix].A), 
                      &(hs_array[whichMatrix].input),
                      &(hs_array[whichMatrix].stat),
                      &(hs_array[whichMatrix].info) );
        } //if(hs_array[whichMatrix].init && flag==1) {

        if(hs_array[whichMatrix].init==0) {
      /* 8.1 hybrid solve: initialize MPI and hybrid solver 
         dpdslin_init( &input, &matrix, &stat, pdslin_comm, argc, argv );
       */
      dpdslin_init( &(hs_array[whichMatrix].input),
                    &(hs_array[whichMatrix].A),
                    &(hs_array[whichMatrix].stat),
                    PETSC_COMM_WORLD,
                    0,
                    NULL );
      hs_array[whichMatrix].input.verbose = PDSLin_YES;
      //myrank=hs_array[whichMatrix].input.proc_id;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
        } //if(hs_array[whichMatrix].init==0) {

      /* 8.2 hybrid solve: set up coefficient matrix in distributed CSR format  *
       * matrix.n  global matrix dimension                        *
       * matrix.m  global matrix dimension                        *
       * matrix.nnz total number of nnz in the global matrix      *
       * matrix.frow first row index of the local matrix          *
       * matrix.mloc local matrix size                            *
       * matrix.lrowptr pointers to the begining of each rows     *
       * matrix.lcolind column indexes of the local matrix        *
       * matrix.lnzval nonzero values of the local matrix         */

        if(flag==1) {
  /* Step 2: get the local dofs */
  getMatrixLocalDofNum_(matrixId, &mloc);
  getMatrixGlobalDofs_(matrixId, &n );
  getMatrixGlobalDofs_(matrixId, &m);
  getMatrixFirstDof_(matrixId, &frow);
  hs_array[whichMatrix].A.n=n;
  hs_array[whichMatrix].A.mloc=mloc;
  hs_array[whichMatrix].A.frow=frow - 1; // shift one from fortran index
  printf("\thybrid_: %d_matrixId_%d dim_%d %d %d\n", myrank, *matrixId, n, mloc, frow); 

  // Step 3: allocate memory for rowptr
  lrowptr = (int *)malloc((unsigned)(mloc+1)*sizeof(int));
     if(!lrowptr) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed lrowptr.\n");
  hs_array[whichMatrix].A.lrowptr = (int *)malloc((unsigned)(mloc+1)*sizeof(int));
     if(!(hs_array[whichMatrix].A.lrowptr)) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed hs_array[].A.lrowptr.\n");

  // Step 4: set up lrowptr and lnnz, global nnz
  lrowptr[0] = 1;
  counter = 1;
  lnnz = 0;
  getMatrixNNZRowSize_(matrixId, &valType, &rowSize);
  printf("\thybrid_: %d_matrixId_%d rowSize_%d\n", myrank, *matrixId, rowSize);
  for(i=1; i<=rowSize; i++) {
     getMatrixNNZRowId_(matrixId, &valType, &i, &rowId);
     getMatrixNNZColSize_(matrixId, &valType, &rowId, &colSize);
     lrowptr[counter] = colSize + lrowptr[counter-1];
     lnnz += colSize;
     //printf("\thybrid_: %d_matrixId_%d row_%d col_%d lrowptr_%d lnnz_%d\n", 
     //       myrank, *matrixId, rowId, colSize, lrowptr[counter], lnnz); 
     counter++;
  }
  hs_array[whichMatrix].lnnz = lnnz;
  /* Combines values from all processes and 
     distributes the result back to all processes */
  MPI_Allreduce(&lnnz, &nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  hs_array[whichMatrix].A.nnz=nnz;
  for(i=0; i<=rowSize; i++) hs_array[whichMatrix].A.lrowptr[i]=lrowptr[i] - 1; // shift on from fortran index
  free(lrowptr);
  PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d nnz_%d\n", *matrixId, nnz); 
        } //if(flag==1) {


        if(flag==1) {
  // Step 5: allocate memory for lnnz
  lnzval = (double *)malloc((unsigned) (hs_array[whichMatrix].lnnz)*sizeof(double));
     if(!lnzval) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed nzval.\n");
  lcolind = (int *)malloc((unsigned) ((hs_array[whichMatrix].lnnz)*sizeof(int)));
     if(!lcolind) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed colind.\n");
  hs_array[whichMatrix].A.lnzval = (double *)malloc((unsigned) (lnnz)*sizeof(double));
     if(!(hs_array[whichMatrix].A.lnzval)) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed hs_array[].A.nzval.\n");
  hs_array[whichMatrix].A.lcolind = (int *)malloc((unsigned) ((lnnz)*sizeof(int)));
     if(!(hs_array[whichMatrix].A.lcolind)) PetscPrintf(PETSC_COMM_WORLD, "Malloc failed hs_array[].A.colind.\n");


  // Step 6: put the computed values to the nzval array
  //getMatrixNNZRowSize_(matrixId, &valType, &rowSize);
  counter = 0;
  for(i=1; i<=rowSize; i++) {
     getMatrixNNZRowId_(matrixId, &valType, &i, &rowId);
     getMatrixNNZColSize_(matrixId, &valType, &rowId, &colSize);
     colId = malloc(colSize*sizeof(int));
     values = malloc(colSize*sizeof(double));
     getMatrixNNZValues_(matrixId, &valType, &rowId, colId, values);
     //memcpy(lnzval+counter, values, colSize); // copy the nonzeros to the nzval
     //memcpy(lcolind+counter, colId, colSize); // copy the colId to the lcolind
     for(j=0; j<colSize; j++) {
        lnzval[counter+j]=values[j];
        lcolind[counter+j]=colId[j];
        //if(hs_array[whichMatrix].init==0 && *matrixId == 5)
        //printf("\thybrid_: %d_matrixId_%d row_%d col_%d_%d lnzval_%e\n", 
        //      myrank, *matrixId, rowId, colSize, colId[j], lnzval[counter+j]); 
     }
     counter += colSize;
     free(colId);
     free(values);
  }
  if(counter != hs_array[whichMatrix].lnnz)
     printf( "The nnz count has difference in %d %d %d ", myrank, counter, hs_array[whichMatrix].lnnz);
  else
     PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d fill colind\n", *matrixId, counter); 
  for(i=0; i<hs_array[whichMatrix].lnnz; i++) hs_array[whichMatrix].A.lnzval[i]=lnzval[i];
  for(i=0; i<hs_array[whichMatrix].lnnz; i++) hs_array[whichMatrix].A.lcolind[i]=lcolind[i]-1; // shift one from fortran index
  free(lnzval);
  free(lcolind); 
        } //if(flag==1) {


  /* Step 8: hybrid solve */
        if(hs_array[whichMatrix].init==0) {
      /* 8.3 hybrid solve: allocate right-hand-side and solution vectors */
      hs_array[whichMatrix].x_loc = pdslin_dmalloc( hs_array[whichMatrix].A.mloc );
      hs_array[whichMatrix].b_loc = pdslin_dmalloc( hs_array[whichMatrix].A.mloc );
  PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d allocate x_loc b_loc (%d)\n", *matrixId, mloc); 
        } //if(hs_array[whichMatrix].init==0) {

        if(flag==1) {
      /* 8.4 hybrid solve: call hybrid solver to compute preconditioner */
      /*hs_array[whichMatrix].input.num_doms = 0;   1: superlu_dist
                                                    0: petsc
                                                    default: pdslin
                                                 */
      /*hs_array[whichMatrix].input.perm_dom = NATURAL;  NATURAL: use the natural ordering 
                                                         MMD_ATA: use minimum degree ordering on structure of A'*A
                                                         MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
                                                         COLAMD: use approximate minimum degree column ordering
                                                         MY_PERMC: use the ordering specified by the user
                                                        */
      hs_array[whichMatrix].input.verbose = PDSLin_VALL;
      hs_array[whichMatrix].input.remove_zero = PDSLin_NO;
      hs_array[whichMatrix].input.job = PDSLin_PRECO;
      //if(*matrixId==5) hs_array[whichMatrix].input.drop_tau0 = 1.e-8;
      dpdslin_solver( hs_array[whichMatrix].b_loc, 
                      hs_array[whichMatrix].x_loc, 
                      &(hs_array[whichMatrix].A), 
                      &(hs_array[whichMatrix].input),
                      &(hs_array[whichMatrix].stat),
                      &(hs_array[whichMatrix].info) ); 
      /* wait for all the processors to finish factorizing */
      MPI_Barrier( MPI_COMM_WORLD ); 
  PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d calculate preconditioner and subdomain info_%d\n", 
                                *matrixId, hs_array[whichMatrix].info); 
        } //if(flag==1) {

      /* 8.5 hybrid solve: setup the right-hand-side */
      rms=0.;
      grms=0.;
      for( i=0; i<hs_array[whichMatrix].A.mloc; i++ ) {
         hs_array[whichMatrix].b_loc[i] = rhs_sol[i];
         rms += rhs_sol[i] * rhs_sol[i] ;
      }
      MPI_Allreduce(&rms, &grms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: b_rms = %d %d %e %e\n", 
                  *matrixId, hs_array[whichMatrix].A.n, grms, sqrt(grms)/(float)hs_array[whichMatrix].A.n);

      /* 8.6 hybrid solve: call hybrid solver to compute solution */
      hs_array[whichMatrix].input.job = PDSLin_SOLVE;
      dpdslin_solver( hs_array[whichMatrix].b_loc, 
                      hs_array[whichMatrix].x_loc, 
                      &(hs_array[whichMatrix].A), 
                      &(hs_array[whichMatrix].input),
                      &(hs_array[whichMatrix].stat),
                      &(hs_array[whichMatrix].info) );
      rms=0.;
      grms=0.;
      for( i=0; i<hs_array[whichMatrix].A.mloc; i++ ) {
         rhs_sol[i] = hs_array[whichMatrix].x_loc[i];
         rms += rhs_sol[i] * rhs_sol[i] ;
         //PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: sol= %d %d %e\n", *matrixId, i, rhs_sol[i]);
      }
      MPI_Allreduce(&rms, &grms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: sol_rms = %d %d %e %e\n", 
                  *matrixId, hs_array[whichMatrix].A.n, grms, sqrt(grms)/(float)hs_array[whichMatrix].A.n);
  
      /* 8.7 hybrid solve: free right-hand-side and solution
      hybrid_free(hs_array[whichMatrix].x_loc);
      hybrid_free(hs_array[whichMatrix].b_loc);
      */
  
      /* 8.8 hybrid solve: print out statistics */
      pdslin_print_stat( &(hs_array[whichMatrix].stat), PETSC_COMM_WORLD );

      /* 8.9 hybrid solve: terminate MPI and hybrid solver 
      dhybrid_finalize(&(hs_array[whichMatrix].input), &(hs_array[whichMatrix].A));
      */
  
  // Step 10: set back the matrix solution to mesh
  setMatrixSoln_(matrixId, &valType, rhs_sol);

  /* Step 11: clean the stored data */
  cleanMatrixValues_(matrixId); 
/*free(hs_array[whichMatrix].A.lrowptr);
  free(hs_array[whichMatrix].A.lnzval);
  free(hs_array[whichMatrix].A.lcolind); 
*/
  hs_array[whichMatrix].init ++;
  PetscPrintf(PETSC_COMM_WORLD, "\thybrid_: matrixId_%d called %d times.\n", *matrixId, hs_array[whichMatrix].init); 

  return 0;
}
#endif 
