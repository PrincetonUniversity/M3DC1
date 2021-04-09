/* 
2018-dec-17
This function is called to prepare options for 3D linear solver.
The original options is give in a filename by "-options_file" as a command line option.
This filename will be rewritten here with added "-hard_" solver options 
for "hard_" "prefix"ed solvers 
to avoid mistakes.

The added options are: lgmres and lgmres_argument.
*/

#include <petscsys.h>

void parse_solver_options_(const int *nplanes, const char *filename)
  {
  PetscErrorCode ierr;
  PetscBool      flg;
  FILE *fptr;
  const char s[2] = " ";
  char *num_of_pc_bjacobi_blocks;
  PetscInt nblocks;
  PetscInt isuperlu=0, imumps=0;

  const size_t len=256;
  char filename_out[len];

// remove heading and trailing gabrage
      char buf[len];
      sscanf(filename, "%s", buf); // Trimming on both sides occurs here
  //PetscPrintf(PETSC_COMM_WORLD,"\nsolver options from the file:%s\n\n", buf);

  char c[1024];
  //if ((fptr = fopen(filename, "r")) == NULL)
  if ((fptr = fopen(buf, "r")) == NULL)
  {
     PetscPrintf(PETSC_COMM_WORLD,"Error! opening file\n");
     exit(1);         
  }

  // reads text until newline 
  while (fgets(c, sizeof(c), fptr) != NULL)
  {
    c[strlen(c) - 1] = '\0'; // eat the newline fgets() stores
    //PetscPrintf(PETSC_COMM_WORLD,"%s\n", c);

    char *token;
    char *pc_bjacobi_blocks="-pc_bjacobi_blocks";
    char *sub_solver_type="-sub_pc_factor_mat_solver_type";
    char *sub_solver_package="-sub_pc_factor_mat_solver_package";
    char *which_solver;
    char *superlu="superlu_dist";
    char *mumps="mumps";
    /* get the first token */
    token = strtok(c, s);
    /* walk through other tokens */
    while( token != NULL ) {
       if(strcmp(token,pc_bjacobi_blocks)==0) {//matched
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", token );

          /* get the second token */
          num_of_pc_bjacobi_blocks = strtok(NULL, s);
          nblocks=atoi(num_of_pc_bjacobi_blocks);
          if(*nplanes<nblocks) {
             PetscPrintf(PETSC_COMM_WORLD, "\nError! %s: The number of Jacobi blocks %d is larger than nplanes % d in file C1input. Please change the block number in your solver options file and resubmit the job.\n", buf, nblocks, *nplanes);
               exit(1);
          }
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", num_of_pc_bjacobi_blocks );
       }
       if(strcmp(token,sub_solver_type)==0) {//matched
          which_solver = strtok(NULL, s);
          if(strcmp(which_solver,superlu)==0) { isuperlu=1; imumps=0; }
          if(strcmp(which_solver,mumps)==0) { isuperlu=0; imumps=1; }
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", which_solver );
       }
       if(strcmp(token,sub_solver_package)==0) {//matched
          which_solver = strtok(NULL, s);
          if(strcmp(which_solver,superlu)==0) { isuperlu=1; imumps=0; }
          if(strcmp(which_solver,mumps)==0) { isuperlu=0; imumps=1; }
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", which_solver );
       }

       token = strtok(NULL, s);
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", token );
    }
  }
          //PetscPrintf(PETSC_COMM_WORLD, "       %d %d %d\n", nblocks,isuperlu,imumps );
          //PetscPrintf(PETSC_COMM_WORLD, "       %s \n", num_of_pc_bjacobi_blocks );
    fclose(fptr);

  //PetscPrintf(PETSC_COMM_WORLD, "\n\n-------------------------------------\n\n");

     /* open the file for writing*/
    //sprintf(filename_out, "%s.out", filename);
    sprintf(filename_out, "%s.out", buf);
    if ((fptr = fopen(filename_out, "w")) == NULL)
    {
       PetscPrintf(PETSC_COMM_WORLD,"Error! opening file again\n");
       exit(1);         
    }
 
    //ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"c-rank%d.gp",rank);CHKERRQ(ierr);
 
    char           tmp[PETSC_MAX_PATH_LEN];
    sprintf(tmp, "%s %s", "-pc_type", "bjacobi");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %d", "-pc_bjacobi_blocks", nblocks);
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-sub_pc_type", "lu");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    if(isuperlu==1) {
#if PETSC_VERSION >= 39
       sprintf(tmp, "%s %s", "-sub_pc_factor_mat_solver_type", "superlu_dist");
#else
       sprintf(tmp, "%s %s", "-sub_pc_factor_mat_solver_package", "superlu_dist");
#endif
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    } else { //imump==1
#if PETSC_VERSION >= 39
       sprintf(tmp, "%s %s", "-sub_pc_factor_mat_solver_type", "mumps");
#else
       sprintf(tmp, "%s %s", "-sub_pc_factor_mat_solver_package", "mumps");
#endif
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
       sprintf(tmp, "%s %s", "-mat_mumps_icntl_14", "50");
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    }
    sprintf(tmp, "%s %s", "-sub_ksp_type", "preonly");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-ksp_type", "fgmres");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-ksp_gmres_restart", "220");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-ksp_max_it", "10000");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-ksp_rtol", "1.e-9");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-ksp_atol", "1.e-20");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s", "-ksp_converged_reason");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);

    sprintf(tmp, "%s %s", "-hard_pc_type", "bjacobi");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %d", "-hard_pc_bjacobi_blocks", nblocks);
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_sub_pc_type", "lu");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
#if PETSC_VERSION >= 39
    sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_type", "superlu_dist");
#else
    sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_package", "superlu_dist");
#endif
    if(isuperlu==1) {
#if PETSC_VERSION >= 39
       sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_type", "superlu_dist");
#else
       sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_package", "superlu_dist");
#endif
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    } else { //imump==1
#if PETSC_VERSION >= 39
       sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_type", "mumps");
#else
       sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_solver_package", "mumps");
#endif
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
       sprintf(tmp, "%s %s", "-hard_sub_pc_factor_mat_mumps_icntl_14", "50");
          //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    }
    sprintf(tmp, "%s %s", "-hard_sub_ksp_type", "preonly");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_type", "fgmres");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_lgmres_augment", "4");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_gmres_restart", "220");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_max_it", "10000");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_rtol", "1.e-9");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s %s", "-hard_ksp_atol", "1.e-20");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);
    sprintf(tmp, "%s", "-hard_ksp_converged_reason");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);

    sprintf(tmp, "%s", "-on_error_abort");
       //PetscPrintf(PETSC_COMM_WORLD,"%s\n", tmp);
          PetscFPrintf(PETSC_COMM_WORLD,fptr,"%s\n", tmp);

   /* close the file*/  
   fclose (fptr);
}
