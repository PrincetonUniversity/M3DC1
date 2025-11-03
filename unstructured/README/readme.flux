******************************************************************************
 NOTE: please add the following runtime options at your srun line:

 -sub_mat_superlu_dist_rowperm norowperm (when using superlu_dist for 3D run)
 -sub_mat_mumps_icntl_14 100 (when using mumps for 3D run)
******************************************************************************

0. documentation
   https://hpc.git.pppl.gov/flux/user-docs/#computing-resources

1. login 
	ssh your_id@flux.pppl.gov

	20 Nodes, each has 2x AMD EPYC 9354 32-Core, 256GB Memory.
	dual 25Gpbs ethernet.

2. code
   git clone https://your_id:password@github.com/PrincetonUniversity/M3DC1.git

3. load modules

   for release version modules:
       module use /home/nferraro/modules

   for development modules:
       setenv M3DC1_CODE_DIR $HOME/src/M3DC1  [for example]
       module use $M3DC1_CODE_DIR/unstructured/modules/flux
       module load m3dc1/devel

       or 
       module load intel/2023.2.0
       module load openmpi/4.1.6-intel-2023.2.0
       module load openmpi-4.1.6/intel-2023.2.0/netcdf-c/main-nzscm
       module load openmpi-4.1.6/intel-2023.2.0/netcdf-fortran/4.6.1-cepgo
       module load simmodsuite/2025.0-250217-dev
       module load netlib-lapack/3.11.0-intel-2023.2.0-pf3ms

4. compile code

   2D real: make OPT=1 RL=1 MAX_PTS=25 ARCH=flux
   2D complex: make OPT=1 COM=1 MAX_PTS=25 ARCH=flux
           - add PAR=1 to run PIC
   3D real: make 3D=1 OPT=1 MAX_PTS=60 ARCH=flux
   3D stellarator : make 3D=1 OPT=1 MAX_PTS=60 ST=1 ARCH=flux
   after compiling, run "make bin"

5. mesh generator
   module load intel/2023.2.0 openmpi/4.1.6-intel-2023.2.0
   module load simmodsuite/2025.0-250217-dev
   module load netlib-lapack/3.11.0-intel-2023.2.0-pf3ms

   /p/tsc/m3dc1/lib/SCORECLib/rhel9/intel2023.2.0-openmpi4.1.6/2025.0-250217-dev/bin/

   mesh utility (partioning, merging, etc.)
   /p/swim/jchen/PETSC/core-trunk/spark-20231214/bin

6. run a job 

   sample batch job scripts are in M3DC1/unstructured/regtest/*/base/batchjob.flux

   regression test

	setenv M3DC1_MPIRUN srun
	setenv M3DC1_VERSION local
	setenv M3DC1_ARCH flux
	make bin ARCH=$M3DC1_ARCH
	cd _flux/bin/
	setenv PATH `pwd`:$PATH
	cd ../../regtest/
	./run $M3DC1_ARCH
	./check $M3DC1_ARCH

7. solver options for mg in the toroidal direction

7.1 specify which solve is the mg target solve, for example, 

    use mg for solve #5
            srun -n 256 ./m3dc1_3d -options_file options_bjacobi.type -mgsolve 5

    use mg for solve #17
            srun -n 256 ./m3dc1_3d -options_file options_bjacobi.type -mgsolve 17

7.2 append the following lines in your solver option file: options_bjacobi.xxx

    -hard_ksp_type fgmres
    -hard_ksp_norm_type unpreconditioned
    -hard_pc_type mg
    -hard_ksp_rtol 1.e-03

    -hard_mg_levels_ksp_type gmres
    -hard_mg_levels_pc_type bjacobi
    -hard_mg_levels_sub_pc_type lu
    -hard_mg_levels_sub_pc_factor_mat_solver_type mumps
    -hard_mg_levels_sub_mat_mumps_icntl_14 100
    -hard_mg_levels_sub_ksp_type preonly

    -hard_mg_coarse_ksp_type gmres
    -hard_mg_coarse_pc_type bjacobi
    -hard_mg_coarse_sub_pc_type lu
    -hard_mg_coarse_sub_pc_factor_mat_solver_type mumps
    -hard_mg_coarse_sub_mat_mumps_icntl_14 100
    -hard_mg_coarse_sub_ksp_type preonly
