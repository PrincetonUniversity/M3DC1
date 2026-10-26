FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=393 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -DPETSC_VERSION=393
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -w -fallow-argument-mismatch -O2 -ffree-line-length-none
# FOPTS  := $(FOPTS) -ggdb3 -Og -w -fallow-argument-mismatch
  CCOPTS := $(CCOPTS) -ggdb3 -Og
else
  FOPTS := $(FOPTS) -g #noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

CC = cc
CPP = CC
F90 = ftn
F77 = ftn
LOADER = ftn
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

#  PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606
  PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609
ifeq ($(COM), 1)
#  PETSC_ARCH=perlmuttergpu-gnu-cplx
#  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu-cplx/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64/stubs -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lkokkoskernels -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lcudart -lnvToolsExt -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lquadmath -lmpifort_gnu_123 -lgfortran -lstdc++
  PETSC_ARCH=perlmuttergpu-gnu-870-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttergpu-gnu-870-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttergpu-gnu-870-cplx/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/math_libs/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/math_libs/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64/stubs -Wl,-rpath,/opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3/lib -L/opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/14 -L/usr/lib64/gcc/x86_64-suse-linux/14 -Wl,-rpath,/usr/x86_64-suse-linux/lib -L/usr/x86_64-suse-linux/lib -lpetsc -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lkokkoskernels -lsuperlu -lflapack -lfblas -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lcudart -lnvtx3interop -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lmpifort_gnu -lmpi_gnu -lxpmem -lgfortran -lm -lhugetlbfs -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath

else
#  PETSC_ARCH=perlmuttergpu-gnu
#  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64/stubs -lpetsc -lHYPRE -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lkokkoskernels -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lcudart -lnvToolsExt -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lquadmath -lmpifort_gnu_123 -lgfortran -lstdc++
  PETSC_ARCH=perlmuttergpu-gnu-870
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttergpu-gnu-870/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttergpu-gnu-870/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/math_libs/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/math_libs/13.2/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/26.5/cuda/13.2/lib64/stubs -Wl,-rpath,/opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3/lib -L/opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/14 -L/usr/lib64/gcc/x86_64-suse-linux/14 -Wl,-rpath,/usr/x86_64-suse-linux/lib -L/usr/x86_64-suse-linux/lib -lpetsc -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lkokkoskernels -lsuperlu -lflapack -lfblas -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lcudart -lnvtx3interop -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lmpifort_gnu -lmpi_gnu -lxpmem -lgfortran -lm -lhugetlbfs -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath

#  /opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64/libcudart.so.12 
endif

#SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/core-240527/upgrade-gnu850-pgpu
SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core-240527/upgrade-gnu870-pgpu
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  ifeq ($(CUDSS_MG), 1)
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_mg_complex
    M3DC1_CUDSS_LIB=-lm3dc1_cudss_mg_complex
  else
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
    M3DC1_CUDSS_LIB=-lm3dc1_cudss_complex
  endif
else
  ifeq ($(CUDSS_MG), 1)
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_mg
    M3DC1_CUDSS_LIB=-lm3dc1_cudss_mg
  else
    M3DC1_SCOREC_LIB=-lm3dc1_scorec
    M3DC1_CUDSS_LIB=-lm3dc1_cudss
  endif
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

# cuDSS block-Jacobi solver (matrix_solve::solve_cudss, runtime option
# -cudsssolve <matrix_id>; CUDSS_MG=1 for the -cudsssolve_mg variant).
# Requires m3dc1_scorec built with ENABLE_CUDSS=ON (or ENABLE_CUDSS_MG=ON) and
# the matching libm3dc1_cudss*.a from m3d-c1_preconditioner_standalone,
# installed alongside m3dc1_scorec.  Vendored cuDSS ships only libcudss.so.0,
# so link it by full path; NCCL_DIR comes from `module load nccl`.
# cudart/cusparse/cublas are already linked via PETSc.
#ifeq ($(CUDSS), 1)
  CUDSS_PKG_DIR ?= $(M3DC1_DIR)/m3d-c1_preconditioner_standalone
  CUDSS_DIR ?= $(CUDSS_PKG_DIR)/deps/cudss.13.2
  CUDSS_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_CUDSS_LIB) \
              -Wl,-rpath,$(CUDSS_DIR)/lib $(CUDSS_DIR)/lib/libcudss.so.0 \
              -Wl,-rpath,$(NCCL_DIR)/lib -L$(NCCL_DIR)/lib -lnccl
#else
#  CUDSS_LIB =
#endif

LIBS = 	\
	-L$(HDF5_DIR)/lib -lhdf5_parallel -lhdf5_hl_parallel -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel \
	$(SCOREC_LIB) \
	$(CUDSS_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \

#	/opt/cray/pe/mpich/9.1.0/gtl/lib/libmpi_gtl_cuda.so.0
#	/opt/cray/pe/mpich/8.1.30/gtl/lib/libmpi_gtl_cuda.so.0
#	/opt/cray/pe/lib64/libmpi_gtl_cuda.so.0

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include -I$(HDF5_DIR)/include 

ifeq ($(ST), 1)
  LIBS += -L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.13/gnu/12.3/lib -lnetcdf -lnetcdff
  INCLUDE += -I/opt/cray/pe/netcdf-hdf5parallel/4.9.0.13/gnu/12.3/include
endif

#ACC?=1
ifeq ($(ACC), 1)
  LDOPTS := $(LDOPTS) #-acc -gpu=cuda11.7 -Minfo=accel
  FOPTS  := $(FOPTS)  #-acc -gpu=cuda11.7 -Minfo=accel
  CCOPTS  := $(CCOPTS) #-acc -gpu=cuda11.7 -Minfo=accel
endif

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP) $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
