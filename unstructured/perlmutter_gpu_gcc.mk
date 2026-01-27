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

  PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606
ifeq ($(COM), 1)
  PETSC_ARCH=perlmuttergpu-gnu-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu-cplx/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64/stubs -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lkokkoskernels -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lcudart -lnvToolsExt -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lquadmath -lmpifort_gnu_123 -lgfortran -lstdc++
else
  PETSC_ARCH=perlmuttergpu-gnu
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttergpu-gnu/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64/stubs -lpetsc -lHYPRE -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lkokkoskernels -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lcudart -lnvToolsExt -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lquadmath -lmpifort_gnu_123 -lgfortran -lstdc++
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/core-240527/upgrade-gnu850-pgpu
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

LIBS = 	\
	-L$(HDF5_DIR)/lib -lhdf5_parallel -lhdf5_hl_parallel -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel \
	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	/opt/cray/pe/lib64/libmpi_gtl_cuda.so.0

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
