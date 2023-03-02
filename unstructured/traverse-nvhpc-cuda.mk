  CPP = mpic++
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort

OPTS := $(OPTS) -DPETSC_VERSION=318 -DUSEBLAS

MPIVER=nvhpc22.5-openmpi4.1.4
PETSC_VER=petsc-3.18.3
PETSCVER=petsc3.18.3

PETSC_DIR=/projects/M3DC1/scorec/petsc/$(PETSC_VER)

ifeq ($(COM), 1)
   PETSC_ARCH=cplx-$(MPIVER)
else
   PETSC_ARCH=real-$(MPIVER)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_BASE_DIR==/projects/M3DC1/scorec/traverse/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif

GSL_DIR=$(GSL_ROOT_DIR)
HDF5_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/usr/local/openmpi/cuda-11.7/4.1.4/nvhpc225/ppc64le/lib64 -L/usr/local/nvhpc/lib64 -L/usr/local/nvhpc/lib64/openmpi -L/usr/local/hdf5/nvhpc-22.5/openmpi-4.1.3/1.10.6/lib64 -L/opt/nvidia/hpc_sdk/Linux_ppc64le/22.5/compilers/lib -L/usr/lib/gcc/ppc64le-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/cuda-11.7/4.1.4/nvhpc225/ppc64le/lib64 -Wl,-rpath,/usr/local/nvhpc/lib64 -Wl,-rpath,/usr/local/nvhpc/lib64/openmpi -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_ppc64le/22.5/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lnvf -lnvomp -latomic -lnvhpcatm -lpthread -lnvcpumath -lnvc -lrt -lgcc_s -lm -lstdc++ -ldl

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
       -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
       -I$(HDF5_DIR)/include \
       -I$(GSL_DIR)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        -L$(GSL_DIR)/lib64 -lgsl -lgslcblas \
	-L$(HDF5_DIR)/lib64 -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz \
        $(PETSC_WITH_EXTERNAL_LIB)

FOPTS = -c -r8 -Mpreprocess $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -fast
  FOPTS  := $(FOPTS)  -fast
  CCOPTS := $(CCOPTS) -fast
else
  FOPTS := $(FOPTS) -g
  CCOPTS := $(CCOPTS) -g
  LDOPTS := $(LDOPTS) -g
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

#ACC?=1
#ifeq ($(ACC), 1)
#  LDOPTS := $(LDOPTS) -acc -gpu=cuda11.3 -Minfo=accel
#  FOPTS  := $(FOPTS)  -acc -gpu=cuda11.3 -Minfo=accel
#  CCOPTS  := $(CCOPTS) -acc -gpu=cuda11.3 -Minfo=accel
#endif

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@

