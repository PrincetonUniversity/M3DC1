  CPP = /projects/M3DC1/openmpi/traverse-nvidia/bin/mpic++
  CC = /projects/M3DC1/openmpi/traverse-nvidia/bin/mpicc
  F90 = /projects/M3DC1/openmpi/traverse-nvidia/bin/mpifort
  F77 = /projects/M3DC1/openmpi/traverse-nvidia/bin/mpifort
  LOADER = /projects/M3DC1/openmpi/traverse-nvidia/bin/mpifort

OPTS := $(OPTS) -DPETSC_VERSION=318 -DUSEBLAS

MPIVER=nvhpc21.5-openmpi4.0.3
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

SCOREC_BASE_DIR=/projects/M3DC1/scorec/traverse/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L/projects/M3DC1/scorec/traverse/$(MPIVER)/$(PETSCVER)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

HDF5_DIR=/projects/M3DC1/hdf5/traverse-nvidia
NETCDFDIR=/projects/M3DC1/netcdf/traverse-nvidia
GSL_DIR=/projects/M3DC1/gsl/traverse
FFTW_DIR=/projects/M3DC1/fftw/traverse-nvidia

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/projects/M3DC1/openmpi/traverse-nvidia/lib -L/opt/nvidia/hpc_sdk/Linux_ppc64le/21.5/compilers/lib -L/usr/lib/gcc/ppc64le-redhat-linux/8 -Wl,-rpath,/projects/M3DC1/openmpi/traverse-nvidia/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_ppc64le/21.5/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lflapack -lfblas -lparmetis -lmetis -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lnvf -lnvomp -latomic -lnvhpcatm -lpthread -lnvcpumath -lnvc -lrt -lm -lgcc_s -lstdc++ -ldl

#PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lsuperlu_dist -lsuperlu -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -llapack -lblas -lstdc++

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
       -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
       -I$(NETCDFDIR)/include \
       -I$(HDF5_DIR)/include \
       -I$(FFTW_DIR)/include \
       -I$(GSL_DIR)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        -L$(NETCDFDIR)/lib -lnetcdff -lnetcdf -lzip -lcurl\
        -L$(HDF5_DIR)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz\
        -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas \
        $(PETSC_WITH_EXTERNAL_LIB)

FOPTS = -c -r8 -Mpreprocess $(OPTS)
R8OPTS = -r8

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

ACC?=0
ifeq ($(ACC), 1)
  LDOPTS := $(LDOPTS) -acc -gpu=cuda11.3 -Minfo=accel
  FOPTS  := $(FOPTS)  -acc -gpu=cuda11.3 -Minfo=accel
  CCOPTS  := $(CCOPTS) -acc -gpu=cuda11.3 -Minfo=accel
endif

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

