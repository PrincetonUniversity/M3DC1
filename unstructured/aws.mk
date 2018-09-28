CCOPTS := -c $(OPTS) -DPETSC_VERSION=38
FOPTS := -c $(OPTS) -cpp -fdefault-real-8 -Wall -DUSEBLAS -DPETSC_VERSION=38

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O3
  CCOPTS := $(CCOPTS) -O3
else
  FOPTS := $(FOPTS) -g
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=../select.tau
  CC     = tau_cc.sh $(TAU_OPTIONS)
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CC = mpicc
  CPP = mpicxx
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90
  LDOPTS := $(LDOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=
ifeq ($(COM), 1)
PETSC_ARCH=complex
HYPRE_LIB=
else
PETSC_ARCH=real
HYPRE_LIB=
endif

BLASLAPACK_LIBS = -lfblas -lflapack

HDF5_DIR=/home/ec2-user/hdf5-1.10.2/hdf5
SCOREC_DIR=$(HOME)/scorec/June2018
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_UTIL_DIR=

SCALAPACK_LIB=-Wl,-rpath,$(SCALAPACK_HOME)/lib -L$(SCALAPACK_HOME) -lscalapack

PETSC_DIR=/home/ec2-user/petsc-3.8.4
ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_ARCH=aws_real
endif

PETSC_LIBS = -Wl,--start-group,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetsc \
        -lsuperlu_dist \
        -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
        $(SCALAPACK_LIB) \
        $(HYPRE_LIB) \
	-lzoltan -lptscotch \
        -lparmetis -lmetis \
        -Wl,--end-group

SCORECLIB= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group

HDF5_LIBS = -Wl,--start-group,-rpath,$(HDF5_DIR)/lib \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz \
	-Wl,--end-group

LIBS = 	\
	$(SCORECLIB) \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	$(HDF5_LIBS) \
	-lfftw3 \
	-lgsl -lgslcblas -lstdc++


INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_DIR)/include

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP) $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
