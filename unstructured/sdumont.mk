ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpiicpc 
  CC = mpiicc 
  F90 = mpiifort 
  F77 = mpiifort 
  LOADER = mpiifort 
endif

#NEWSOLVERDEVELOPMENT needs more tests.
OPTS := $(OPTS) -DPETSC_VERSION=37 #-DNEWSOLVERDEVELOPMENT

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB = m3dc1_scorec_complex
  PETSC_DIR=
  PETSC_ARCH=
  HYPRE_LIB = 
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif
  ifeq ($(OMP), 1)
    PETSC_DIR=
    PETSC_ARCH=
  else
    PETSC_DIR=/scratch/ntm/software/petsc-3.7.6
    PETSC_ARCH=real-intel
  endif
  HYPRE_LIB = -lHYPRE
endif

PETSC_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
     -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
     $(HYPRE_LIB) \
     -lsuperlu -lsuperlu_dist \
     -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl -lstdc++


SCOREC_UTIL_DIR=
SCOREC_DIR=/scratch/ntm/software/scorec/3.7.6/May2018
ZOLTAN_DIR=/scratch/ntm/software/scorec

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

# Include option to use ADIOS
#OPTS := $(OPTS) -DUSEADIOS
#ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
#ADIOS_FLIB_V1 = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
#             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
             -L/usr/lib64/ -llustreapi

MKL_LIB = $(MKLROOT)/lib/intel64_lin/libmkl_blas95_lp64.a -L$(MKLROOT)/lib/intel64_lin -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \

LIBS := \
        $(LIBS) \
        $(SCOREC_LIBS) \
        -L$(ZOLTAN_DIR)/lib -lzoltan \
        $(PETSC_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB_V1) \
	$(MKL_LIB)

# Optimization flags
FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)
CCOPTS  = -c $(OPTS)
ifeq ($(OPT), 0)
  FOPTS := $(FOPTS) -g -O0 -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -openmp 
  FOPTS  := $(FOPTS)  -openmp 
  CCOPTS := $(CCOPTS) -openmp 
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
