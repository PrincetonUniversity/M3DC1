FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DLATESTSCOREC -DUSEBLAS -DPETSC_VERSION=39
CCOPTS  = -c -DPETSC_VERSION=39

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -g -check-uninit -debug all
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
  LDOPTS := $(LDOPTS) -cxxlib
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
PETSC_DIR=/home/jinchen/LIB/petsc-3.9.3
ifeq ($(COM), 1)
PETSC_ARCH=eddy-intel-openmpi-cplx
HYPRE_LIB=
else
PETSC_ARCH=eddy-intel-openmpi-real
HYPRE_LIB=
endif
BLASLAPACK_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
BLASLAPACK_LIBS =-Wl,-rpath,$(BLASLAPACK_DIR)/lib -L$(BLASLAPACK_DIR)/lib -lflapack -lfblas

SCALAPACK_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
SCALAPACK_LIB=-Wl,-rpath,$(SCALAPACK_DIR)/lib -L$(SCALAPACK_DIR) -lscalapack

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
      -lpetsc \
      -lptscotcherr -lscotch -lptscotcherrexit -lscotcherr -lptesmumps -lscotcherrexit -lptscotch \
      -lsuperlu_dist -lsuperlu \
      -lcmumps -ldmumps -lsmumps -lesmumps -lzmumps -lmumps_common -lpord \
      $(SCALAPACK_LIB) \
      -lfftw3 -lfftw3_mpi \
      $(HYPRE_LIB) \
      -lparmetis -lmetis \
      -Wl,--end-group

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB=-lm3dc1_scorec
  endif
endif

SCOREC_BASE_DIR=/home/jinchen/LIB/scorec/intel16.0-openmpi1.10.2/petsc-3.9.3
SCOREC_UTIL_DIR=/home/jinchen/LIB/scorec/intel16.0-openmpi1.10.2/bin
SCOREC_DIR=/home/jinchen/LIB/m3dc1_scorec/build

ZOLTAN_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

LIBS = 	\
        -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
	$(SCOREC_LIBS) \
        $(TRILINOS_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5DIR)/lib64 -lhdf5_fortran -lhdf5 -lz \
	-L$(GSL_ROOT_DIR)/lib64 -lgsl -lgslcblas \
	-lX11

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5DIR)/include \
        -I$(GSL_ROOT_DIR)/include

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
