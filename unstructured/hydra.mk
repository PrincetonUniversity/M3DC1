FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)
CCOPTS  = -c -O

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -fast
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
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
  F90 = mpiifort
  F77 = mpiifort
  LOADER = mpiifort -cxxlib
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

SUPERLU_DIST_HOME = $(PETSC_DIR)
INCLUDE = -I$(MPIHOME)/include \
        -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SUPERLU_DIST_HOME)/include \
        -I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib \
        -I$(GSL_HOME)/include \
        -I$(FFTW_HOME)/include

PETSC_LIBS = -L$(PETSC_DIR)/lib \
   -lpetsc \

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu_4.3 \
        -L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_3.3 \

PARMETIS_HOME=$(PETSC_DIR)
PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
        -Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

SCORECDIR = /hydra/u/m3dc1/scorec/Dec2015/lib
INCLUDE := -I/hydra/u/m3dc1/scorec/Dec2015/include \
        $(INCLUDE)

ifeq ($(COM), 1)
  SCOREC_LIBS =-L$(SCORECDIR) -Wl,--start-group -lapf -lgmi -lmds -lpcu \
              -lspr -lapf_zoltan -lma -lparma -lph -lm3dc1_scorec_complex \
              -Wl,--end-group -lzoltan
else
  SCOREC_LIBS =-L$(SCORECDIR) -Wl,--start-group -lapf -lgmi -lmds -lpcu \
              -lspr -lapf_zoltan -lma -lparma -lph -lm3dc1_scorec \
              -Wl,--end-group -lzoltan
endif

LIBS =  $(PETSC_LIBS) \
        $(SUPERLU_LIBS) \
        $(SCOREC_LIBS) \
        $(PARMETIS_LIBS) \
        -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
        -L$(FFTW_HOME)/lib -lfftw3 -lfftw3_mpi -lfftw3_threads \
        -L$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 \
        -Wl,-rpath -Wl,$(HDF5_HOME)/lib \
        -L$(ZLIB_HOME) -lz \
        -L$(GSL_HOME)/lib -lgsl \
        -L/usr/X11R6/lib -lX11

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
