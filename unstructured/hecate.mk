FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DKSPITS $(OPTS)
CCOPTS  = -c -O -DPetscDEV -DKSPITS -DPetscOLD #-DCJ_MATRIX_DUMP -DUSEHYBRID 

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -vec-report0 # -fast
#  FOPTS  := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
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
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90 -cxxlib
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)


# define petsc and other libs
	PETSC_DIR = /home/jinchen/lib/petsc-3.5.3
ifeq ($(COM), 1)
	PETSC_ARCH = hecate-intel-mt-complex
else
	PETSC_ARCH = hecate-intel-mt
endif

INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5DIR)/include -I$(HDF5DIR)/lib64 \

ifeq ($(COM), 1)
PETSCLIBS = \
        -lpetsc \
        -lzoltan \
        -lfftw3_mpi -lfftw3 \
        -lmetis -lparmetis \
        -lsuperlu_4.3 -lsuperlu_dist_3.3 \
        -ldmumps -lmumps_common -lcmumps -lzmumps -lpord
else
PETSCLIBS = \
        -lpetsc \
        -lHYPRE \
        -lzoltan \
        -lfftw3_mpi -lfftw3 \
        -lmetis -lparmetis \
        -lsuperlu_4.3 -lsuperlu_dist_3.3 \
        -ldmumps -lmumps_common -lcmumps -lzmumps -lpord
endif

BLASLAPACKLIBS = -L$(MKLROOT)/lib/intel64 -Wl,--start-group -lmkl_blacs_sgimpt_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_cdft_core -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group $(PETSCLIBS) -Wl,--end-group

SCORECDIR = /home/jinchen/lib/scorec/Dec2015
INCLUDE := -I$(SCORECDIR)/include $(INCLUDE)
ifeq ($(COM), 1)
  SCORECLIB=-lapf -lgmi -lm3dc1_scorec_complex -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr
else
  SCORECLIB=-lapf -lgmi -lm3dc1_scorec -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr
endif
SCOREC_LIBS =-L$(SCORECDIR)/lib -Wl,--start-group $(SCORECLIB) -Wl,--end-group

LIBS = \
	$(BLASLAPACKLIBS) \
       $(SCOREC_LIBS) \
       $(PETSC_LIBS) \
        -L$(HDF5DIR)/lib64 -lhdf5_fortran -lhdf5 \
        -Wl,-rpath -Wl,$(HDF5_HOME)/lib \
        -L$(ZLIB_HOME) -lz \
        -L/usr/X11R6/lib -lX11 \
        -L/usr/lib64 -lssl -lgsl -lgslcblas

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
