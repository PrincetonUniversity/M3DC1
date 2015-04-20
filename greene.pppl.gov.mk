#FOPTS = -c -r8 -implicitnone -fpp -warn all -DxPetscDEV -DPETSC_31 $(OPTS)
FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DKSPITS -DxCJ_MATRIX_DUMP  $(OPTS) #-pg
CCOPTS  = -c -O -DPetscDEV -DKSPITS -DxCJ_MATRIX_DUMP -DxUSEHYBRID #-DPETSC_31 -pg

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
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90 -cxxlib #-pg
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)


# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/pdslin_0.0
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

INCLUDE = -I$(MPIHOME)/include \
	  -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	  -I$(PARMETIS_HOME)/include \
	  -I$(FFTWHOME)/include \
	  -I$(SUPERLU_DIST_HOME)/include -I$(SUPERLUHOME)/include \
	  -I$(BLACS_HOME)/include \
	  -I$(Zoltan_HOME)/include \
	  -I$(GSLHOME)/include \
	  -I$(HDF5_HOME)/include \
	  -I$(NCARG_ROOT)/include

#	-L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord \
#	-L$(SCALAPACK_HOME) -lscalapack \
#	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacs -lmpiblacsCinit -lmpiblacs

LIBS = 	-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lmumps_common -ldmumps -lcmumps -lzmumps -lsmumps -lpord -lpetsc \
        -L$(PARMETIS_HOME)/lib -lparmetis -lmetis \
	-L$(FFTWHOME)/lib -lfftw3 -lfftw3_mpi -lfftw3_threads \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_3.3 -L$(SUPERLUHOME)/lib -lsuperlu_4.3 \
	-L$(SCALAPACK_HOME)/lib -lscalapack -L$(BLACS_HOME)/lib -lmpiblacs -lmpiblacsCinit -lmpiblacsF77init \
        -L$(Zoltan_HOME)/lib -lzoltan \
	-L$(GSLHOME)/lib -lgsl \
        -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
        -L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L/usr/X11R6/lib -lX11

  SCORECDIR=/p/tsc/m3dc1/lib/SCORECLib/greene/Mar2015
ifeq ($(COM), 1)
  SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib \
             -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lm3dc1_scorec_complex \
             -Wl,--start-group
else
  SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib \
             -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lm3dc1_scorec \
             -Wl,--start-group
endif

  LIBS := $(SCORECLIB) \
	$(LIBS) 
  INCLUDE := -I$(SCORECDIR)/include \
        $(INCLUDE)

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
