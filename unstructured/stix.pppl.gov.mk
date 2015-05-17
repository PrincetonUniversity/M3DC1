FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DKSPITS $(OPTS)
CCOPTS  = -c -O -DPetscDEV -DKSPITS -DPetscOLD #-DCJ_MATRIX_DUMP -DUSEHYBRID 

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -vec-report0 -fast
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


# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/pdslin_0.0
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

INCLUDE = \
        -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SUPERLU_HOME)/include -I$(SUPERLU_DIST_HOME)/include \
	-I$(SCALAPACK_HOME)/include -I$(BLACS_HOME)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib \
	-I$(FFTWHOME)/include \
	-I$(HYBRID_HOME)/include \
	-I$(GSLHOME)/include

#        -L$(PETSC_DIR)/lib -lpetsc -ldmumps -lmumps_common -lcmumps -lpord \

ifeq ($(COM), 1)
PETSC_DIR= /usr/pppl/intel/11-pkgs/vSMPICH2-pkgs/petsc-3.4.5/
PETSC_ARCH = intel-vsmp2-complex
PETSC_LIBS = \
        -Wl,--start-group,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lmumps_common -lpetsc -lpord lsmumps -lzmumps -Wl,--end-group \
	-L$(SCALAPACK_HOME) -lscalapack \
	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacs -lmpiblacsCinit -lmpiblacs
SCORECLIBS = -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lm3dc1_scorec_complex -lapf_zoltan 
else
PETSC_DIR= /usr/pppl/intel/11-pkgs/vSMPICH2-pkgs/petsc-3.4.5-real/
PETSC_ARCH= intel-vsmp2-real
PETSC_LIBS = \
        -Wl,--start-group,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lmumps_common -lpetsc -lpord lsmumps -lzmumps -Wl,--end-group \
	-L$(SCALAPACK_HOME) -lscalapack \
	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacs -lmpiblacsCinit -lmpiblacs
SCORECLIBS = -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lm3dc1_scorec -lapf_zoltan
endif

#SUPERLU_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
#SUPERLU_DIST_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu_4.3 \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_3.3 \

SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/stix/May2015
INCLUDE := -I$(SCORECDIR)/include $(INCLUDE)

SCOREC_LIBS = -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib \
              $(SCORECLIBS) -Wl,--start-group

LIBS = 	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
        $(SCOREC_LIBS) \
        -L/usr/pppl/intel/11-pkgs/vSMPICH2-pkgs/ParMetis-4.0.3/lib -lmetis -lparmetis \
        -L/usr/pppl/intel/11-pkgs/vSMPICH2-pkgs/zoltan-3.81/lib -lzoltan \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-L$(FFTWHOME)/lib -lfftw3 \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack \
	-L$(CCHOME)/lib/intel64 -lguide \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
	-L$(GSLHOME)/lib -lgsl \
	-L/usr/X11R6/lib -lX11

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
