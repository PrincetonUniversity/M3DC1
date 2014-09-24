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


# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/pdslin_0.0
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib \
	-I$(HYBRID_HOME)/include \
	-I$(GSLHOME)/include

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetsc \
	-lpromfei -lprometheus \
	-lHYPRE \
	-L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord \
	-L$(SCALAPACK_HOME) -lscalapack \
	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacs -lmpiblacsCinit -lmpiblacs

SUPERLU_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
SUPERLU_DIST_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu_4.1 \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_2.5 \

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis


LIBS = 	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(PARMETIS_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-L$(FFTWHOME)/lib -lfftw3 \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
	-L$(GSLHOME)/lib -lgsl \
	-L/usr/X11R6/lib -lX11

ifeq ($(USESCOREC), 1)

#  SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/latest/
  SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/04082013/
  INCLUDE := -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/093011 \
        $(INCLUDE)


  SCOREC_ARCH=x86_64_linux-icc
  SCOREC_LIBS = \
	-L$(SCORECDIR) \
	-Wl,-rpath,$(SCORECDIR) \
        -lPPPLFusion \
        -lMeshAdapt -lFMDB -lGMI -lGMIMeshModel -lSCORECUtil -lipcomman   

#  SCOREC_LIBS = \
#	-L$(SCORECDIR) \
#	-Wl,-rpath,$(SCORECDIR) \
#	-lFMDB-mpich2$(SCORECOPT) \
#	-lSCORECModel-mpich2$(SCORECOPT) \
#	-lSCORECUtil-mpich2$(SCORECOPT) \
#	-lField-mpich2$(SCORECOPT) \
#	-lCore-mpich2$(SCORECOPT) \
#	-lmeshAdapt-mpich2$(SCORECOPT) \
#	-lmeshTools-mpich2$(SCORECOPT) \
#	-lSolver-mpich2$(SCORECOPT) \
#	-lPPPL-mpich2$(SCORECOPT) \
#	-lipcomman-mpich2$(SCORECOPT)
##	-lPPPLPetscDEV-mpich2$(SCORECOPT) \

  LIBS := $(SCOREC_LIBS) $(LIBS)

endif   # on USESCOREC

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
