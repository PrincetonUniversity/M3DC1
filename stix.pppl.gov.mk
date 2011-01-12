H5_VERSION = 169

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-DH5_VERSION=$(H5_VERSION) 
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs #\
#	-g -check all -check noarg_temp_created -debug all -ftrapuv
CCOPTS  = -c

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CC     = tau_cc.sh $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CC = mpicc -c
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90
  FOPTS := $(FOPTS) -DRANDOM_NUM='drand(0)'
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)


# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test

INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib \
	-I$(HYBRID_HOME)/include -I$(SUPERLU_DIST_HOME)/include


PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist

MUMPS_LIBS = -L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord

BLACS_LIBS = -L$(BLACS_HOME)/lib -lmpiblacs -lmpiblacsF77init -lmpiblacsCinit -lmpiblacs

SCALAPACK_LIBS = -L$(SCALAPACK_HOME)/lib -lscalapack

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

NAG_LIBS = -L$(NAG_ROOT)/lib -lnag

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib -lautopack-O

HYBRID_LIBS = $(HYBRID_HOME)/lib/libhsolver.a

LIBS = 	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(MUMPS_LIBS) \
	$(SCALAPACK_LIBS) \
	$(BLACS_LIBS) \
	$(PARMETIS_LIBS) \
	$(HYBRID_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack -lmkl_ipr \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
        -L/usr/X11R6/lib -lX11


ifeq ($(USESCOREC), 1)

#  ifeq ($(USE3D), 1)

    # 3D libraries
    ifndef SCORECDIR
      SCORECDIR = /p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.test/lib
    endif
    INCLUDE := -I/p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.test/include \
	$(INCLUDE)

    SCOREC_ARCH=x86_64_linux-icc
    SCOREC_LIBS = \
        -L$(SCORECDIR) \
	-Wl,-rpath,$(SCORECDIR) \
	-lFMDB-mpich2$(SCORECOPT) \
	-lSCORECModel-mpich2$(SCORECOPT) \
	-lSCORECUtil-mpich2$(SCORECOPT) \
	-lField-mpich2$(SCORECOPT) \
	-lCore-mpich2$(SCORECOPT) \
	-lmeshAdapt-mpich2$(SCORECOPT) \
	-lmeshTools-mpich2$(SCORECOPT) \
	-lSolver-mpich2$(SCORECOPT) \
	-lPPPL-mpich2$(SCORECOPT) \
	-lipcomman-mpich2$(SCORECOPT)
#  else
#    # 2D Libraries
#
#    ifndef SCORECDIR
#      SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/092210
#    endif
#    INCLUDE := -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/092210 $(INCLUDE)
#
#    SCOREC_LIBS = \
#        -L$(SCORECDIR) \
#	-Wl,-rpath,$(SCORECDIR) \
#	-lFMDB-mpich2$(SCORECOPT) \
#	-lSCORECModel-mpich2$(SCORECOPT) \
#	-lSCORECUtil-mpich2$(SCORECOPT) \
#	-lField-mpich2$(SCORECOPT) \
#	-lCore-mpich2$(SCORECOPT) \
#	-lmeshAdapt-mpich2$(SCORECOPT) \
#	-ltemplateRefine-mpich2$(SCORECOPT) \
#	-lmeshTools-mpich2$(SCORECOPT) \
#	-lSolver-mpich2$(SCORECOPT) \
#	-lPPPL-mpich2$(SCORECOPT)
#
#  endif # on USE3D

  LIBS := $(SCOREC_LIBS) $(AUTOPACK_LIBS) $(LIBS)

endif   # on USESCOREC


%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
