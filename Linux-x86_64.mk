LOADER = mpif90
F90    = mpif90
F77    = mpif90
CC     = mpicc

# define where you want to locate the mesh adapt libraries

INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib

H5_VERSION = 169

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist

MUMPS_LIBS = -L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord

BLACS_LIBS = -L$(BLACS_HOME)/lib -lmpiblacs -lmpiblacsF77init

SCALAPACK_LIBS = -L$(SCALAPACK_HOME)/lib -lscalapack

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

NAG_LIBS = -L$(NAG_ROOT)/lib -lnag

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib -lautopack-O

LIBS = 	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(MUMPS_LIBS) \
	$(SCALAPACK_LIBS) \
	$(BLACS_LIBS) \
	$(PARMETIS_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack -lmkl_ipr \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
        -L/usr/X11R6/lib -lX11


ifeq ($(USESCOREC), 1)

  ifeq ($(USE3D), 1)

    # 3D libraries
    ifndef SCORECDIR
      SCORECDIR = /p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.test/lib
    endif
    INCLUDE := -I/p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.test/include $(INCLUDE)

    SCOREC_ARCH=x86_64_linux-icc

    SCOREC_LIBS = -L$(SCORECDIR) \
	-Wl,-rpath,$(SCORECDIR) \
	-lFMDB-mpich2$(SCORECOPT) \
	-lSCORECModel-mpich2$(SCORECOPT) \
	-lSCORECUtil-mpich2$(SCORECOPT) \
	-lField-mpich2$(SCORECOPT) \
	-lCore-mpich2$(SCORECOPT) \
	-lmeshAdapt-mpich2$(SCORECOPT) \
	-ltemplateRefine-mpich2$(SCORECOPT) \
	-lmeshTools-mpich2$(SCORECOPT) \
	-lSolver-mpich2$(SCORECOPT) \
	-lPPPL-mpich2$(SCORECOPT) \
	-lipcomman-mpich2$(SCORECOPT)

#    SCOREC_LIBS = \
#	-L$(SCORECDIR)/FMDB/FMDB/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/FMDB/FMDB/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/mctk/Field/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/mctk/Field/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/mctk/Core/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/mctk/Core/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/mctk/Solver/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/mctk/Solver/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)/ipcomman/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)/ipcomman/lib/$(SCOREC_ARCH) \
#	-lFMDB-mpich2$(SCORECOPT) \
#	-lSCORECModel-mpich2$(SCORECOPT) \
#	-lSCORECUtil-mpich2$(SCORECOPT) \
#	-lField-mpich2$(SCORECOPT) \
#	-lCore-mpich2$(SCORECOPT) \
#	-lmeshAdapt-mpich2$(SCORECOPT) \
#	-ltemplateRefine-mpich2$(SCORECOPT) \
#	-lmeshTools-mpich2$(SCORECOPT) \
#	-lSolver-mpich2$(SCORECOPT) \
#	-lPPPL-mpich2$(SCORECOPT) \
#	-lipcomman-mpich2$(SCORECOPT)
  else
    # 2D Libraries

    ifndef SCORECDIR
      SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/092210
#      SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/10212010
#      SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/latest
    endif
    INCLUDE := -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/092210 $(INCLUDE)
#    INCLUDE := -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/10212010 $(INCLUDE)
#    INCLUDE := -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/latest $(INCLUDE)

    SCOREC_LIBS = \
	-L$(SCORECDIR) \
	-Wl,-rpath,$(SCORECDIR) \
	-lFMDB-mpich2$(SCORECOPT) \
	-lSCORECModel-mpich2$(SCORECOPT) \
	-lSCORECUtil-mpich2$(SCORECOPT) \
	-lField-mpich2$(SCORECOPT) \
	-lCore-mpich2$(SCORECOPT) \
	-lmeshAdapt-mpich2$(SCORECOPT) \
	-ltemplateRefine-mpich2$(SCORECOPT) \
	-lmeshTools-mpich2$(SCORECOPT) \
	-lSolver-mpich2$(SCORECOPT) \
	-lPPPL-mpich2$(SCORECOPT)

  endif # on USE3D

  LIBS := $(SCOREC_LIBS) $(AUTOPACK_LIBS) $(LIBS)

endif   # on USESCOREC


FOPTS = -c -r8 -implicitnone -fpp -warn all $(INCLUDE) $(OPTS) \
	-DH5_VERSION=$(H5_VERSION) -DRANDOM_NUM='drand(0)' \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs #\
#	-g -check all -check noarg_temp_created -debug all -ftrapuv
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)
CCOPTS = -c $(INCLUDE)


%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) -fpic $< -o $@
