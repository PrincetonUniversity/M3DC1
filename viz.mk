LOADER = ifort -Wl,--warn-unresolved-symbols
F90    = ifort
F77    = ifort
CC     = icc

HYBRID_HOME = /p/swim/jchen/hybrid.test

INCLUDE = -I$(NTCCHOME)/mod -I$(LIBDIR) \
	-I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HYBRID_HOME)/include

FOPTS = -c -r8 -implicitnone -fpp -warn all $(INCLUDE) $(OPTS) -O #\
#	-g -check all -check noarg_temp_created 

F90OPTS = $(FOPTS)
F77OPTS = $(FOPTS)
CCOPTS = -c $(INCLUDE)

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib/ia64-sgi \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib/ia64-sgi -lautopack-O

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc \

SUPERLU_LIBS = -L$(SUPERLU_HOME) -lsuperlu_3.0 \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

LIBS = 	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(PARMETIS_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-L$(CCHOME)/lib -lipr -lstdc++ \
	-L$(MKLHOME)/lib/64 -lguide -lmkl_lapack -lmkl_ipf \
	-L$(LIBDIR) -lhdf5_fortran -lhdf5 -lz \
	-Wl,-rpath -Wl,$(LIBDIR) \
        -L/usr/X11R6/lib -lX11 -lmpi


ifeq ($(USESCOREC), 1)

  SCOREC_ARCH = ia64_linux

  ifndef SCORECDIR
    SCORECDIR = /p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.viz
  endif

  SCORECDIR = /p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.viz/lib
  INCLUDE := -I/p/tsc/m3dc1/lib/develop.petsc3.Fan/develop.viz/include \
             $(INCLUDE)

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


#ifeq ($(USE3D), 1)
#
## For new libraries =============================
#ifndef SCORECDIR
#  SCORECDIR = /p/tsc/m3dc1/lib/develop.petscGlob2/
#endif
#FOPTS := $(FOPTS) -DUSERW
#SCORECINCLUDE = -I$(SCORECDIR)/mctk/Examples/PPPL/PPPL
#SCOREC_LIBS = \
#	-L$(SCORECDIR)FMDB/FMDB/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)FMDB/FMDB/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)mctk/Field/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)mctk/Field/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)mctk/Core/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)mctk/Core/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)mctk/Solver/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)mctk/Solver/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
#	-L$(SCORECDIR)meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
#	-Wl,-rpath,$(SCORECDIR)meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
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
## ===============================================
#
#else  # on 3D
#
## For old libraries =============================
#ifndef SCORECDIR
#  SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Viz/022310
#endif
#SCORECINCLUDE = -I/p/tsc/m3dc1/lib/SCORECLib/include/Viz/022310
#SCOREC_LIBS = \
#	-L$(SCORECDIR) \
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
## ===============================================
#
#endif  # on 3D

INCLUDE := $(SCORECINCLUDE) $(INCLUDE)
LIBS := $(SCOREC_LIBS) $(AUTOPACK_LIBS) $(LIBS)

endif  # on USESCOREC


%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $< -o $@
