LOADER = ifort
F90    = ifort
F77    = ifort
CC     = icc

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib
endif

INCLUDE = -I$(NTCCHOME)/mod -I$(LIBDIR) \
	-I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCORECDIR)/include/Viz/latest

H5_VERSION = 166

FOPTS = -c -r8 -implicitnone -fpp -warn all $(INCLUDE) $(OPTS) -DH5_VERSION=$(H5_VERSION) # -g -check all -check noarg_temp_created
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

SCOREC_LIBS = \
	-L$(SCORECDIR)/lib/Viz/latest \
	-Wl,-rpath,$(SCORECDIR)/lib/Viz/latest \
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


LIBS = $(SCOREC_LIBS) \
	$(AUTOPACK_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	$(PARMETIS_LIBS) \
	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-L$(CCHOME)/lib -lipr -lstdc++ \
	-L$(MKLHOME)/lib/64 -lguide -lmkl_lapack -lmkl_ipf \
	-L$(LIBDIR) -lhdf5_fortran -lhdf5 -lz \
	-Wl,-rpath -Wl,$(LIBDIR) \
        -L/usr/X11R6/lib -lX11 -lmpi


%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $< -o $@
