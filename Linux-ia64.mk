LOADER = ifort
F90    = ifort
F77    = ifort
CC     = icc

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = /u/xluo/develop.petsc3/
endif

INCLUDE = -I$(COMMONDIR) -I$(NTCCHOME)/mod -I$(LIBDIR) \
	-I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCORECDIR)

H5_VERSION = 166

FOPTS = -c -r8 -implicitnone -save -fpp -warn unused $(INCLUDE) $(OPTS) -DH5_VERSION=$(H5_VERSION) # -check all -check noarg_temp_created
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
	-L$(SCORECDIR)FMDB/FMDB/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/FMDB/lib/ia64_linux \
	-L$(SCORECDIR)FMDB/SCORECModel/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECModel/lib/ia64_linux \
	-L$(SCORECDIR)FMDB/SCORECUtil/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECUtil/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Examples/PPPL/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Examples/PPPL/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Field/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Field/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Core/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Core/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Solver/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Solver/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/meshAdapt/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshAdapt/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/meshTools/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshTools/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/templateRefine/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/templateRefine/lib/ia64_linux \
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
	$(F90) $(F90OPTS) -fpic $< -o $@
