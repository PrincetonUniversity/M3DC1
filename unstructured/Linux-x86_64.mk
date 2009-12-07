LOADER = mpif90
F90    = mpif90
F77    = mpif90
CC     = mpicc

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = /u/xluo/develop.petsc3.stix.intel/
endif

INCLUDE = -I$(SCORECDIR) \
	-I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib

H5_VERSION = 169

FOPTS = -c -r8 -implicitnone -fpp -warn unused $(INCLUDE) $(OPTS) -DH5_VERSION=$(H5_VERSION) # -g -check all -check noarg_temp_created
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)
CCOPTS = -c $(INCLUDE)

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib -lautopack-O

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc \

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

SCOREC_ARCH=x86_64_linux-icc

SCOREC_LIBS = \
	-L$(SCORECDIR)FMDB/FMDB/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)FMDB/FMDB/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECModel/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECUtil/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)mctk/Examples/PPPL/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)mctk/Field/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)mctk/Field/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)mctk/Core/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)mctk/Core/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)mctk/Solver/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)mctk/Solver/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshAdapt/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshTools/lib/$(SCOREC_ARCH) \
	-L$(SCORECDIR)meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/templateRefine/lib/$(SCOREC_ARCH) \
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

LIBS = 	$(SCOREC_LIBS) \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack -lmkl_ipr \
	$(AUTOPACK_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	$(PARMETIS_LIBS) \
	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
        -L/usr/X11R6/lib -lX11

%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) -fpic $< -o $@
