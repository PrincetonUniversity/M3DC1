LOADER = mpiCC
F90    = mpif90
F77    = mpif90
CC     = mpicc

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = /u/xluo/develop.petsc3.stix/
endif

#INCLUDE = -I$(COMMONDIR) -I$(SCORECDIR) -I$(INCLUDE_PATH)

INCLUDE = -I$(COMMONDIR) -I$(SCORECDIR) \
	-I/usr/pppl/pathscale/3.2-pkgs/vsmpich-1.2.7/include \
	-I/usr/pppl/pathscale/3.2-pkgs/vsmpich-pkgs/PETSc-3.0.0-p6/path-vsmp/include \
	-I/usr/pppl/pathscale/3.2-pkgs/hdf5-1.8.1-serial/include


FOPTS = -c -r8 -cpp $(INCLUDE) $(OPTS) # -check all -check noarg_temp_created
F90OPTS = $(FOPTS)
F77OPTS = $(FOPTS)
CCOPTS = -c $(INCLUDE)

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib/x86_64-linux \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib/x86_64-linux -lautopack-ethernet-O

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc \

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

SCOREC_LIBS = \
	-L$(SCORECDIR)FMDB/FMDB/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/FMDB/lib/x86_64_linux \
	-L$(SCORECDIR)FMDB/SCORECModel/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECModel/lib/x86_64_linux \
	-L$(SCORECDIR)FMDB/SCORECUtil/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB/SCORECUtil/lib/x86_64_linux \
	-L$(SCORECDIR)mctk/Examples/PPPL/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Examples/PPPL/lib/x86_64_linux \
	-L$(SCORECDIR)mctk/Field/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Field/lib/x86_64_linux \
	-L$(SCORECDIR)mctk/Core/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Core/lib/x86_64_linux \
	-L$(SCORECDIR)mctk/Solver/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk/Solver/lib/x86_64_linux \
	-L$(SCORECDIR)meshAdapt/meshAdapt/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshAdapt/lib/x86_64_linux \
	-L$(SCORECDIR)meshAdapt/meshTools/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/meshTools/lib/x86_64_linux \
	-L$(SCORECDIR)meshAdapt/templateRefine/lib/x86_64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt/templateRefine/lib/x86_64_linux \
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
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz \
	-Wl,-rpath -Wl,$(HDF5_HOME) \
        -L/usr/X11R6/lib -lX11 \
	-L$(ACML_HOME)/pathscale64/lib -lacml


%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) -fpic $< -o $@
