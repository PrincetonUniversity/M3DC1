LOADER = openmpif90
F90    = openmpif90
F77    = openmpif77
CC     = openmpicc

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = 
endif

H5_VERSION = 184

MPIHOME = /opt/local/include/openmpi
HDF5_HOME = /usr/local/hdf5-1.8.5/hdf5

INCLUDE = -I$(SCORECDIR)/mctk/Examples/PPPL/PPPL \
	-I$(MPIHOME) -I/opt/local/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include

FOPTS = -c -fdefault-real-8 -Wall $(INCLUDE) $(OPTS) \
	-DH5_VERSION=$(H5_VERSION) -DPETSC_31 # -g -fbounds-check
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)
CCOPTS = -c $(INCLUDE)

NAG_LIBS = # -L$(NAG_ROOT)/lib -lnag

ifeq ($(USECOMPLEX), 1)
  PETSC_ARCH = arch-osx-10.6-complex
endif

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
	-lmumps_common -ldmumps -lzmumps -lpord \
	-lblas -lblacs -llapack -lscalapack

SUPERLU_LIBS = -lsuperlu_dist_2.3

PARMETIS_LIBS = -lmetis -lparmetis

HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-lhdf5hl_fortran -lhdf5_hl -lz

ifeq ($(USESCOREC), 1)
AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib -lautopack-O

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

ZOLTAN_LIBS = -L$(Zoltan_HOME)/lib -lzoltan

endif


LIBS = 	-L/opt/local/lib \
	$(SCOREC_LIBS) \
	$(AUTOPACK_LIBS) \
	$(ZOLTAN_LIBS) \
	$(PARMETIS_LIBS) \
	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(HDF5_LIBS) \
	-L/usr/X11R6/lib -lX11

AUX = i1mach.o r1mach.o d1mach.o

%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	cp $< $(*F).tmp.F90
	$(F90) $(F90OPTS) -fpic $(*F).tmp.F90 -o $@
	rm $(*F).tmp.F90
