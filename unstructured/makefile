SHELL=/bin/bash

COMMONDIR = ../common/

# define where you want to locate the mesh adapt libraries
# defult is /u/xluo/develop
ifndef SCORECDIR
SCORECDIR = /u/xluo/develop.newCompiler.constraint/
endif

INCLUDE = -I$(COMMONDIR) -I$(NTCCHOME)/mod -I$(LIBDIR) \
	-I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) \
	-I$(SCORECDIR)

LOADER = ifort
F90    = ifort -c
F77    = ifort -c
CC     = icc -c

OPTS =

# define the complex version
ifeq ($(COM), 1)
OPTS := $(OPTS) -Dvectype=complex -DUSECOMPLEX
BIN_POSTFIX := $(BIN_POSTFIX)-complex
endif

ifeq ($(RL), 1)
OPTS := $(OPTS) -Dvectype=real
BIN_POSTFIX := $(BIN_POSTFIX)-real
endif

# specify whether debug or optimization 
ifeq ($(OPT), 1)
 OPTS := $(OPTS) -O
 SCORECOPT = -O
 BIN_POSTFIX := $(BIN_POSTFIX)-opt
else
 OPTS := $(OPTS) -g 
 SCORECOPT =
endif


# Define the size of sampling point arrays.
# This sets the upper limit for number of points used
# in numerical integrations
ifdef MAX_PTS
OPTS := $(OPTS) -DMAX_PTS=$(MAX_PTS)
BIN_POSTFIX := $(BIN_POSTFIX)-$(MAX_PTS)
else
OPTS := $(OPTS) -DMAX_PTS=79
endif

BIN = gonewp${BIN_POSTFIX}

#FOPTS = -r8 -implicitnone -fpp $(INCLUDE) ${OPTS}
FOPTS = -r8 -implicitnone -save -fpp $(INCLUDE) $(OPTS) # -check all -check noarg_temp_created
#FOPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE) ${OPTS}
F90OPTS = ${FOPTS}
F77OPTS = ${FOPTS}
CCOPTS = -c $(INCLUDE)

NEWOBJS = $(COMMONDIR)subp.o $(COMMONDIR)dbesj0.o $(COMMONDIR)dbesj1.o \
        $(COMMONDIR)fdump.o \
	control.o M3Dmodules.o nintegrate_mod.o metricterms_new.o \
	newvar.o diagnostics.o gradshafranov.o transport.o \
	hdf5_output.o time_step.o newpar.o \
	fin.o part_fin.o ludef_t.o \
	boundary.o unknown.o restart.o \
	acbauer.o metricterms.o \
	init_conds.o PETScInterface.o

LDRNEW = \
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
	-lPPPL-mpich2$(SCORECOPT) \
	-L$(AUTOPACK_HOME)/lib/ia64-sgi -Wl,-rpath,$(AUTOPACK_HOME)/lib/ia64-sgi -lautopack-O \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(PARMETIS_HOME)/lib -Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis \
	-L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscksp -lpetscmat -lpetscvec -lpetsc \
	-L$(SUPERLU_HOME) -lsuperlu_3.0 \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-L$(CCHOME)/lib -lipr -lstdc++ \
	-L$(MKLHOME)/lib/64 -lguide -lmkl_lapack -lmkl_ipf \
	-L${LIBDIR} -lhdf5_fortran -lhdf5 -lz \
	-Wl,-rpath -Wl,${LIBDIR} \
        -L/usr/X11R6/lib -lX11 -lmpi



${BIN}: $(NEWOBJS)
	$(LOADER) $(NEWOBJS) $(LDRNEW) -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) -fpic $< -o $@

clean:
	rm -f ${BIN}
	rm -f $(NEWOBJS)
	rm -f *.o 
	rm -f *.mod 
	rm -f *~

fullclean:
	rm -f ${BIN}
	rm -r lib*so
	rm -f *.o 
	rm -f *.mod 
	rm -f out.* fort* new.* restartout PI* core*
	rm -f setup* 
	rm -f *~

