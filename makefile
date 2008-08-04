SHELL=/bin/bash

COMMONDIR = ../common/

INCLUDE = -I$(COMMONDIR) -I$(NTCCHOME)/mod -I$(LIBDIR) \
	-I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH)

LOADER = ifort
F90    = ifort -c
F77    = ifort -c
CC     = icc -c

# define where you want to locate the mesh adapt libraries
# defult is /u/xluo/develop
ifndef SCORECDIR
SCORECDIR = /u/xluo/develop.newCompiler/
#SCORECDIR = /u/nferraro/C1/src/SCOREC/
endif

# define the version of mesh adapt : real or complex version
# ---- Note, this is a temporary solution. Eventually, there should
# ---- be only one version
ifndef SCORECVERS
SCORECVERS = 
endif

# specify whether debug or optimization 
ifeq ($(OPT), 1)
 COMPLEX = -O
 SCORECOPT = -O
else
 COMPLEX = -g 
 SCORECOPT =
endif

# define the complex version
ifeq ($(COM), 1)
COMPLEX := $(COMPLEX) -Dvectype=complex -DUSECOMPLEX
ifeq ($(OPT), 1)
BIN_POSTFIX = -complex-opt
else
BIN_POSTFIX = -complex
endif

endif

# define the real version
ifeq ($(RL), 1)
COMPLEX := $(COMPLEX) -Dvectype=real

ifeq ($(OPT), 1)
BIN_POSTFIX = -real-opt
else
BIN_POSTFIX = -real
endif

endif

BIN = gonewp${BIN_POSTFIX}

#FOPTS = -r8 -implicitnone -fpp $(INCLUDE) ${COMPLEX}
FOPTS = -r8 -implicitnone -save -fpp $(INCLUDE) ${COMPLEX}
#FOPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE) ${COMPLEX}
F90OPTS = ${FOPTS}
F77OPTS = ${FOPTS}
CCOPTS = -c $(INCLUDE)

NEWOBJS = control.o M3Dmodules.o nintegrate_mod.o metricterms_new.o \
	newvar.o diagnostics.o gradshafranov.o \
	$(COMMONDIR)subp.o \
	$(COMMONDIR)dbesj0.o $(COMMONDIR)dbesj1.o \
        $(COMMONDIR)fdump.o hdf5_output.o time_step.o newpar.o \
	fin.o part_fin.o ludef_t.o \
	boundary.o unknown.o restart.o \
	acbauer.o metricterms.o \
	init_conds.o PETScInterface.o

LDRNEW = \
        -L$(SCORECDIR)FMDB$(SCORECVERS)/FMDB/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB$(SCORECVERS)/FMDB/lib/ia64_linux \
	-L$(SCORECDIR)FMDB$(SCORECVERS)/SCORECModel/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB$(SCORECVERS)/SCORECModel/lib/ia64_linux \
	-L$(SCORECDIR)FMDB$(SCORECVERS)/SCORECUtil/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)FMDB$(SCORECVERS)/SCORECUtil/lib/ia64_linux \
	-L$(SCORECDIR)mctk$(SCORECVERS)/Examples/PPPL/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk$(SCORECVERS)/Examples/PPPL/lib/ia64_linux \
	-L$(SCORECDIR)mctk$(SCORECVERS)/Field/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk$(SCORECVERS)/Field/lib/ia64_linux \
	-L$(SCORECDIR)mctk$(SCORECVERS)/Core/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk$(SCORECVERS)/Core/lib/ia64_linux \
	-L$(SCORECDIR)mctk$(SCORECVERS)/Solver/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)mctk$(SCORECVERS)/Solver/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt$(SCORECVERS)/meshAdapt/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt$(SCORECVERS)/meshAdapt/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt$(SCORECVERS)/meshTools/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt$(SCORECVERS)/meshTools/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt$(SCORECVERS)/templateRefine/lib/ia64_linux \
	-Wl,-rpath,$(SCORECDIR)meshAdapt$(SCORECVERS)/templateRefine/lib/ia64_linux \
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

