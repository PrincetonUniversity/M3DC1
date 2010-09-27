CPP = mpicxx
LOADER = mpif90

LOPTS = -Wl,--warn-unresolved-symbols

INCLUDE = -I/p/tsc/m3dc1/lib/SCORECLib/include/Stix/latest

SCORECDIR = /p/tsc/m3dc1/lib/SCORECLib/lib/Stix/092210

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

AUTOPACK_LIBS = -L$(AUTOPACK_HOME)/lib \
	-Wl,-rpath,$(AUTOPACK_HOME)/lib -lautopack-O

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
	-lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc

SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist

MUMPS_LIBS = -L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord

BLACS_LIBS = -L$(BLACS_HOME)/lib -lmpiblacs

SCALAPACK_LIBS = -L$(SCALAPACK_HOME)/lib -lscalapack

PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

LIBS := $(SCOREC_LIBS) \
	$(AUTOPACK_LIBS) \
	$(PETSC_LIBS) \
	$(SUPERLU_LIBS) \
	$(MUMPS_LIBS) \
	$(SCALAPACK_LIBS) \
	$(BLACS_LIBS) \
	$(PARMETIS_LIBS) \
	-L$(Zoltan_HOME)/lib -lzoltan \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack -lmkl_ipr \
	-L$(ZLIB_HOME) -lz \
        -L/usr/X11R6/lib -lX11 \
	-lifcore

all : structMesh.x main

main : main.o
	$(CPP) $(LOPTS) $(LIBS) $< -o $@

structMesh.x : structMesh.o
	$(CPP) $< -o $@

%.o : %.cc
	$(CPP) -c $< $(INCLUDE) -o $@

clean : 
	rm -fr main structMesh.x *.o