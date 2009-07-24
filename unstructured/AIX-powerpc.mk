LOADER = mpCC_r -blpdata
PP     = xlc_r -E -C
F90    = mpxlf90_r -c
F77    = mpxlf_r -c
CC     = mpCC_r -c

# define where you want to locate the mesh adapt libraries
ifndef SCORECDIR
  SCORECDIR = /project/projectdirs/mp288/scorec/develop.newCompiler.constraint/
endif


INCLUDE = -I$(COMMONDIR) \
	-I$(INCLUDE_PATH) $(HDF5_INCLUDE) \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) \
	-I$(SCORECDIR)

PPOPTS = -DPETSC_AVOID_MPIF_H -Dmpi -DNEW_VELOCITY \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/rs6000_64_O \
	$(OPTS)
F90OPTS = -qautodbl=dbl4 -qsave $(QEXTNAME) -qstrict -O2 $(INCLUDE)
F77OPTS = -qautodbl=dbl4 -qsave $(QEXTNAME) -O2 $(INCLUDE)
CCOPTS = -O2 $(INCLUDE)
DBESJDBESJOPTS = -O2 -qstrict -qautodbl=dbl4 $(INCLUDE) $(QEXTNAME)
TV80OPTS=  -qsave -O2 -qstrict $(INCLUDE) $(QEXTNAME)

LIBS = \
        -L$(SCORECDIR)FMDB/FMDB/lib/ia64_linux \
	-L$(SCORECDIR)FMDB/SCORECModel/lib/ia64_linux \
	-L$(SCORECDIR)FMDB/SCORECUtil/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Examples/PPPL/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Field/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Core/lib/ia64_linux \
	-L$(SCORECDIR)mctk/Solver/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/meshAdapt/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/meshTools/lib/ia64_linux \
	-L$(SCORECDIR)meshAdapt/templateRefine/lib/ia64_linux \
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
	-L$(AUTOPACK_HOME)/lib/IBM -lautopack-O \
	$(ZOLTAN) \
	$(METIS) \
	-L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscksp -lpetscmat -lpetscvec -lpetsc \
	$(SUPERLU) $(SUPERLU_DIST) \
	$(LAPACK) \
	$(NCAR) \
	$(HDF5) \
	$(NAG) \
	-lxlfpmt4 -L/usr/vacpp/lib \
	-lC -lpthreads -lm -lc -lxlopt -lxlf -lxlomp_ser -lxlf90

%.o : %.c
	$(CC)  $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(PP) $(PPOPTS) $< > $<.i
	$(F77) $(F77OPTS) -qsuffix=f=i $<.i -o $@
	rm $<.i

%.o: %.f90
	$(PP) $(PPOPTS) $< > $<.i
	$(F90) $(F90OPTS) -qsuffix=f=i $<.i -o $@
	rm $<.i
