SHELL=/bin/bash

COMMONDIR = ../common/

INCLUDE = -I$(COMMONDIR) -I$(NTCCHOME)/mod -I$(LIBDIR) \
	  -I$(SUPERLU_DIST_HOME) -I$(HDF5_HOME)/include

LOADER = ifort
F90    = ifort -c
F77    = ifort -c
CC     = icc -c

F90OPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE) -DNEW_VELOCITY
F77OPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE) -DNEW_VELOCITY
#F90OPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE)
#F77OPTS = -r8 -save -Dmpi -ftz -fpp $(INCLUDE)



NEWOBJS1 = M3Dmodules.o nintegrate_mod.o metricterms_n.o metricterms_new.o \
	diagnostics.o gradshafranov.o newvar.o control.o \
	$(COMMONDIR)tv80lib.o $(COMMONDIR)subp.o \
	$(COMMONDIR)dbesj0.o $(COMMONDIR)dbesj1.o \
        $(COMMONDIR)fdump.o hdf5_output.o

NEWOBJS2 = fin.o part_fin.o ludef_t.o \
	  part_fin3.o boundary.o unknown.o restart.o \
	  acbauer.o metricterms.o compare.o \
	  init_conds.o  output.o 

SCORECDIR = /l/mhd/acbauer/develop/
SCORECVERS =-stable5
SCORECOPT = -O

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
	-lFMDB-mpich2$(SCORECOPT) -lSCORECModel-mpich2$(SCORECOPT) -lSCORECUtil-mpich2$(SCORECOPT) -lField-mpich2$(SCORECOPT) -lCore-mpich2$(SCORECOPT) \
	-lmeshAdapt-mpich2$(SCORECOPT) -ltemplateRefine-mpich2$(SCORECOPT) -lmeshTools-mpich2$(SCORECOPT) -lSolver-mpich2$(SCORECOPT) -lPPPL-mpich2$(SCORECOPT) \
	-L$(AUTOPACK_HOME)/lib/ia64-sgi -Wl,-rpath,$(AUTOPACK_HOME)/lib/ia64-sgi -lautopack-O \
	-L$(PARMETIS_HOME)/lib -Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis \
        -L$(NTCCHOME)/lib -lezcdf \
        -L$(NETCDFHOME)/lib -lnetcdf \
        -L$(SUPERLU_DIST_HOME)/lib -lsuperlu \
        -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
        -L$(F90HOME)/lib -lifport -lifcore -lifcoremt -lunwind \
	-Wl,-rpath,$(F90HOME)/lib \
        -L$(MKLHOME)/lib/64 -lguide -lmkl_lapack -lmkl_ipf -lpthread \
        -L${LIBDIR} -lhdf5 -lhdf5_fortran \
	-Wl,--rpath -Wl,${LIBDIR} \
	-L$(CCHOME)/lib -lipr \
	-Wl,-rpath,$(CCHOME)/lib \
        -L/usr/X11R6/lib -lX11 -lmpi -lcprts -lcxa


gonewp: $(NEWOBJS1) newpar.o newpar-lib.o $(NEWOBJS2)
	ifort -shared -o libnewpar.so $(NEWOBJS1) newpar-lib.o $(NEWOBJS2) 
	$(LOADER) $(NEWOBJS1) newpar.o  $(NEWOBJS2) $(LDRNEW) -o $@

newpar-lib.o: newpar.f90
	$(F90) $(F90OPTS) -DIS_LIBRARY $< -o $@

$(COMMONDIR)tv80lib.o: $(COMMONDIR)tv80lib.f
	$(F77) $< -o $@ 

%.o : %.c
	$(CC) $(CCOPTS) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) -fpic $< -o $@

clean:
	rm -f gonewp*
	rm -f $(NEWOBJS1)
	rm -f $(NEWOBJS2)
	rm -f *.o 
	rm -f *.mod 
	rm -f *~
	rm -f *.so

fullclean:
	rm -f gonewp* 
	rm -r lib*so
	rm -f *.o 
	rm -f *.mod 
	rm -f out.* fort* new.* restartout PI* core*
	rm -f setup* 
	rm -f *~







