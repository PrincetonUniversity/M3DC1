SHELL=/bin/bash

CFOPTS = -r8 -save -Dmpi -ftz -fpp  -I${LIBDIR}
TVOPTS =

NEWDIR = ./
NEWINC = ./

COMMONDIR = ../common/

NTCCMOD = $(NTCCHOME)/mod
NEWFLAGS =  $(CFOPTS)  -I$(NEWINC) -I$(COMMONDIR) -I$(NTCCMOD)
TV80FLAGS = $(TVOPTS)  -I$(NEWINC) -I$(NTCCMOD)

F90 = ifort -c
COMPNEW = $(F90) -c $(NEWFLAGS) $<
COMPNEWLIB = $(F90) -c $(NEWFLAGS) -o newpar-lib.o -DIS_LIBRARY $(NEWDIR)newpar.f90

COMPTV80 = ifort -c $(TV80FLAGS)  $(COMMONDIR)$(@F:.o=.f)
COMPBESJ = ifort -c $(NEWFLAGS)  $(COMMONDIR)$(@F:.o=.f)
LOADER = ifort
CCOMPILE = icc -c

NEWOBJS1 = M3Dmodules.o nintegrate_mod.o metricterms_n.o newvar.o \
	$(COMMONDIR)tv80lib.o $(COMMONDIR)subp.o \
	$(COMMONDIR)dbesj0.o $(COMMONDIR)dbesj1.o \
        $(COMMONDIR)fdump.o hdf5_output.o \
#        $(COMMONDIR)writeHDF5.o 

NEWOBJS2 = fin.o part_fin.o ludef_t.o \
	  part_fin3.o boundary.o unknown.o restart.o \
	  acbauer.o sort.o metricterms.o errorcalc.o compare.o \
	  gradshafranov.o init_conds.o  output.o 

SCORECDIR = /l/mhd/acbauer/develop/
SCORECVERS =-stable4

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
	-lFMDB-mpich2-O -lSCORECModel-mpich2-O -lSCORECUtil-mpich2-O -lField-mpich2-O -lCore-mpich2-O \
	-lmeshAdapt-mpich2-O -ltemplateRefine-mpich2-O -lmeshTools-mpich2-O -lSolver-mpich2-O -lPPPL-mpich2-O \
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
	$(LOADER) -o $@ $(NEWOBJS1) newpar.o  $(NEWOBJS2) $(LDRNEW)
	rm -rf C1restartout check.txt check1.txt check0.txt

newpar-lib.o: $(NEWDIR)newpar.f90
	$(COMPNEWLIB)

newpar.o: $(NEWDIR)newpar.f90
	$(COMPNEW)
	$(COMPNEWLIB)

acbauer.o: $(NEWDIR)acbauer.F
	$(COMPNEW)

compare.o: $(NEWDIR)compare.F
	$(COMPNEW)

errorcalc.o: $(NEWDIR)errorcalc.F
	$(COMPNEW)

metricterms.o: $(NEWDIR)metricterms.F
	$(COMPNEW)

fin.o: $(NEWDIR)fin.F
	$(COMPNEW)

ludef.o: $(NEWDIR)ludef.F
	$(COMPNEW)

sort.o: $(NEWDIR)sort.F
	$(COMPNEW)

restart.o: $(NEWDIR)restart.F
	$(COMPNEW)

part_fin.o: $(NEWDIR)part_fin.F
	$(COMPNEW)

part_fin2.o: $(NEWDIR)part_fin2.F
	$(COMPNEW)

part_fin3.o: $(NEWDIR)part_fin3.F
	$(COMPNEW)

unknown.o: $(NEWDIR)unknown.F
	$(COMPNEW)

sparse_params.o: $(NEWDIR)sparse_params.F
	$(COMPNEW)

sparse_matrix.o: $(NEWDIR)sparse_matrix.F
	$(COMPNEW)

$(COMMONDIR)tv80lib.o: $(COMMONDIR)tv80lib.f
	$(COMPTV80) -o $@ 

$(COMMONDIR)dbesj0.o: $(COMMONDIR)dbesj0.f
	$(COMPBESJ) -o $@

$(COMMONDIR)dbesj1.o: $(COMMONDIR)dbesj1.f
	$(COMPBESJ) -o $@

$(COMMONDIR)fdump.o: $(COMMONDIR)fdump.f
	$(COMPBESJ) -o $@

$(COMMONDIR)subp.o: $(COMMONDIR)subp.f90
	ifort -c $(NEWFLAGS) $< -o $@

basic_mod.o: $(NEWDIR)basic_mod.f90
	ifort -c $(NEWFLAGS) $< -o $@

$(COMMONDIR)superlu_mod.o: $(COMMONDIR)superlu_mod.f90
	ifort -c $(NEWFLAGS) $< -o $@

supralu_dist_mod.o: supralu_dist_mod.f90
	ifort -c $(NEWFLAGS) $< -o $@

$(COMMONDIR)superlu_c2f_wrap.o: $(COMMONDIR)superlu_c2f_wrap.c
	$(CCOMPILE)  -I${SUPERLU_DIST_HOME} \
	$< -c -o $@

$(COMMONDIR)dcreate_dist_matrix.o: $(COMMONDIR)dcreate_dist_matrix.c
	$(CCOMPILE)  -I${SUPERLU_DIST_HOME} \
	$< -c -o $@

c_fortran_dgssv.o: $(NEWDIR)c_fortran_dgssv.c
	$(CCOMPILE) $< -c -o $@

$(COMMONDIR)writeHDF5.o: $(COMMONDIR)writeHDF5.c
	$(CCOMPILE)  -I$(HDF5_HOME)/include $< -o $@


%.o: %.c
	-rm -f $(*F)_err
	-date > $(*F)_err
	$(CCOMPILE)  -o $@ $<  >>& $(*F)_err

%.o: %.F
	$(F90) $(CFOPTS) $< -o $@

%.o: %.f90
	$(F90) $(CFOPTS) $< -o $@

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







