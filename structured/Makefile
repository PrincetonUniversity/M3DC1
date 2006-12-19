SHELL=/bin/bash

include make.sgi

COMMONDIR = ../common/

OBJS =  $(COMMONDIR)tv80lib.o $(COMMONDIR)subp.o $(COMMONDIR)dbesj0.o \
	$(COMMONDIR)dbesj1.o $(COMMONDIR)fdump.o \
	$(COMMONDIR)superlu_c2f_wrap.o superlupara.o superlu_mod.o \
	$(COMMONDIR)dcreate_dist_matrix.o  $(COMMONDIR)writeHDF5.o \
	supralu_dist_mod.o \
	inout_mod.o newpar.o opdef_mod.o newvar_mod.o init_conds_mod.o \
        mesh_mod.o  gradshafranov.o
	


LDRNEW = \
        -L$(NTCCHOME)/lib -lezcdf \
	-L$(NETCDFHOME)/lib -lnetcdf \
	$(SLULIB) \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	$(FFLIBS) \
	$(LAPACKLIBS) \
	$(HDF5LIB) \
	$(MPILIB) \
	$(XLIB) \
	$(IPRLIB) \
#	-lpthread -lwrap \

all : $(BIN)


$(BIN) : $(OBJS)
	$(LOADER) -o $@ $(OBJS) $(LDRNEW)

clean :
	rm -f $(OBJS) *.mod $(BIN) *~

$(COMMONDIR)tv80lib.o : $(COMMONDIR)tv80lib.f
	$(F77) $(TV80OPTS) $< -o $@

%.o: %.f
	$(F77) $(CFOPTS) -I$(COMMONDIR) $< -o $@

%.o: %.F
	$(F77) $(CFOPTS) -I$(COMMONDIR) $< -o $@

%.o: %.f90
	$(F90) $(CFOPTS) -I$(COMMONDIR) -I$(NETCDFHOME)/include $< -o $@

%.o: %.c
	$(CC) $< -o $@






