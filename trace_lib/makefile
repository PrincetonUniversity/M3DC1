include $(ROOT_DIR)/make.inc

LIBS :=  $(LIBS) $(HDF5_LIBS) \
	-L$(ROOT_DIR)/m3dc1_lib -lm3dc1
INCLUDE := $(INCLUDE) $(HDF5_INCLUDE) \
	-I$(ROOT_DIR)/m3dc1_lib

OBJS = interpolate.o trace_source.o interpolation_source.o \
	m3dc1_source.o geqdsk_source.o coil_source.o \
	diiid_coils.o trace_integrator.o

trace : main.o libtrace.a
	$(CC) $(CFLAGS) main.o -L. -ltrace $(LIBS) -o $@

libtrace.a : $(OBJS) $(HEADERS)
	$(ARCH) $@ $(OBJS) 
#	$(DYLIB) $(LIBS) $(OBJS) -o $@

install : trace
	cp trace $(INSTALL_DIR)/

clean : 
	rm -f trace *.o *.so *.a *~

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@
