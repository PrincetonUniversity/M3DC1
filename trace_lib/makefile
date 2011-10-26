include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

LIBS := -L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS) 

INCLUDE := -I$(M3DC1_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

OBJS = interpolate.o trace_source.o interpolation_source.o \
	m3dc1_source.o geqdsk_source.o coil_source.o \
	diiid_coils.o trace_integrator.o

trace : main.o libtrace.a
	$(F90) $(LDFLAGS) $< -L. -ltrace $(LIBS) -o $@

libtrace.a : $(OBJS)
	$(ARCH) $@ $(OBJS) 

install : trace
	cp trace $(M3DC1_INSTALL_DIR)

clean : 
	rm -fr *.o *.so *.a *~ _$(M3DC1_ARCH)

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@
