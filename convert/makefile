include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

LIBS := -L$(M3DC1_ROOT)/trace_lib/_$(M3DC1_ARCH) -ltrace \
	-L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(M3DC1_ROOT)/trace_lib \
	-I$(M3DC1_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

OBJS = main.o

m3dc1_convert : $(OBJS)
	$(LD) $(LDFLAGS) $(CFLAGS) $+ $(LIBS) -o $@

install : m3dc1_convert
	cp m3dc1_convert $(INSTALL_DIR)/

clean : 
	rm -fr *.o *.so *.a *~ _$(M3DC1_ARCH)

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@
