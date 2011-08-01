include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

#LIBS := -L$(M3DC1_ROOT)/trace_lib/_$(M3DC1_ARCH) -ltrace \
#	-L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
#	$(HDF5_LIBS) $(LIBS)

LIBS := -L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1_fortran \
	-L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS)


INCLUDE := -I$(M3DC1_ROOT)/trace_lib \
	-I$(M3DC1_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

#OBJS = main.o
OBJS = convert.o

all : m3dc1_convert # libm3dc1_read_lib.a m3dc1_read_test 

m3dc1_convert : $(OBJS)
	$(LD) $(LDFLAGS) $+ $(LIBS) -o $@

#libm3dc1_read_lib.a : m3dc1_read_lib.o
#	$(ARCH) $@ $+
#
#m3dc1_read_test : read_test.o libm3dc1_read_lib.a
#	$(LD) $(LDFLAGS) $< -L. -lm3dc1_read_lib $(LIBS) -o $@

install : m3dc1_convert
	cp m3dc1_convert $(M3DC1_INSTALL_DIR)/

clean : 
	rm -fr *.o *.so *.a *~ _$(M3DC1_ARCH)

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) -fdefault-real-8 $< -c $(F90FLAGS) -o $@
