include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

LIBS := $(HDF5_LIBS) $(LIBS)
INCLUDE := $(HDF5_INCLUDE) $(INCLUDE)

OBJS = m3dc1_mesh.o m3dc1_field.o m3dc1_file.o

all : libm3dc1.a

libm3dc1.a : $(OBJS)
	$(ARCH) $@ $(OBJS)
#	$(DYLIB) $(LIBS) $(OBJS) -o $@

clean : 
	rm -fr *.o *.so *.a *~ _$(M3DC1_ARCH)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@