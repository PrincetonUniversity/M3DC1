include $(ROOT_DIR)/make.inc

LIBS := $(HDF5_LIBS) $(LIBS)
INCLUDE := $(HDF5_INCLUDE) $(INCLUDE)

OBJS = m3dc1.o m3dc1_file.o

libm3dc1.a : $(OBJS)
	$(ARCH) $@ $(OBJS)
#	$(DYLIB) $(LIBS) $(OBJS) -o $@

clean : 
	rm -f *.o *.so *.a *~

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@