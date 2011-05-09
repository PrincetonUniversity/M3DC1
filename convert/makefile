include $(ROOT_DIR)/make.inc

LIBS := -L$(ROOT_DIR)/trace_lib -ltrace \
	-L$(ROOT_DIR)/m3dc1_lib -lm3dc1 \
	$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(ROOT_DIR)/trace_lib \
	-I$(ROOT_DIR)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

OBJS = main.o

m3dc1_convert : $(OBJS)
	$(LD) $(LDFLAGS) $(CFLAGS) $+ $(LIBS) -o $@

install : m3dc1_convert
	cp m3dc1_convert $(INSTALL_DIR)/

clean : 
	rm -f m3dc1_convert *.o *.so *.a *~

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@
