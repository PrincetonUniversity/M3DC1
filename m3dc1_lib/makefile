include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

LIBS := $(HDF5_LIBS) $(LIBS)
INCLUDE := $(HDF5_INCLUDE) $(INCLUDE)

OBJS = m3dc1_mesh.o m3dc1_field.o m3dc1_file.o

all : libm3dc1.a libm3dc1_fortran.a m3dc1_fortran_test

libm3dc1.a : $(OBJS)
	$(ARCH) $@ $(OBJS)

libm3dc1_fortran.a : m3dc1_fortran.o
	$(ARCH) $@ $<

m3dc1_fortran_test : test.o libm3dc1_fortran.a libm3dc1.a
	$(F90) $< -L. -lm3dc1_fortran -lm3dc1 $(LIBS) -lmpi_cxx -lstdc++ -o $@

install : m3dc1_fortran_test
	cp m3dc1_fortran_test $(M3DC1_INSTALL_DIR)

clean : 
	rm -fr *.o *.so *.a *~ _$(M3DC1_ARCH)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@