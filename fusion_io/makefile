#######################################################
# if we're not in the build directory, build it
#######################################################
ifeq (,$(filter _%,$(notdir $(CURDIR))))

include $(M3DC1_ROOT)/install/common_header.mk

#######################################################
# if we are in the build directory, do normal build
#######################################################
else

# set source directory
VPATH=$(SRCDIR)

# include platform-specific options
include $(M3DC1_ROOT)/install/make.inc.$(M3DC1_ARCH)

LIBS := -L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS) 

INCLUDE := -I$(M3DC1_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

LIB = libfusionio.a
DLIB = libfusionio.so
OBJS = 	interpolate.o \
	options.o \
	fusion_io_species.o \
	fusion_io_field.o \
	fio_operations.o \
	compound_field.o \
	m3dc1_source.o \
	m3dc1_field.o \
	geqdsk_source.o \
	geqdsk_field.o \
	gato_source.o \
	gato_field.o \
	fusion_io.o \
	c_interface.o \
	fortran_interface.o \
	fusion_io_fortran.o

PUSH_LIB = libfio_push.a
PUSH_OBJS = push_interface.o

F90FLAGS := $(F90FLAGS) -DFORTRAN

all : $(LIB)

shared : $(DLIB)

push : $(PUSH_LIB)

$(DLIB) : $(OBJS)
	$(LDD) $(LDOPTS) $(OBJS) $(LIBS) -o $@

$(LIB) : $(OBJS)
	$(ARCH) $@ $(OBJS)

$(PUSH_LIB) : $(PUSH_OBJS)
	$(ARCH) $@ $(PUSH_OBJS)

python : fio_setup.py python_interface.cpp
	ARCHFLAGS='$(LIBS)' CXX='$(CXX)' CC='$(CXX)' \
	CFLAGS='$(CFLAGS) $(INCLUDE) $(LIBS)' \
	$(PYTHON) $< build $(PYFLAGS)

.PHONY: install
install :
	mkdir -p $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	-cp $(LIB) $(DLIB) $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	mkdir -p $(M3DC1_INSTALL_DIR)/include/
	-cp $(SRCDIR)/fusion_io.h $(SRCDIR)/fusion_io_defs.h $(M3DC1_INSTALL_DIR)/include/
	mkdir -p $(M3DC1_INSTALL_DIR)/include/_$(M3DC1_ARCH)
	-cp fusion_io.mod FUSION_IO.mod fio_push.mod FIO_PUSH.mod $(M3DC1_INSTALL_DIR)/include/_$(M3DC1_ARCH)
#	CFLAGS='$(CFLAGS)' ARCHFLAGS='$(LIBS)' \
#	$(PYTHON) $(SRCDIR)/fio_setup.py install

include $(M3DC1_ROOT)/install/common_footer.mk

endif