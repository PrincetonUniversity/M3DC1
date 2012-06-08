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
	fusion_io.o \
	c_interface.o \
	fortran_interface.o \
	fusion_io_fortran.o

LIB = libfusionio.a
EXAMPLE_BIN = example
EXAMPLE_OBJS = example.o
EXAMPLE_FORTRAN_BIN = example_fortran
EXAMPLE_FORTRAN_OBJS = example_fortran.o
F90FLAGS := $(F90FLAGS) -DFORTRAN

all : $(LIB) $(EXAMPLE_BIN) $(EXAMPLE_FORTRAN_BIN)

$(LIB) : $(OBJS)
	$(ARCH) $@ $(OBJS)

$(EXAMPLE_BIN) : $(EXAMPLE_OBJS) $(LIB)
	$(LD) $(LDFLAGS) $< -L. -lfusionio $(LIBS) -o $@
$(EXAMPLE_FORTRAN_BIN) : $(EXAMPLE_FORTRAN_OBJS) $(LIB)
	$(F90) $(LDFLAGS) $< -L. -lfusionio $(LIBS) -o $@


.PHONY: install
install : $(LIB) $(EXAMPLE_BIN)
	mkdir -p $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	cp $(LIB) $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	mkdir -p $(M3DC1_INSTALL_DIR)/include/
	cp $(SRCDIR)/fusion_io.h fusion_io.mod $(M3DC1_INSTALL_DIR)/include/

include $(M3DC1_ROOT)/install/common_footer.mk

endif