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

LIBS := $(HDF5_LIBS) $(LIBS)
INCLUDE := $(HDF5_INCLUDE) $(INCLUDE)

OBJS = m3dc1_mesh.o m3dc1_field.o m3dc1_file.o

all : libm3dc1.a libm3dc1_fortran.a libm3dc1.so

libm3dc1.so : $(OBJS)
	$(LDD) $(OBJS) -o $@

libm3dc1.a : $(OBJS)
	$(ARCH) $@ $(OBJS)

libm3dc1_fortran.a : m3dc1_fortran.o
	$(ARCH) $@ $<

m3dc1_fortran_test : test.o libm3dc1_fortran.a libm3dc1.a
	$(F90) $< -L. -lm3dc1_fortran -lm3dc1 $(LDFLAGS) $(LIBS) -o $@

.PHONY: install
install : libm3dc1.a libm3dc1.so libm3dc1_fortran.a
	mkdir -p $(M3DC1_INSTALL_DIR)/bin/_$(M3DC1_ARCH)
	mkdir -p $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	cp libm3dc1.a libm3dc1.so libm3dc1_fortran.a \
	 $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)

include $(M3DC1_ROOT)/install/common_footer.mk

endif
