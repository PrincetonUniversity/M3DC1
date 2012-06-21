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

LIBS := -L$(M3DC1_ROOT)/fusion_io/_$(M3DC1_ARCH) -lfusionio \
	-L$(M3DC1_ROOT)/m3dc1_lib/_$(M3DC1_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS) 

INCLUDE := -I$(M3DC1_ROOT)/fusion_io -I$(M3DC1_ROOT)/fusion_io/_$(M3DC1_ARCH) \
	-I$(M3DC1_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

EXAMPLE_BIN = example
EXAMPLE_OBJS = example.o
EXAMPLE_FORTRAN_BIN = example_fortran
EXAMPLE_FORTRAN_OBJS = example_fortran.o
EXAMPLE_PUSH_BIN = example_push
EXAMPLE_PUSH_OBJS = example_push.o 
F90FLAGS := $(F90FLAGS) -DFORTRAN

all : examples

examples: $(EXAMPLE_BIN) $(EXAMPLE_FORTRAN_BIN) $(EXAMPLE_PUSH_BIN)

$(EXAMPLE_BIN) : $(EXAMPLE_OBJS)
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

$(EXAMPLE_FORTRAN_BIN) : $(EXAMPLE_FORTRAN_OBJS)
	$(F90) $(LDFLAGS) $< $(LIBS) -o $@

$(EXAMPLE_PUSH_BIN) : $(EXAMPLE_PUSH_OBJS)
	$(F90) $(LDFLAGS) $< -lfio_push $(LIBS) -o $@

include $(M3DC1_ROOT)/install/common_footer.mk

endif