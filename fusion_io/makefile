#######################################################
# if we're not in the build directory, build it
#######################################################
ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

OBJDIR := _$(M3DC1_ARCH)

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/makefile \
	SRCDIR=$(CURDIR) $(MAKECMDGOALS) 

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: clean
clean : 
	rm -fr *~ $(OBJDIR)


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

OBJS = options.o \
	compound_field.o \
	m3dc1_field.o \
	m3dc1_source.o 

LIB = libfusionio.a
EXAMPLE_BIN = example
EXAMPLE_OBJS = example.o

all : $(LIB) $(EXAMPLE_BIN)

$(LIB) : $(OBJS)
	$(ARCH) $@ $(OBJS)

$(EXAMPLE_BIN) : $(EXAMPLE_OBJS) $(LIB)
	$(LD) $(LDFLAGS) $< -L. -lfusionio $(LIBS) -o $@

.PHONY: install
install : $(LIB) $(EXAMPLE_BIN)
	mkdir -p $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	cp $(LIB) $(M3DC1_INSTALL_DIR)/lib/_$(M3DC1_ARCH)
	mkdir -p $(M3DC1_INSTALL_DIR)/include/
	cp fusion_io.h fusion_io.mod $(M3DC1_INSTALL_DIR)/include/

%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CXX) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@

endif