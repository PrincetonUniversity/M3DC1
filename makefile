ifeq (,$(filter _%,$(notdir $(CURDIR))))

# The following section determines the machine and code version
# and prepares a directory for object files and modules
#--------------------------------------------------------------

# specify whether real or complex
ifeq ($(COM), 1)
OPTS := $(OPTS) -Dvectype=complex -DUSECOMPLEX
BIN_POSTFIX := $(BIN_POSTFIX)-complex
else
OPTS := $(OPTS) -Dvectype=real
BIN_POSTFIX := $(BIN_POSTFIX)-real
endif

# specify whether debug or optimization 
ifeq ($(OPT), 1)
OPTS := $(OPTS) -O
SCORECOPT = -O
BIN_POSTFIX := $(BIN_POSTFIX)-opt
else
OPTS := $(OPTS) -g 
endif

# Define the size of sampling point arrays.
# This sets the upper limit for number of points used
# in numerical integrations
ifdef MAX_PTS
OPTS := $(OPTS) -DMAX_PTS=$(MAX_PTS)
BIN_POSTFIX := $(BIN_POSTFIX)-$(MAX_PTS)
else
OPTS := $(OPTS) -DMAX_PTS=79
endif

export OPTS
export SCORECOPT

include target.mk

else

# The machine-independent parts of the makefile are specified in this section.
# Machine-dependent parts should be put in $HOSTNAME.mk
# ----------------------------------------------------------------------------

VPATH=$(SRCDIR)

BIN = m3dc1

COMMONDIR = $(SRCDIR)/../common/

READGATO_OBJS = polar.o readgato.o
READJSOLVER_OBJS = polar.o read_jsolver_exec.o

OBJS = $(COMMONDIR)subp.o $(COMMONDIR)dbesj0.o $(COMMONDIR)dbesj1.o \
        $(COMMONDIR)fdump.o interpolate.o control.o M3Dmodules.o \
	nintegrate_mod.o metricterms_new.o newvar.o diagnostics.o \
	coils.o gradshafranov.o transport.o hdf5_output.o time_step.o \
	newpar.o fin.o part_fin.o ludef_t.o boundary.o unknown.o \
	restart.o acbauer.o metricterms.o readgeqdsk.o read_dskbal.o \
	read_jsolver.o init_conds.o PETScInterface.o

all : $(BIN)

include $(SRCDIR)/$(ARCH).mk

$(BIN): $(OBJS)
	$(LOADER) $(OBJS) $(LIBS) -o $@

make_polar : make_polar.c
	$(CC) $< -lm -o $@

readgato :  $(READGATO_OBJS)
	$(F90) $(READGATO_OBJS) -L$(NTCCHOME)/lib -lpspline -o $@

endif