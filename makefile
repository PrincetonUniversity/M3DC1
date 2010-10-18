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

ifeq ($(USEPETSC), 1)
  OPTS := $(OPTS) -DUSEPETSC -Dmesh_mod=basic_mesh_mod \
	-Dvector_mod=petsc_vector_mod -Dmatrix_mod=petsc_matrix_mod \
	-Dmatrix_type=petsc_matrix -Dvector_type=petsc_vector
  V_OBJ := basic_mesh.o petsc_vector.o petsc_matrix.o
  BIN_POSTFIX := $(BIN_POSTFIX)-petsc
else
  USESCOREC = 1
  OPTS := $(OPTS) -DUSESCOREC -Dmesh_mod=scorec_mesh_mod \
	-Dvector_mod=scorec_vector_mod -Dmatrix_mod=scorec_matrix_mod \
	-Dmatrix_type=scorec_matrix -Dvector_type=scorec_vector
  V_OBJ := scorec_mesh.o scorec_vector.o scorec_matrix.o PETScInterface.o
endif

# specify whether debug or optimization 
ifeq ($(OPT), 1)
  OPTS := $(OPTS) -O
  SCORECOPT = -O
  BIN_POSTFIX := $(BIN_POSTFIX)-opt
else
  SCORECOPT =
endif

ifeq ($(3D), 1)
  OPTS := $(OPTS) -DUSE3D
  BIN_POSTFIX := $(BIN_POSTFIX)-3d
  USE3D = 1
  ifndef MAX_PTS
    MAX_PTS = 125
  endif
else
  USE3D = 0
  ifndef MAX_PTS
    MAX_PTS = 25
  endif
endif

# Define the size of sampling point arrays.
# This sets the upper limit for number of points used
# in numerical integrations
OPTS := $(OPTS) -DMAX_PTS=$(MAX_PTS)
BIN_POSTFIX := $(BIN_POSTFIX)-$(MAX_PTS)

OPTS := $(OPTS) -DPETSC_FORTRAN_PETSCTRUTH_INT # -DxCJ_MATRIX_DUMP

export OPTS
export SCORECOPT
export V_OBJ
export USESCOREC
export USE3D

include target.mk

else

# The machine-independent parts of the makefile are specified in this section.
# Machine-dependent parts should be put in $HOSTNAME.mk
# ----------------------------------------------------------------------------

VPATH=$(SRCDIR)

include $(SRCDIR)/$(ARCH).mk

BIN = m3dc1

READGATO_OBJS = polar.o readgato.o
READJSOLVER_OBJS = polar.o read_jsolver_exec.o

OBJS := $(AUX) subp.o \
	math.o interpolate.o control.o \
	element.o $(V_OBJ) field.o nintegrate_mod.o \
	M3Dmodules.o \
	m3dc1_nint.o vacuum_interface.o boundary.o \
	harned_mikic.o metricterms_new.o biharmonic.o \
	electrostatic_potential.o newvar.o diagnostics.o \
	coils.o coil_sets.o gradshafranov.o transport.o \
	time_step.o hdf5_output.o output.o \
	newpar.o input.o ludef_t.o \
	restart.o readgeqdsk.o read_dskbal.o \
	read_jsolver.o output.o \
	ic_resistive_wall.o \
	init_conds.o

S2V_OBJS = math.o element.o scorec_mesh.o struct2vac.o

all : $(BIN)

$(BIN): $(OBJS)
	$(LOADER) $(OBJS) $(LIBS) -o $@

make_polar : make_polar.c
	$(CC) $< -lm -o $@

readgato :  $(READGATO_OBJS)
	$(F90) $(READGATO_OBJS) -L$(NTCCHOME)/lib -lpspline -o $@

read_jsolver : $(READJSOLVER_OBJS)
	$(F90) $(READJSOLVER_OBJS) -L$(NTCCHOME)/lib -lpspline -o $@

read_jsolver_exec.o : read_jsolver.f90
	$(F90) $(F90OPTS) -DREAD_JSOLVER $< -o $@

struct2vac : $(S2V_OBJS)
	$(LOADER) -Wl,--warn-unresolved-symbols $(S2V_OBJS) $(LIBS) -o $@

endif
