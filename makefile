ifeq (,$(filter _%,$(notdir $(CURDIR))))

# The following section determines the machine and code version
# and prepares a directory for object files and modules
#--------------------------------------------------------------

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
  ifeq ($(TRILINOS), 1)
    V_OBJ := scorec_mesh.o scorec_vector.o scorec_matrix.o
  else
    V_OBJ := scorec_mesh.o scorec_vector.o scorec_matrix.o PETScInterface.o
  endif
  ifeq ($(RW), 1)
    USERW = 1
    OPTS := $(OPTS) -DUSERW
    BIN_POSTFIX := $(BIN_POSTFIX)-rw
  else
    OPTS := $(OPTS) -Dglobalinsertval=insertval -Dglobalentdofs=entdofs
  endif
endif

# determine whether 2d, 3d, or 2d-complex
ifeq ($(3D), 1)
  OPTS := $(OPTS) -DUSE3D -Dvectype=real
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

  # specify whether real or complex
  ifeq ($(COM), 1)
    OPTS := $(OPTS) -Dvectype=complex -DUSECOMPLEX
    BIN_POSTFIX := $(BIN_POSTFIX)-complex
    USECOMPLEX = 1
  else
    OPTS := $(OPTS) -Dvectype=real
    USECOMPLEX = 0
  endif
endif

# specify whether debug or optimization 
ifeq ($(OPT), 1)
  SCORECOPT = -O
  BIN_POSTFIX := $(BIN_POSTFIX)-opt
else
  SCORECOPT =
endif

ifeq ($(TAU), 1)
  BIN_POSTFIX := $(BIN_POSTFIX)-tau
endif

# Define the size of sampling point arrays.
# This sets the upper limit for number of points used
# in numerical integrations
OPTS := $(OPTS) -DMAX_PTS=$(MAX_PTS)
BIN_POSTFIX := $(BIN_POSTFIX)-$(MAX_PTS)

ifeq ($(TRILINOS),1)
  OPTS := $(OPTS) -DM3DC1_TRILINOS
  BIN_POSTFIX := $(BIN_POSTFIX)-trilinos
endif

ifeq ($(USEADIOS), 1)
  BIN_POSTFIX := $(BIN_POSTFIX)-adios
endif

ifeq ($(OMP), 1)
  BIN_POSTFIX := $(BIN_POSTFIX)-omp
endif
OPTS := $(OPTS) -DPETSC_FORTRAN_PETSCTRUTH_INT #-DCJ_MATRIX_DUMP


# add date stamp
OPTS := $(OPTS) -DNOUSE -DDATE_BUILT="'$(shell date)'" -DRELEASE_VERSION="'$(shell cat release_version)'" #-DBUILD_INFO="'$(shell svn info)'"

export OPT
export OPTS
export SCORECOPT
export V_OBJ
export USESCOREC
export USECOMPLEX
export USE3D
export USERW
export TAU
export HPCTK

include target.mk

else

# The machine-independent parts of the makefile are specified in this section.
# Machine-dependent parts should be put in $HOSTNAME.mk
# ----------------------------------------------------------------------------

VPATH=$(SRCDIR)

include $(SRCDIR)/$(ARCH).mk
BIN = m3dc1

ifeq ($(3D), 1)
  BIN := $(BIN)_3d
else
  BIN := $(BIN)_2d
endif
ifeq ($(TRILINOS),1)
  BIN := $(BIN)_trilinos
endif
ifeq ($(COM), 1)
  BIN := $(BIN)_complex
endif
ifeq ($(USEADIOS), 1)
  BIN := $(BIN)_adios
endif

READGATO_OBJS = polar.o readgato.o
READJSOLVER_OBJS = polar.o read_jsolver_exec.o

OBJS := $(AUX) fftw_fortran.o read_namelist.o gsl_wrapper.o \
	subp.o random.o spline.o \
	math.o read_ascii.o interpolate.o control.o \
	iterdb.o read_gyro.o read_neo.o radiation.o \
	element.o $(V_OBJ) field.o nintegrate_mod.o \
	M3Dmodules.o resistive_wall.o \
	m3dc1_nint.o boundary.o gyroviscosity.o bootstrap.o \
	metricterms_new.o two_fluid.o harned_mikic.o biharmonic.o \
	electric_field.o pellet.o \
	temperature_plots.o \
	electrostatic_potential.o newvar.o diagnostics.o \
	read_schaffer_field.o  neutral_beam.o \
	coils.o coil_sets.o model.o \
	fit_magnetics.o \
	init_common.o \
	gradshafranov.o rmp.o \
	readgeqdsk.o read_jsolver.o read_dskbal.o \
	init_rwm.o init_solovev.o init_circle.o init_basicj.o \
	init_tilt.o init_taylor.o init_force_free.o init_gem.o init_wave.o \
	init_gmode.o init_strauss.o init_mri.o init_rotating_cylinder.o \
        init_eqdsk.o init_dskbal.o init_jsolver.o \
	init_3dwave.o init_3ddiffusion.o init_frs.o init_ftz.o init_eigen.o \
	init_intkink.o init_lz.o init_kstar.o init_basicq.o \
	transport.o \
	auxiliary_fields.o  \
	time_step_split.o time_step_unsplit.o \
	time_step.o hdf5_output.o output.o \
	particle.o \
	error_estimate.o adapt.o newpar.o input.o ludef_t.o \
	restart.o \
	init_conds.o \
	get_pc_skip_count.o

S2V_OBJS = subp.o math.o element.o scorec_mesh.o struct2vac.o
A2CC_OBJS = readaeqdsk.o a2cc.o

all : $(BIN) m3dc1

$(BIN): $(OBJS)
	$(LOADER) $(LDOPTS) $(OBJS) $(LIBS) -o $@

m3dc1 : $(BIN)
	ln -s $< $@

readgato :  $(READGATO_OBJS)
	$(F90) $(READGATO_OBJS) -L$(PSPLINE_DIR)/lib -lpspline -o $@

read_jsolver : $(READJSOLVER_OBJS)
	$(F90) $(READJSOLVER_OBJS) -L$(PSPLINE_DIR)/lib -lpspline -o $@

read_jsolver_exec.o : read_jsolver.f90
	$(F90) $(F90OPTS) -DREAD_JSOLVER $< -o $@

struct2vac : $(S2V_OBJS)
	$(LOADER) -Wl,--warn-unresolved-symbols $(S2V_OBJS) $(LIBS) -o $@

a2cc : $(A2CC_OBJS)
	$(LOADER) $(A2CC_OBJS) -o $@

endif
