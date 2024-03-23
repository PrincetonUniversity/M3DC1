ifeq (,$(filter _%,$(notdir $(CURDIR))))

# The following section determines the machine and code version
# and prepares a directory for object files and modules
#--------------------------------------------------------------

ifeq ($(USEPETSC), 1)
  OPTS := -DUSEPETSC -Dmesh_mod=basic_mesh_mod \
	-Dvector_mod=petsc_vector_mod -Dmatrix_mod=petsc_matrix_mod \
	-Dmatrix_type=petsc_matrix -Dvector_type=petsc_vector
  V_OBJ := basic_mesh.o petsc_vector.o petsc_matrix.o
  BIN_POSTFIX := $(BIN_POSTFIX)-petsc
else
  USESCOREC = 1
  OPTS := -DUSESCOREC -Dmesh_mod=scorec_mesh_mod \
	-Dvector_mod=scorec_vector_mod -Dmatrix_mod=scorec_matrix_mod \
	-Dmatrix_type=scorec_matrix -Dvector_type=scorec_vector
  ifeq ($(TRILINOS), 1)
    V_OBJ := scorec_mesh.o scorec_vector.o scorec_matrix.o
  else
    V_OBJ := scorec_mesh.o scorec_vector.o scorec_matrix.o #PETScInterface.o
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

ifeq ($(ST), 1)
  OPTS := $(OPTS) -DUSEST 
  BIN_POSTFIX := $(BIN_POSTFIX)-st
  USEST = 1
else 
  USEST = 0 
endif

ifeq ($(REORDERED), 1)
  OPTS := $(OPTS) -DREORDERED
endif

ifeq ($(SCORECVER), reordered)
  OPTS := $(OPTS) -DREORDERED
endif

ifdef SCORECVER
  BIN_POSTFIX := $(BIN_POSTFIX)-$(SCORECVER)
endif

ifeq ($(PAR), 1)
  BIN_POSTFIX := $(BIN_POSTFIX)-pic
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

ifeq ($(OMP), 1)
  BIN_POSTFIX := $(BIN_POSTFIX)-omp
endif

# ifeq ($(ACC), 1)
#   BIN_POSTFIX := $(BIN_POSTFIX)-acc
# endif

OOPTS := $(OPTS) -DPETSC_FORTRAN_PETSCTRUTH_INT #-DCJ_MATRIX_DUMP

# add date stamp
OPTS := $(OPTS) -DDATE_BUILT='"$(shell date)"' -DRELEASE_VERSION='"$(shell cat release_version)"' #-DBUILD_INFO="'$(shell svn info)'"

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

USEADAS = 0

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
ifeq ($(OMP), 1)
  BIN := $(BIN)_omp
endif
ifeq ($(ST), 1)
  BIN := $(BIN)_st
endif

READGATO_OBJS = polar.o readgato.o
READJSOLVER_OBJS = polar.o read_jsolver_exec.o

ifeq ($(ADAS), 1)
  ifeq  ($(USEADAS), 0)
    USEADAS = 1
    LIBS += -L$(SRCDIR)/_$(ARCH) -ladaslib
  endif
endif

ifeq ($(USEADAS), 1)
  OPTS := $(OPTS) -DUSEADAS
  ADAS_OBJ = adas_m3dc1.o
endif

OBJS := fftw_fortran.o read_namelist.o gsl_wrapper.o \
	random.o spline.o region.o \
	math.o read_ascii.o interpolate.o control.o \
	iterdb.o read_gyro.o read_neo.o radiation.o \
	element.o read_vmec.o physical_mesh.o $(V_OBJ) field.o nintegrate_mod.o \
	M3Dmodules.o resistive_wall.o $(ADAS_OBJ) kprad.o \
	 m3dc1_nint.o boundary.o gyroviscosity.o bootstrap.o \
	metricterms_new.o two_fluid.o harned_mikic.o biharmonic.o \
	electric_field.o pellet.o parallel_heat_flux.o \
	temperature_plots.o electrostatic_potential.o \
	newvar.o kprad_m3dc1.o diagnostics.o runaway.o \
	read_schaffer_field.o  neutral_beam.o \
	coils.o coil_sets.o model.o hypervisc.o \
	fit_magnetics.o wall.o \
	geometry.o\
	init_common.o init_vmec.o\
	gradshafranov.o rmp.o  \
	readgeqdsk.o read_jsolver.o read_dskbal.o \
	init_rwm.o init_solovev.o init_circle.o init_basicj.o \
	init_tilt.o init_taylor.o init_force_free.o init_gem.o init_wave.o \
	init_gmode.o init_strauss.o init_mri.o init_rotating_cylinder.o \
        init_eqdsk.o init_dskbal.o init_jsolver.o \
	init_3dwave.o init_3ddiffusion.o init_frs.o init_ftz.o init_eigen.o \
	init_intkink.o init_lz.o init_kstar.o init_basicq.o init_cyl.o \
	transport.o \
	auxiliary_fields.o  \
	time_step_split.o time_step_unsplit.o \
        hdf5_output.o particle.o   \
	time_step.o output.o restart_hdf5.o \
	error_estimate.o scorec_adapt.o adapt.o newpar.o input.o ludef_t.o \
	init_conds.o \
	get_pc_skip_count.o \
	parse_solver_options.o

ifneq (,$(findstring gpu,$(ARCH)))
	OBJS := $(OBJS:metricterms_new.o=metricterms_new_gpu.o)
	OBJS := $(OBJS:ludef_t.o=ludef_t_gpu.o)
endif

S2V_OBJS = subp.o math.o element.o scorec_mesh.o struct2vac.o
A2CC_OBJS = readaeqdsk.o a2cc.o

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

a2cc : $(A2CC_OBJS)
	$(LOADER) $(LDOPTS) $(A2CC_OBJS) -o $@

mackenbach_profiles : mackenbach_profiles.f90
	$(F90) $(R8OPTS) $< -o $@

kprad : kprad.f90 kprad_test.f90
	$(F90) $(F90OPTS) $< -o $@

test_region : region.o test_region.o
	$(F90) $^ -o $@

OA_DATA = https://open.adas.ac.uk/download/adf11
OA_OBJS = xxrptn.o xfelem.o i4unit.o i4fctn.o xxword.o xxcase.o \
	xxslen.o xxdata_11.o
.PHONY : download_open_adas
download_open_adas :
	rm -rf adf11
	mkdir -p adf11/scd85 adf11/plt89 adf11/scd89 adf11/plt96 adf11/scd96
	wget https://open.adas.ac.uk/code/xxdata_11.tar.gz
	tar -zxf xxdata_11.tar.gz -C $(SRCDIR) --strip-components=1 --wildcards '*.for' --transform 's/\.for/.f/'
	rm -rf xxdata_11.tar.gz
	cd adf11/scd85 && \
		wget $(OA_DATA)/scd85/scd85_ar.dat
	cd adf11/scd89 && \
		wget $(OA_DATA)/scd89/scd89_b.dat
	cd adf11/plt89 && \
		wget $(OA_DATA)/plt89/plt89_b.dat && \
		wget $(OA_DATA)/plt89/plt89_ar.dat
	cd adf11/scd96 && \
		wget $(OA_DATA)/scd96/scd96_he.dat && \
		wget $(OA_DATA)/scd96/scd96_be.dat && \
		wget $(OA_DATA)/scd96/scd96_c.dat && \
		wget $(OA_DATA)/scd96/scd96_ne.dat
	cd adf11/plt96 && \
		wget $(OA_DATA)/plt96/plt96_he.dat && \
		wget $(OA_DATA)/plt96/plt96_be.dat && \
		wget $(OA_DATA)/plt96/plt96_c.dat && \
		wget $(OA_DATA)/plt96/plt96_ne.dat

.PHONY : open_adas
open_adas : download_open_adas $(OA_OBJS)
	ar rcs libadaslib.a $(OA_OBJS)

.PHONY : bin
bin : a2cc mackenbach_profiles
	rm -fr bin
	mkdir -p bin
	if [ -d $(SCOREC_UTIL_DIR)/create_smb ]; then \
		ln -s $(SCOREC_UTIL_DIR)/create_smb/create_smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/create_smb/seed0.smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/split_smb/split_smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/m3dc1_meshgen/m3dc1_meshgen bin ; \
		ln -s $(SCOREC_UTIL_DIR)/m3dc1_meshgen/convert_sim_sms bin ; \
	else \
		ln -s $(SCOREC_UTIL_DIR)/create_smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/seed0.smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/split_smb bin ; \
		ln -s $(SCOREC_UTIL_DIR)/collapse bin ; \
		ln -s $(SCOREC_UTIL_DIR)/m3dc1_meshgen bin ; \
		ln -s $(SCOREC_UTIL_DIR)/m3dc1_mfmgen bin ; \
		ln -s $(SCOREC_UTIL_DIR)/convert_sim_sms bin ; \
	fi
	if [ -f $(SCOREC_UTIL_DIR)/split_smb/make_model ]; then \
		ln -s $(SCOREC_UTIL_DIR)/split_smb/make_model bin ; \
	elif [ -f $(SCOREC_UTIL_DIR)/create_smb/make_model ]; then \
		ln -s $(SCOREC_UTIL_DIR)/create_smb/make_model bin ; \
	else \
		ln -s $(SCOREC_UTIL_DIR)/make_model bin ; \
	fi
	ln -s $(SRCDIR)/_$(ARCH)/a2cc bin
	ln -s $(SRCDIR)/_$(ARCH)/mackenbach_profiles bin
	ln -s $(SRCDIR)/_$(ARCH)-opt-25/m3dc1_2d bin
	ln -s $(SRCDIR)/_$(ARCH)-opt-25-omp/m3dc1_2d_omp bin
	ln -s $(SRCDIR)/_$(ARCH)-complex-opt-25/m3dc1_2d_complex bin
	ln -s $(SRCDIR)/_$(ARCH)-complex-opt-25-omp/m3dc1_2d_complex_omp bin
	ln -s $(SRCDIR)/_$(ARCH)-3d-opt-60/m3dc1_3d bin
	ln -s $(SRCDIR)/_$(ARCH)-3d-opt-60-omp/m3dc1_3d_omp bin
	ln -s $(SRCDIR)/_$(ARCH)-3d-st-opt-125/m3dc1_3d_st bin
	ln -s $(SRCDIR)/_$(ARCH)-3d-st-opt-125-omp/m3dc1_3d_omp_st bin
	ln -s $(SRCDIR)/sbin/extract_profiles.sh bin
	ln -s $(SRCDIR)/sbin/m3dc1_units.sh bin
	ln -s $(SRCDIR)/sbin/create_mesh.sh bin
	ln -s $(SRCDIR)/sbin/create_fixed_mesh.sh bin
	ln -s $(SRCDIR)/sbin/part_mesh.sh bin

endif
