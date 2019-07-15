ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = CC
  CC = cc
  F90 = ftn
  F77 = ftn
  LOADER = ftn
endif

OPTS := $(OPTS) -xMIC-AVX512 -DPETSC_VERSION=39 -DUSEADIOS -DUSEBLAS

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

PETSC_DIR = /global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.9.3
ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB = m3dc1_scorec_complex
  PETSC_ARCH = cori-knl-mpich773-cplx-nomkl-510
else
  M3DC1_SCOREC_LIB = m3dc1_scorec
  PETSC_ARCH = cori-knl-mpich773-real-nomkl-510
endif

PETSC_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
     -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
     -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lptesmumps \
     -lpord -lsuperlu -lsuperlu_dist -lstrumpack \
     -lparmetis -lmetis -lpthread -ldl -lstdc++  \
     -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis -lscotch -lscotcherr -lscotcherrexit

ifeq ($(REORDERED), 1)
  SCORECVER=reordered
endif

SCOREC_BASE_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.7.3/knl-petsc3.9.3
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)

ZOLTAN_LIB=-L$(SCOREC_DIR)/lib -lzoltan

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
             -L/usr/lib64/ -llustreapi

FFTW_DIR=/opt/cray/pe/fftw/3.3.8.2/mic_knl
GSL_DIR=/usr/common/software/gsl/2.1/intel
HDF5_DIR=/opt/cray/pe/hdf5-parallel/1.10.2.0/INTEL/16.0

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
           -I$(FFTW_DIR)/include \
           -I$(HDF5_DIR)/include \
           -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
           -I$(GSL_DIR)/include

LIBS := $(LIBS) \
        -L$(SCOREC_DIR)/lib -l$(M3DC1_SCOREC_LIB) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB)\
        $(PETSC_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz \
        -L$(FFTW_DIR)/lib -lfftw3 \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas -lhugetlbfs \
        $(ADIOS_FLIB)

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2 -O2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2 -O2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2 -O2
else

# Optimization flags
ifeq ($(OPT), 1)
#  LDOPTS := $(LDOPTS) -dynamic -ipo -qopt-report  -qopt-report-phase=vec #-h profile_generate 
#  FOPTS  := $(FOPTS)  -O3 -ipo -qopt-report  -qopt-report-phase=vec #-h profile_generate 
#  CCOPTS := $(CCOPTS) -O3 -ipo -qopt-report  -qopt-report-phase=vec #-h profile_generate 
  LDOPTS := $(LDOPTS) -static -qopt-report=5 -qopt-report-phase=vec,loop
  FOPTS  := $(FOPTS)  -qopt-report=5 -qopt-report-phase=vec,loop
  CCOPTS := $(CCOPTS) -qopt-report=5 -qopt-report-phase=vec,loop
else
  LDOPTS := $(LDOPTS) -static
  FOPTS := $(FOPTS) -g -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -openmp 
  FOPTS  := $(FOPTS)  -openmp 
  CCOPTS := $(CCOPTS) -openmp 
endif

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
