ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpiicc
  CC = mpiicc
  F90 = mpiifort -xMIC-AVX512
  F77 = mpiifort -xMIC-AVX512
  LOADER = mpiifort -xMIC-AVX512
endif

OPTS := $(OPTS) -DUSEBLAS

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

SCOREC_DIR = /global/project/projectdirs/mp288/carl/scorec/Nov2016-impi5.1.3

ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif
endif

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

ifeq ($(TRILINOS),1)
TRILINOS_LIBS = -Wl,--start-group,-rpath,$(CRAY_TRILINOS_PREFIX_DIR)/lib -L$(CRAY_TRILINOS_PREFIX_DIR)/lib \
                -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm_intel -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
else
TRILINOS_LIBS =
endif

ifeq ($(COM), 1)
      PETSC_DIR = /project/projectdirs/mp288/carl/petsc-3.5.4
      PETSC_ARCH = intel_complex
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_dist_3.3 -lparmetis -lmetis
      PETSC_LIB =  -lpetsc
      ZOLTAN_LIB =
else
      PETSC_DIR = /project/projectdirs/mp288/carl/petsc-3.5.4
      PETSC_ARCH = intel_real
ifeq ($(TRILINOS), 1)
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lparmetis -lmetis
      PETSC_LIB =
else
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lesmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu_dist -lparmetis -lmetis -lpthread -lssl 
      PETSC_LIB = -lpetsc
endif
      ZOLTAN_LIB =
endif

ifeq ($(USEADIOS), 1)
  OPTS := $(OPTS) -DUSEADIOS
  ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L/usr/common/usg/minixml/2.7/lib -lm -lmxml \
             -L/usr/lib64/ -llustreapi
else
  ADIOS_FLIB =
endif

HDF5_DIR := /global/project/projectdirs/mp288/carl/hdf5-1.8.18

INCLUDE := $(INCLUDE) -I$(HDF5_DIR)/include -I$(MKLROOT)/include/fftw \
        -I$(SCOREC_DIR)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include/petsc \
        -I$(PETSC_DIR)/include \
	-I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include

LIBS := $(LIBS) \
        $(ZOLTAN_LIB) \
        $(TRILINOS_LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(PETSC_LIB) \
        $(PETSC_EXTERNAL_LIB_BASIC) -lstdc++ \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz \
        -L$(GSL_DIR)/lib -lgsl \
	-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 \
        $(ADIOS_FLIB)

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)  -O3
  CCOPTS := $(CCOPTS) -O3
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -check noarg-temp-created -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)  
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -openmp -Wl,-z,muldefs 
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
