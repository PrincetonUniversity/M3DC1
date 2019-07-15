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
  LOADER = ftn -static
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

OPTS := $(OPTS) -DPETSC_VERSION=37 -DUSEBLAS

SCOREC_UTIL_DIR=/global/project/projectdirs/mp288/edison/scorec/mpich7.6.0/bin

SCOREC_DIR = /global/project/projectdirs/mp288/edison/scorec/mpich7.6.0/July2017
ZOLTAN_LIB = -L/project/projectdirs/mp288/jinchen/PETSC/petsc-3.7.7/edison-intel1801163-mpich773/lib -lzoltan

ifeq ($(TRILINOS), 1)
    TRILINOS_LIBS = -Wl,--start-group,-rpath,$(CRAY_TRILINOS_PREFIX_DIR)/lib -L$(CRAY_TRILINOS_PREFIX_DIR)/lib \
                -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm_intel -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
    PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(CRAY_TPSL_PREFIX_DIR)/lib -L$(CRAY_TPSL_PREFIX_DIR)/lib \
                -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_trilinos
else
  ifeq ($(COM), 1)
    HYPRE_LIB=
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
    PETSC_DIR=/opt/cray/pe/petsc/3.7.6.0/complex/INTEL/16.0/sandybridge
    PETSC_LIB = -lcraypetsc_intel_complex
  else
    HYPRE_LIB=-lHYPRE
    M3DC1_SCOREC_LIB=-lm3dc1_scorec
   #PETSC_DIR=/opt/cray/pe/petsc/3.7.6.0/real/INTEL/16.0/sandybridge
   #PETSC_LIB = -lcraypetsc_intel_real
    PETSC_DIR=/project/projectdirs/mp288/jinchen/PETSC/petsc-3.7.7
    PETSC_ARCH=edison-intel1801163-mpich773
    PETSC_LIB = -lpetsc
  endif

  PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
                $(PETSC_LIB) $(HYPRE_LIB) \
                -lsuperlu -lsuperlu_dist -lflapack -lfblas -lscalapack \
                -lcmumps -ldmumps -lesmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord \
                -lparmetis -lmetis -lptscotch -lscotch -lptscotcherr -lscotcherr \
                -lpthread -ldl
endif

SCOREC_LIBS=-Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
            $(M3DC1_SCOREC_LIB) -lpumi -lcrv -lph -lsam -lspr -lma \
            -lapf_zoltan -lparma -lmds -lapf -llion -lmth -lgmi -lpcu \
            -Wl,--end-group

# Use Adios
OPTS := $(OPTS) -DUSEADIOS
ADIOS_FLIB = ${ADIOSREAD_FLIB_V1} \
	-L/usr/lib64/ -llustreapi

INCLUDE := $(INCLUDE) $(FFTW_INCLUDE_OPTS) \
        -I$(SCOREC_DIR)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_EXTERNAL_LIB_BASIC) \
        -lhdf5_fortran -lhdf5 -lz \
        $(FFTW_POST_LINK_OPTS) -lfftw3 \
        -L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
        $(ADIOS_FLIB) \
	-lstdc++

ifeq ($(TRILINOS),1)
  LIBS := $(LIBS) $(TRILINOS_LIBS) \
         -L/usr/lib64/gcc/x86_64-suse-linux/4.3 -lstdc++
endif

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)
CCOPTS = -c $(OPTS)

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)
  CCOPTS := $(CCOPTS)
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
