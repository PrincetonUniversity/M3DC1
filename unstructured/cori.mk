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

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

SCOREC_UTIL_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/bin
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/Aug2017/
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/Nov2017/
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/debug
SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/Dec2017/
ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
    ZOLTAN_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.8.2/cori-hsw-mpich760-cplx
    ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif
    ZOLTAN_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.8.2/cori-hsw-mpich760-real
    ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan
endif

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

ifeq ($(COM), 1)
      PETSC_DIR = /global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.8.2
      PETSC_ARCH = cori-hsw-mpich760-cplx
      HYPRE_LIB = 
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(HYPRE_LIB) \
       $(HYPRE_LIB) \
       -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu -lsuperlu_dist \
       -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl -lstdc++  \
       -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis -lscotch -lscotcherr -lscotcherrexit \
       -lflapack -lfblas
      PETSC_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
      OPTS := $(OPTS) -DPETSC_VERSION=38
else
      PETSC_DIR = /global/homes/j/jinchen/project/PETSC/petsc-3.8.2
      PETSC_ARCH = cori-hsw-mpich760-real
      HYPRE_LIB = -lHYPRE
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
        $(HYPRE_LIB) \
       -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu -lsuperlu_dist \
       -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl -lstdc++  \
       -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis -lscotch -lscotcherr -lscotcherrexit \
       -lflapack -lfblas 
#       -lstrumpack_sparse \

      PETSC_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
      OPTS := $(OPTS) -DPETSC_VERSION=38
endif

# Include option to use adios
OPTS := $(OPTS) -DUSEADIOS

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
#ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
             -L/usr/lib64/ -llustreapi


OPTS := $(OPTS) -DUSEBLAS -DNEWSOLVERDEVELOPMENT

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
           $(FFTW_INCLUDE_OPTS) \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIB) $(PETSC_EXTERNAL_LIB_BASIC) \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz \
	$(FFTW_POST_LINK_OPTS) -lfftw3 \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB)
#        $(HYBRID_LIBS) \

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
        -Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -dynamic -qopt-report
  FOPTS  := $(FOPTS)  -qopt-report
  CCOPTS := $(CCOPTS) -qopt-report
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -fopenmp 
  FOPTS  := $(FOPTS)  -fopenmp 
  CCOPTS := $(CCOPTS) -fopenmp 
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
