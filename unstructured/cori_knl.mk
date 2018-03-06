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

OPTS := $(OPTS) -xMIC-AVX512 -DUSEBLAS -DPETSC_VERSION=38 -DNEWSOLVERDEVELOPMENT

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB = m3dc1_scorec_complex
  PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.8.2/
  PETSC_ARCH=cori-hsw-mpich760-cplx
  HYPRE_LIB = 
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif
  ifeq ($(OMP), 1)
    PETSC_DIR=/global/homes/j/jinchen/project/PETSC/master
    PETSC_ARCH=cori-hsw-knl-mpich760-omp-strumpack
    OPTS := $(OPTS) -DSTRUMPACK
    STRUMPACK_LIB = -lstrumpack_sparse
  else
    PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.8.2
    PETSC_ARCH = cori-hsw-mpich760-real
  endif
  HYPRE_LIB = -lHYPRE
endif

PETSC_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
     -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
     $(HYPRE_LIB) \
     -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lptesmumps \
     -lpord -lsuperlu -lsuperlu_dist $(STRUMPACK_LIB) \
     -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl -lstdc++  \
     -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis -lscotch -lscotcherr -lscotcherrexit #\
#       -lflapack -lfblas
#       -lstrumpack_sparse \


SCOREC_UTIL_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/knl/Aug2017/bin
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/knl/Aug2017
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/knl/Nov2017
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/haswell/debug
#SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/knl/Dec2017/
SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.6.0/knl/March2018

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

# Include option to use ADIOS
OPTS := $(OPTS) -DUSEADIOS
#
##only define them if adios-1.3 is used; otherwise use hopper default
##ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
##ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
##ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB_V1 = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
#             -L/usr/lib64/ -llustreapi

MKL_LIB = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := \
        $(LIBS) \
        $(SCOREC_LIBS) \
        -lzoltan \
        $(PETSC_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB_V1) \
	$(MKL_LIB)
#        $(HYBRID_LIBS) \

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
        -Dglobalinsertval=insertval -Dglobalentdofs=entdofs
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
  LDOPTS := $(LDOPTS) -dynamic -qopt-report=5 -qopt-report-phase=vec,loop
  FOPTS  := $(FOPTS)  -qopt-report=5 -qopt-report-phase=vec,loop
  CCOPTS := $(CCOPTS) -qopt-report=5 -qopt-report-phase=vec,loop
else
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
