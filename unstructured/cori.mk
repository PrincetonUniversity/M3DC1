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
 
OPTS := $(OPTS) -DUSEADIOS -DPETSC_VERSION=39 -DUSEBLAS #-DNEWSOLVERDEVELOPMENT

SCOREC_BASE_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.7.6/hsw-petsc3.9.3
SCOREC_UTIL_DIR=/global/project/projectdirs/mp288/cori/scorec/mpich7.7.6/hsw-bin

ifeq ($(REORDERED), 1)
  SCORECVER=reordered
endif

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifdef SCORECVER
  PETSC_DIR=/global/project/projectdirs/mp288/cori/petsc/petsc-3.9.3
  ifeq ($(COM), 1)
    PETSC_ARCH=cplx-intel-mpi7.7.6-hsw
  else
    PETSC_ARCH=real-intel-mpi7.7.6-hsw
  endif  
  MKL_LIB =
  PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,/opt/cray/pe/hdf5-parallel/1.10.2.0/INTEL/16.0/lib -L/opt/cray/pe/hdf5-parallel/1.10.2.0/INTEL/16.0/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lstrumpack -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lrt -lm -lpthread -lz -ldl -lstdc++
else
  PETSC_DIR=/global/homes/j/jinchen/project/PETSC/petsc-3.9.3
  ifeq ($(COM), 1)
    PETSC_ARCH=cori-hsw-mpich776-cplx-nomkl-510
  else
    PETSC_ARCH=cori-hsw-mpich776-real-nomkl-510
  endif

  MKL_LIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
  PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lrt -lm -lpthread -lz -ldl -lstdc++
endif
#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
#ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
             -L/usr/lib64/ -llustreapi

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        -L$(SCOREC_DIR)/lib -l$(M3DC1_SCOREC_LIB) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB) \
	$(MKL_LIB)
#        $(HYBRID_LIBS) \

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -static -qopt-report
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
