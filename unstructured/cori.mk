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

OPTS := $(OPTS) -DUSEADIOS -DPETSC_VERSION=39 -DUSEBLAS

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609
ifeq ($(COM), 1)
  PETSC_ARCH=corihsw-PrgEnvintel6010-craympich7719-master-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-cplx/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lgsl -lgslcblas -lquadmath -lstdc++ -ldl
M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
  M3DC1_SCOREC_LIBA=libm3dc1_scorec_complex.a
  SCOREC_COMPLEX=ON
else
  ifeq ($(ST), 1)
    PETSC_ARCH=corihsw-PrgEnvintel6010-craympich7719-master-real-st
    PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-real-st/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-real-st/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lz -lgsl -lgslcblas -lquadmath -lstdc++ -ldl
  else
    PETSC_ARCH=corihsw-PrgEnvintel6010-craympich7719-master-real
    PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-real/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/corihsw-PrgEnvintel6010-craympich7719-master-real/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lgsl -lgslcblas -lquadmath -lstdc++ -ldl
  endif
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  M3DC1_SCOREC_LIBA=libm3dc1_scorec.a
  SCOREC_COMPLEX=OFF
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core/upgrade-intel6610-craympich7719-hsw2
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

MKL_LIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L$(ADIOS_DIR)/src/mxml -lm -lmxml \
             -L/usr/lib64/ -llustreapi

INCLUDE := $(INCLUDE) -I$(SCOREC_BASE_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include

LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
        $(SCOREC_LIB) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB) \
	$(MKL_LIB)

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L/global/homes/j/jinchen/project/NETCDF/buildhsw/lib -Wl,-rpath,/global/homes/j/jinchen/project/NETCDF/buildhsw/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lnetcdf -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I/global/cfs/cdirs/mp288/jinchen/NETCDF/buildhsw/include
endif

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
  FOPTS := $(FOPTS) -g -Mbounds -check noarg_temp_created -fpe0 -warn -traceback -debug extended
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
