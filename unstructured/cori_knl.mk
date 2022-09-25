addifeq ($(TAU), 1)
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

OPTS := $(OPTS) -xMIC-AVX512 -DUSEBLAS -DPETSC_VERSION=39

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609
ifeq ($(COM), 1)
  PETSC_ARCH = coriknl-PrgEnvintel6010-craympich7719-master-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/coriknl-PrgEnvintel6010-craympich7719-master-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/coriknl-PrgEnvintel6010-craympich7719-master-cplx/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lquadmath -lstdc++ -ldl
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
  M3DC1_SCOREC_LIBA=libm3dc1_scorec_complex.a
  SCOREC_COMPLEX=ON
else
  PETSC_ARCH = coriknl-PrgEnvintel6010-craympich7719-master-real
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/coriknl-PrgEnvintel6010-craympich7719-master-real/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/coriknl-PrgEnvintel6010-craympich7719-master-real/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lquadmath -lstdc++ -ldl
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  M3DC1_SCOREC_LIBA=libm3dc1_scorec.a
  SCOREC_COMPLEX=OFF
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core/upgrade-intel6610-craympich7719-knl2
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

MKL_LIB =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

INCLUDE += -I$(SCOREC_BASE_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include

LIBS += $(SCOREC_LIB) \
        $(ZOLTAN_LIB)\
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(MKL_LIB)

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L/global/homes/j/jinchen/project/NETCDF/buildhsw/lib -Wl,-rpath,/global/homes/j/jinchen/project/NETCDF/buildhsw/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lnetcdf -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I/global/cfs/cdirs/mp288/jinchen/NETCDF/buildhsw/include
else
  LIBS += -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz
endif

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
  LDOPTS := $(LDOPTS) -static -qopt-report=5 -qopt-report-phase=vec,loop
  FOPTS  := $(FOPTS)  -qopt-report=5 -qopt-report-phase=vec,loop
  CCOPTS := $(CCOPTS) -qopt-report=5 -qopt-report-phase=vec,loop
else
  LDOPTS := $(LDOPTS) -static
  FOPTS := $(FOPTS) -g -Mbounds -check noarg_temp_created -fpe0 -warn -traceback -debug extended
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
