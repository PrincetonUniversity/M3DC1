ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpic++ 
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif
 
OPTS := $(OPTS) -DPETSC_VERSION=37 -DUSEBLAS #-DNEWSOLVERDEVELOPMENT

SCOREC_BASE_DIR=/gpfs/wolf/gen127/proj-shared/scorec/gcc8.1-cuda10.1-mpi110.3
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ifeq ($(REORDERED), 1)
  SCOREC_DIR=$(SCOREC_BASE_DIR)/reordered
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

#ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

PETSC_DIR=/gpfs/wolf/gen127/proj-shared/petsc/petsc-3.7.6
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-gcc-cuda-mpi10.3.0
else
  PETSC_ARCH=real-gcc-cuda-mpi10.3.0
endif


PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(OLCF_PARMETIS_ROOT)/lib -L$(OLCF_PARMETIS_ROOT)/lib -Wl,-rpath,$(OLCF_METIS_ROOT)/lib -L$(OLCF_METIS_ROOT)/lib -Wl,-rpath,$(OLCF_NETLIB_SCALAPACK_ROOT)/lib -L$(OLCF_NETLIB_SCALAPACK_ROOT)/lib -Wl,-rpath,$(OLCF_NETLIB_LAPACK_ROOT)/lib64 -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64 -lpetsc -lsuperlu_dist -lparmetis -lmetis -lsuperlu -lscalapack -llapack -lblas -lX11 -lhwloc -lm -ldl -lstdc++

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(OLCF_GSL_ROOT)/include # \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(OLCF_HDF5_ROOT)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
        -L$(OLCF_FFTW_ROOT)/lib -lfftw3_mpi -lfftw3 \
	-L$(OLCF_GSL_ROOT)/lib -lgsl -lgslcblas #-lhugetlbfs 
#        $(HYBRID_LIBS) \

#FIXME: please check if FOPTS is correct
# See https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html#Option-Summary
FOPTS = -c -fdefault-real-8 -fimplicit-none -cpp -Wall$(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
# FIXME
#  LDOPTS := $(LDOPTS) -static -qopt-report
#  FOPTS  := $(FOPTS)  -qopt-report
#  CCOPTS := $(CCOPTS) -qopt-report
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -fpe0 -traceback -debug extended
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
