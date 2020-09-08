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
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif
 
OPTS := $(OPTS) -DNOTUSEADIOS -DPETSC_VERSION=39 -DUSEBLAS #-DNEWSOLVERDEVELOPMENT

PETSC_DIR=/global/homes/j/jinchen/project/PETSC/petsc
ifeq ($(COM), 1)
  PETSC_ARCH=corigpu-pgi1910-mvapich2.232-cuda10289-master-cplx2
else
  PETSC_ARCH=corigpu-pgi1910-mvapich2.232-cuda10289-master-real2
endif
PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lsuperlu -lsuperlu_dist -lparmetis -lmetis -lrt -lm -lpthread -lz -ldl -lstdc++ 


FFTW_LIB=-L$(FFTW_DIR) -lfftw3_mpi -lfftw3 

SCOREC_DIR=/global/project/projectdirs/mp288/cori/scorec/mvapich2.3.2/cuda10.2-pgi19.10-petsc3.12.4
SCOREC_UTIL_DIR=$(SCOREC_DIR)/bin

ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

ZOLTAN_DIR=$(SCOREC_DIR)
ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lzoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
#ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
#ADIOS_DIR=/global/homes/j/jinchen/project/LIB/adios-1.13.0/build-mpi
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L/usr/common/software/minixml/2.9/hsw//gnu/lib -lm -lmxml \
             -L/usr/lib64/ -llustreapi

MKL_LIB = -L/usr/common/software/pgi/linux86-64/18.10/lib -lblas -llapack #-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

CUDA_LIB:=-L/usr/common/software/cuda/10.2.89/lib64 -lcusolver_static -lcublas_static -lcusparse_static -lcudart_static -lcudadevrt -lcublasLt_static -lculibos -L/usr/common/software/pgi/19.10/linux86-64-llvm/19.10/lib -lblas -lcudaforblas -llapack
#-L/usr/local/cuda-10.2/lib64 -lcudart_static -lcusparse_static -lcusolver_static -lculibos 
#-lcufft_static

#-lcufftw_static
#-lcufft_static_nocallback
#-lcurand_static
#-lcudadevrt

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(FFTW_INC) \
	   -I$(GSL_DIR)/include  \
           -I$(HDF5_DIR)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	$(CUDA_LIB) \
	$(FFTW_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
        -lz \
	-lX11

FOPTS = -c -r8 -Mpreprocess $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) #-static
  FOPTS  := $(FOPTS) -Kieee
  CCOPTS := $(CCOPTS)
else
  FOPTS := $(FOPTS) -g -Mbounds #-check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -fopenmp 
  FOPTS  := $(FOPTS)  -fopenmp  -Kieee
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

%.o: %.F90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
