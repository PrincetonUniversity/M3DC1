FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DUSEBLAS -DPETSC_VERSION=312
CCOPTS  = -c -DPETSC_VERSION=312 -DSTELLAR

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -g -check-uninit -debug all
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=../select.tau
  CC     = tau_cc.sh $(TAU_OPTIONS)
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CC = mpicc
  CPP = mpicxx
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90
  LDOPTS := $(LDOPTS) -cxxlib
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
PETSCVER=petsc3.12.0
PETSC_VER=petsc-3.12.0
PETSC_DIR=/projects/M3DC1/PETSC/petsc.stellar
ifeq ($(COM), 1)
PETSC_ARCH=cplx-stella-openmpi-2021
else
PETSC_ARCH=real-stella-openmpi-2021
endif

MPI_DIR=/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi
INTEL_DIR=/opt/intel/compilers_and_libraries_2018.3.222/linux
GCC_DIR=/usr/lib/gcc/x86_64-redhat-linux/4.8.5

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH} -L/usr/local/openmpi/4.1.0/intel20211/lib64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/4.1.0/intel20211/lib64 -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lflapack -lfblas -lzoltan -lparmetis -lmetis -lz -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_BASE_DIR=${PETSC_DIR}/${PETSC_ARCH}
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

LIBS = 	\
	$(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW3DIR)/lib -lfftw3 -lfftw3f_mpi -lfftw3l_mpi -lfftw3_mpi \
	-L$(HDF5DIR)/lib64 -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran \
	-L$(GSL_ROOT_DIR)/lib64 -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5DIR)/include \
        -I$(GSL_ROOT_DIR)/include

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP) $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
