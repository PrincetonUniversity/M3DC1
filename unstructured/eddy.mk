FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DLATESTSCOREC -DUSEBLAS -DPETSC_VERSION=37
CCOPTS  = -c -DPETSC_VERSION=37

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
PETSCVER=petsc3.7.6
PETSC_VER=petsc-3.7.6
PETSC_DIR=/home/jinchen/LIB/$(PETSC_VER)
ifeq ($(COM), 1)
PETSC_ARCH=cplx-intel18.0-openmpi3.0.0
else
PETSC_ARCH=real-intel18.0-openmpi3.0.0
endif

MPI_DIR=/usr/local/openmpi/3.0.0/intel180/x86_64
INTEL_DIR=/opt/intel/compilers_and_libraries_2018.3.222/linux
GCC_DIR=/usr/lib/gcc/x86_64-redhat-linux/4.8.5

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/usr/local/intel/lib64 -L/usr/local/intel/lib64/openmpi -L$(MPI_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L$(INTEL_DIR)/tbb/lib/intel64/gcc4.7 -L$(INTEL_DIR)/compiler/lib/intel64 -L$(INTEL_DIR)/mkl/lib/intel64 -L$(INTEL_DIR)/compiler/lib/intel64_lin -L$(GCC_DIR) -Wl,-rpath,/usr/local/intel/lib64 -Wl,-rpath,/usr/local/intel/lib64/openmpi -Wl,-rpath,$(MPI_DIR)/lib64 -lpetsc -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lsuperlu -lscalapack -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhwloc -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lifport -lifcoremt_pic -lrt -lm -lpthread -lz -ldl -L/usr/local/intel/lib64 -L/usr/local/intel/lib64/openmpi -L$(MPI_DIR)/lib64 -lmpi -L/usr/local/intel/lib64 -L/usr/local/intel/lib64/openmpi -L$(MPI_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L$(INTEL_DIR)/tbb/lib/intel64/gcc4.7 -L$(INTEL_DIR)/compiler/lib/intel64 -L$(INTEL_DIR)/mkl/lib/intel64 -L/usr/local/intel/lib64 -L$(INTEL_DIR)/compiler/lib/intel64_lin -L$(GSL_ROOT_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L/usr/local/intel/lib64 -L/usr/local/intel/lib64 -L$(GCC_DIR) -L$(GSL_ROOT_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L$(INTEL_DIR)/tbb/lib/intel64/gcc4.7 -L$(INTEL_DIR)/compiler/lib/intel64 -L$(INTEL_DIR)/mkl/lib/intel64 -L/usr/local/intel/lib64 -L/usr/local/intel/lib64 -Wl,-rpath,/usr/local/intel/lib64 -Wl,-rpath,/usr/local/intel/lib64/openmpi -Wl,-rpath,$(MPI_DIR)/lib64 -limf -lsvml -lirng -lm -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lpthread -lirc_s -L/usr/local/intel/lib64 -L/usr/local/intel/lib64/openmpi -L$(MPI_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L$(INTEL_DIR)/tbb/lib/intel64/gcc4.7 -L$(INTEL_DIR)/compiler/lib/intel64 -L$(INTEL_DIR)/mkl/lib/intel64 -L/usr/local/intel/lib64 -L$(INTEL_DIR)/compiler/lib/intel64_lin -L$(GSL_ROOT_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L/usr/local/intel/lib64 -L/usr/local/intel/lib64 -L$(GCC_DIR) -L$(GSL_ROOT_DIR)/lib64 -L$(GSL_ROOT_DIR)/lib64 -L$(INTEL_DIR)/tbb/lib/intel64/gcc4.7 -LL$(INTEL_DIR)/compiler/lib/intel64 -L$(INTEL_DIR)/mkl/lib/intel64 -L/usr/local/intel/lib64 -L/usr/local/intel/lib64 -ldl

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_BASE_DIR=/home/jinchen/LIB/scorec/intel18.0-openmpi3.3.0/
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(PETSCVER)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(PETSCVER)
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
