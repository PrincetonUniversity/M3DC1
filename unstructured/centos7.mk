FOPTS = $(OPTS) -DPETSC_VERSION=313 -c -r8 -implicitnone -fpp -warn all -DUSEBLAS
CCOPTS  = -c -DPETSC_VERSION=313

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

ifeq ($(OMP), 1)
  FOPTS := $(FOPTS) -qopenmp
  CCOPTS := $(CCOPTS) -qopenmp
  LDOPTS := $(LDOPTS) -qopenmp
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
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver
PETSC_VER=petsc-3.13.5
PETSCVER=petsc3.13.5

PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/PETSC/petsc-3.13.5
ifeq ($(COM), 1)
PETSC_ARCH=cplx-rhel7-intel2019u3-openmpi4.0.3
else
PETSC_ARCH=real-rhel7-intel2019u3-openmpi4.0.3
endif

SCOREC_BASE_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel7/intel2019u3-openmpi4.0.3/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,/usr/pppl/intel/2019-pkgs/openmpi-4.0.3-pkgs/fftw-3.3.8/lib -L/usr/pppl/intel/2019-pkgs/openmpi-4.0.3-pkgs/fftw-3.3.8/lib -Wl,-rpath,/usr/pppl/intel/2019-pkgs/openmpi-4.0.3-pkgs/hdf5-parallel-1.10.5/lib -L/usr/pppl/intel/2019-pkgs/openmpi-4.0.3-pkgs/hdf5-parallel-1.10.5/lib -L/usr/pppl/intel/2019-pkgs/openmpi-4.0.3/lib -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/mpi/intel64/libfabric/lib -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/ipp/lib/intel64 -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2019.u3/compilers_and_libraries_2019.3.199/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -Wl,-rpath,/usr/pppl/intel/2019-pkgs/openmpi-4.0.3/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lX11 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl

LIBS = 	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-lgsl \
	-lX11

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_HOME)/include 

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
