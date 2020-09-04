FOPTS = -c -r8 -implicitnone -fpp -warn all -DNEXTPetscDEV -DPETSC_VERSION=38 $(OPTS)
CCOPTS  = -c -O -DPETSC_VERSION=38 -DNEXTPetscDEV

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 # -fast
#  FOPTS  := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
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
  LOADER = mpif90 -cxxlib
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define petsc and other libs
PETSC_DIR=/projects/M3DC1/PETSC/petsc
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-intelmpi
else
  PETSC_ARCH=real-intelmpi
endif

MPI_DIR=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/lib64
MKL_DIR=/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl
COMPILER_DIR=/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64
COMPILER_LIN_DIR=/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin

PETSC_WITH_EXTERNAL_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(MPI_DIR) -L$(COMPILER_DIR) -L$(MKL_DIR) -L$(COMPILER_LIN_DIR) -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -Wl,-rpath,$(MPI_DIR) -lpetsc -lsuperlu -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack  -lzoltan -L$(FFTW3DIR) -lfftw3_mpi -lfftw3 -L$(MKLROOT)/lib/intel64_lin -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lmkl_cdft_core -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lX11 -lifport -lifcoremt_pic -lmpicxx -lintlc -lrt -lmpi -lmpigi -lrt -lpthread -lm -lpthread -lz -ldl -lstdc++ -lmpi -lmpifort -limf -lsvml -lirng -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lpthread -lirc_s -lifport -lifcoremt_pic -lmpicxx -lintlc -lrt -lm -lpthread -lz -ldl -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/lib/release_mt -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/lib -lmpifort -lmpi -lmpigi -lrt -lpthread

PETSCVER=petsc3.12.4
SCOREC_BASE_DIR=/projects/M3DC1/scorec/intel18.0-mpi2018.3/petsc3.12.4/
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
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

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

INCLUDE = -I$(SCOREC_DIR)/include \
          -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

LIBS = $(SCOREC_LIBS) \
       $(ZOLTAN_LIB) \
       $(PETSC_WITH_EXTERNAL_LIB) \
       -L$(GSL_ROOT_DIR)/lib64  -lssl -lgsl -lgslcblas

LDOPTS := $(LDOPTS) -lmpicxx -lmpifort -limf -lmpigi -lrt -lpthread $(PETSC_KSP_LIB)
FOPTS  := $(FOPTS)  
CCOPTS := $(CCOPTS) 

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
