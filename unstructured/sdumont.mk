ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpiicpc 
  CC = mpiicc 
  F90 = mpiifort 
  F77 = mpiifort 
  LOADER = mpiifort 
endif

#NEWSOLVERDEVELOPMENT needs more tests.
ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=m3dc1_scorec_complex
  PETSC_ARCH=cplx-intel-psxe2019
else
  M3DC1_SCOREC_LIB=m3dc1_scorec
  PETSC_ARCH=real-intel-psxe2019
endif

ifeq ($(PETSCVER), 37)
  PETSCVER=37
  PETSC_DIR=/scratch/ntm/software/petsc/petsc-3.7.6
  SCOREC_DIR=/scratch/ntm/software/scorec/intel-psxe2019/petsc3.7.6
  PETSC_WITH_EXTERNAL_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group -L$(MKLROOT)/lib/intel64_lin -Wl,--end-group -lpetsc -lsuperlu_dist -lparmetis -lmetis -lsuperlu -lmkl_blas95_lp64 -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -lfftw3_mpi -lfftw3 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lrt -lm -lpthread -lz -ldl -lstdc++
else
  PETSCVER=39
  PETSC_DIR=/scratch/ntm/software/petsc/petsc-3.9.3
  SCOREC_DIR=/scratch/ntm/software/scorec/intel-psxe2019/petsc3.9.3
  PETSC_WITH_EXTERNAL_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group -L$(MKLROOT)/lib/intel64_lin -Wl,--end-group -lpetsc -lmkl_blas95_lp64 -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lrt -lm -lpthread -lz -ldl -lstdc++
endif

GSL_DIR=/scratch/app/gsl/2.5_gnu
ZOLTAN_DIR=$(SCOREC_DIR)
SCOREC_UTIL_DIR=$(SCOREC_DIR)/bin

PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_LIB = -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           -l$(M3DC1_SCOREC_LIB) $(PUMI_LIB) -Wl,--end-group

LIBS =  $(SCOREC_LIB) \
        -L$(ZOLTAN_DIR)/lib -lzoltan \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas

INCLUDE = -I$(SCOREC_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(GSL_DIR)/include

# Optimization flags
FOPTS = -c -r8 -implicitnone -fpp -warn all \
        -DPETSC_VERSION=$(PETSCVER) $(OPTS)
CCOPTS  = -c \
        -DPETSC_VERSION=$(PETSCVER) $(OPTS)
ifeq ($(OPT), 0)
  FOPTS := $(FOPTS) -g -O0 -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
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
