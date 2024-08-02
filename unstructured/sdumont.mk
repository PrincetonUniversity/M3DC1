FOPTS = -c -r8 -implicitnone -fpp -warn all -DPETSC_VERSION=39 $(OPTS)
CCOPTS  = -c -DPETSC_VERSION=321 $(OPTS)
ifeq ($(OPT), 0)
  FOPTS := $(FOPTS) -g -O0 -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -openmp
  FOPTS  := $(FOPTS)  -openmp
  CCOPTS := $(CCOPTS) -openmp
endif

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

MPIVER=intel-psxe2020
PETSCVER=3.21.2
PETSC_DIR=/scratch/ntm/software/scorec/petsc/petsc-$(PETSCVER)

SCOREC_BASE_DIR=/scratch/ntm/software/scorec/$(MPIVER)/petsc$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
  PETSC_ARCH=cplx-$(MPIVER)
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_ARCH=real-$(MPIVER)
endif

PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin -fortlib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -Wl,--start-group -lmkl_blas95_lp64 -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -Wl,--end-group -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -Wl,--start-group -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -Wl,--end-group -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lquadmath -lstdc++

# -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lX11 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz 
#
#PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/scratch/ntm/software/scorec/petsc/petsc-3.21.2/real-intel-psxe2020/lib -L/scratch/ntm/software/scorec/petsc/petsc-3.21.2/real-intel-psxe2020/lib -Wl,-rpath,/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin -fortlib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -Wl,--start-group -lmkl_blas95_lp64 -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -Wl,--end-group -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -Wl,--start-group -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -Wl,--end-group -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lquadmath -lstdc++

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

GSL_DIR=/scratch/app/gsl/2.7_intel_2020

LIBS =  $(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas


INCLUDE = -I$(SCOREC_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(GSL_DIR)/include

ifeq ($(ST), 1)
NETCDF_DIR=/scratch/app/netcdf/4.7_intel_2020
LIBS += -Wl,--start-group -L$(NETCDF_DIR)/lib -Wl,-rpath,$(NETCDF_DIR)/lib -lnetcdf -lnetcdff -lz -Wl,--end-group
INCLUDE += -I$(NETCDF_DIR)/include
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
