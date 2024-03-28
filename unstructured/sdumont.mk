FOPTS = -c -r8 -implicitnone -fpp -warn all -DPETSC_VERSION=318 $(OPTS)
CCOPTS  = -c -DPETSC_VERSION=318 $(OPTS)

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

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = /scratch/app/openmpi/icc/4.1.4/bin/mpicxx
  CC = /scratch/app/openmpi/icc/4.1.4/bin/mpicc 
  F90 = /scratch/app/openmpi/icc/4.1.4/bin/mpif90
  F77 = /scratch/app/openmpi/icc/4.1.4/bin/mpif90
  LOADER = /scratch/app/openmpi/icc/4.1.4/bin/mpif90
endif

#NEWSOLVERDEVELOPMENT needs more tests.
ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

MPIVER=intel-psxe2020-openmpiicc4.1.4
PETSCVER=3.18.2
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$(PETSCVER)
GSL_DIR=/scratch/app/gsl/2.7_intel_2020

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
  PETSC_ARCH=cplx-$(MPIVER)
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_ARCH=real-$(MPIVER)
endif

SCOREC_BASE_DIR=/scratch/ntm/software/scorec/$(MPIVER)/petsc-$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/scratch/app/openmpi/icc/4.1.4/lib -L/scratch/app/ucx/1.13/lib -L/scratch/app/xpmem/2.6.5/lib -L/scratch/app/openpmix/4.2.1/lib -L/opt/intel/parallel_studio_xe_2020/clck/2019.10/lib/intel64 -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64 -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8 -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.4 -L/opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.8 -L/scratch/app/gcc/9.3/lib/gcc/x86_64-pc-linux-gnu/9.3.0 -L/scratch/app/gcc/9.3/lib64 -L/scratch/app/gcc/9.3/lib -Wl,-rpath,/scratch/app/openmpi/icc/4.1.4/lib -Wl,-rpath,/scratch/app/ucx/1.13/lib -Wl,-rpath,/scratch/app/xpmem/2.6.5/lib -Wl,-rpath,/scratch/app/openpmix/4.2.1/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -ldl -lstdc++ -lmpi_cxx -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -ldl -lstdc++

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

LIBS =  $(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas

INCLUDE = -I$(SCOREC_DIR)/include \
        -I$(PUMI_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(GSL_DIR)/include

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
