FOPTS = -c -r8 -implicitnone -fpp -warn all -DPETSC_VERSION=39 $(OPTS)
CCOPTS  = -c -DPETSC_VERSION=39 $(OPTS)
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
  CPP = mpicxx
  CC = mpicc 
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90
endif

#NEWSOLVERDEVELOPMENT needs more tests.
ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

MPIVER=intel-psxe2019-openmpiicc4.0.4
PETSCVER=3.9.4
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$(PETSCVER)

GSL_DIR=/scratch/app/gsl/2.5_gnu
SCOREC_DIR=/scratch/ntm/software/scorec/intel2019-openmpi4.0.4/petsc$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_DIR)/bin

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=m3dc1_scorec_complex
  PETSC_ARCH=cplx-$(MPIVER)
  MUMPS_LIBS=-lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord
else
  M3DC1_SCOREC_LIB=m3dc1_scorec
  PETSC_ARCH=real-$(MPIVER)
  MUMPS_LIBS=
endif

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(MPI_ROOT)/lib -L$(I_MPI_ROOT)/intel64/libfabric/lib -L/opt/intel/parallel_studio_xe_2019/intelpython3/lib/libfabric -L/opt/intel/parallel_studio_xe_2019/clck/2019.2.1/lib/intel64 -L$(IPPROOR)/lib/intel64 -L/opt/intel/parallel_studio_xe_2019/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64_lin -L$(MKLROOT)/lib/intel64_lin -L$(TBBROOT)/lib/intel64/gcc4.7 -L$(DAALROOT)/lib/intel64_lin -L$(TBBROOT)/lib/intel64_lin/gcc4.4 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -Wl,-rpath,$(MPI_ROOT)/lib -lpetsc $(MUMPS_LIBS) -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lzoltan -lparmetis -lmetis -ldl -lstdc++ -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lstdc++

PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_LIB = -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           -l$(M3DC1_SCOREC_LIB) $(PUMI_LIB) -Wl,--end-group

LIBS =  $(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas

INCLUDE = -I$(SCOREC_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(GSL_DIR)/include

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
