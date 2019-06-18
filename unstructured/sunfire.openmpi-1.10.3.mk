FOPTS = $(OPTS) -DPETSC_VERSION=37 -c -r8 -implicitnone -fpp -warn all -DLATESTSCOREC -DUSEBLAS
# FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DLATESTSCOREC -DUSEPARTICLES
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
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver
PETSC_VER=petsc-3.7.6
PETSCVER=petsc3.7.6

PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/$(PETSC_VER)

ifeq ($(COM), 1)
PETSC_ARCH=cplx-intel2015-openmpi1.10.3-gcc4.4.7
HYPRE_LIB=
else
PETSC_ARCH=real-intel2015-openmpi1.10.3-gcc4.4.7
HYPRE_LIB=-lHYPRE
endif

SCOREC_BASE_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/intel2015-openmpi1.10.3-gcc4.4.7
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan

ifeq ($(REORDERED), 1)
  PUMI_DIR=$(SCOREC_DIR)/$(PETSCVER)/reordered
else
  PUMI_DIR=$(SCOREC_DIR)/$(PETSCVER)
endif

PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group

COMP_LIB_DIR=/usr/pppl/intel/2015.u1/composer_xe_2015.1.133/compiler/lib/intel64
MPI_LIB_DIR=/usr/pppl/intel/2015-pkgs/openmpi-1.10.3/lib
GCC_HOME=/usr/lib/gcc/x86_64-redhat-linux/4.4.7

PETSC_LIBS =-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(MPI_LIB_DIR) -L$(COMP_LIB_DIR) -L$(GCC_HOME) -Wl,-rpath,$(MPI_LIB_DIR) -lpetsc -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lsuperlu $(HYPRE_LIB) -lscalapack -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhwloc -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lifport -lifcore -lmpi_cxx -lintlc -ldl -lstdc++ -L$(MPI_LIB_DIR) -lmpi -L$(COMP_LIB_DIR) -L$(GCC_HOME) -Wl,-rpath,$(MPI_LIB_DIR) -limf -lsvml -lirng -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lpthread -lirc_s -L$(COMP_LIB_DIR) -L$(GCC_HOME) -ldl -lstdc++

LIBS = 	\
	$(SCOREC_LIB) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIBS) \
        -L$(HDF5_HOME)/lib  -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 \
	-L$(GSL_HOME)/lib -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_HOME)/include \
        -I$(GSL_HOME)/include

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
