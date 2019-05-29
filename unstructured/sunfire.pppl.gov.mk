FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DLATESTSCOREC -DUSEBLAS -DPETSC_VERSION=37
CCOPTS  = -c -DPETSC_VERSION=37

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
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

PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/petsc-3.7.6
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-intel2015-openmpi1.10.3-gcc4.4.7
  HYPRE_LIB=
else
  PETSC_ARCH=real-intel2015-openmpi1.10.3-gcc4.4.7
  HYPRE_LIB=-lHYPRE
endif

BLASLAPACK_LIBS =-Wl,-rpath,$(MKLROOT) -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
HDF5_DIR=$(HDF5_HOME)
SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/Aug2017/openmpi-1.10.3/debug
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_UTIL_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/intel2015-openmpi1.10.3-gcc4.4.7/bin

ZOLTAN_LIB=-L/p/tsc/m3dc1/lib/SCORECLib/rhel6/openmpi-1.10.3/lib -lzoltan
SCALAPACK_LIB=-Wl,-rpath,$(SCALAPACK_HOME)/lib -L$(SCALAPACK_HOME) -lscalapack
PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
        -lsuperlu_dist -lsuperlu \
        -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
        $(SCALAPACK_LIB) \
        -lfftw3 -lfftw3_mpi \
        $(HYPRE_LIB) \
        -lparmetis -lmetis \
        -Wl,--end-group


ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif
SCORECLIB= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group

LIBS = 	\
	$(SCORECLIB) \
        $(TRILINOS_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz \
	-L$(GSL_HOME)/lib -lgsl -lgslcblas \
	-lX11

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_DIR)/include \
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
