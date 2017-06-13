FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DPETSC_31 -DKSPITS $(OPTS) -DLATESTSCOREC -DUSEBLAS
# FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DPETSC_31 -DKSPITS $(OPTS) -DLATESTSCOREC -DUSEPARTICLES
CCOPTS  = -c -O -DPetscDEV -DPETSC_31 -DPetscOLD #-DCJ_MATRIX_DUMP -DUSEHYBRID 

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
  LOADER = mpif90 -cxxlib
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/petsc-3.5.4
ifeq ($(COM), 1)
PETSC_ARCH=complex-openmpi-1.10.3
HYPRE_LIB=
else
PETSC_ARCH=real-openmpi-1.10.3
HYPRE_LIB=-lHYPRE
endif

BLASLAPACK_LIBS =-Wl,-rpath,$(MKLROOT) -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
HDF5_DIR=/usr/pppl/intel/2015-pkgs/openmpi-1.10.3-pkgs/hdf5-parallel-1.8.17
SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/June2017/openmpi-1.10.3/debug
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifeq ($(TRILINOS), 1)
  TRILINOS_DIR=/usr/pppl/intel/2015-pkgs/openmpi-1.10.3-pkgs/trilinos-11.12.1
  ZOLTAN_LIB=-L$(TRILINOS_DIR)/lib -lzoltan
  TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_DIR)/lib -L$(TRILINOS_DIR)/lib \
        -lstdc++  -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra -lepetraext \
        -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics -lteuchosremainder
  PETSC_LIBS=
else
  ZOLTAN_LIB=-L/p/tsc/m3dc1/lib/SCORECLib/rhel6/openmpi-1.10.3/lib -lzoltan
  TRILINOS_LIBS=
  SCALAPACK_LIB=-Wl,-rpath,$(SCALAPACK_HOME)/lib -L$(SCALAPACK_HOME) -lscalapack
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
        -lsuperlu_dist_3.3 -lsuperlu_4.3 \
        -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
        $(SCALAPACK_LIB) \
        -lfftw3 -lfftw3_mpi \
        $(HYPRE_LIB) \
        -lparmetis -lmetis \
        -Wl,--end-group
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB=-lm3dc1_scorec
  endif
endif

SCORECLIB= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group

LIBS = 	\
	$(SCORECLIB) \
        $(TRILINOS_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 \
	-L$(ZLIB_DIR) -lz \
	-L$(GSLHOME)/lib -lgsl -lgslcblas \
	-lX11

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_DIR)/include \
        -I$(GSLHOME)/include

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
