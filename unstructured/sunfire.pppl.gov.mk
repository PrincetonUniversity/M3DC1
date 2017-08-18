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

PETSC_DIR = /p/swim/jchen/PETSC/petsc-3.5.3/
ifeq ($(COM), 1)
PETSC_ARCH = portalr6-intel-openmpi-1.8.4-complex
HYPRE_LIB=
else
PETSC_ARCH = portalr6-intel-openmpi-1.8.4
HYPRE_LIB=-lHYPRE
endif

SCOREC_UTIL_DIR = /p/tsc/m3dc1/lib/SCORECLib/rhel6/openmpi-1.8.4/utilities

PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
	-ldmumps -lmumps_common -lpord -lcmumps -lsmumps -lzmumps \
	-lfftw3 -lfftw3_mpi \
        $(HYPRE_LIB) \
	-lparmetis -lmetis \
	-lscalapack \
	-lsuperlu_dist_3.3 -lsuperlu_4.3 \
	-Wl,--end-group

BLASLAPACKLIBS = -L$(MKLROOT)/lib/intel64 -Wl,--start-group  \
	-lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_cdft_core -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core \
	-Wl,--end-group
HDF5_DIR=/usr/pppl/intel/2015-pkgs/openmpi-1.8-pkgs/hdf5-1.8.14-parallel
#SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/May2017-openmpi-1.8.4
#PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/Dec2015
ifeq ($(PAR), 1)
  SCOREC_DIR= /p/tsc/m3dc1/lib/SCORECLib/rhel6/Dec2016_PIC
  PUMI_LIB = -lapf -lapf_zoltan -lapf_omega_h -lgmi -llion -lma -lmds -lmth \
    -lomega_h -lparma -lpcu -lph -lsam -lspr -lzoltan
else
  SCOREC_DIR= /p/tsc/m3dc1/lib/SCORECLib/rhel6/Dec2015
  PUMI_LIB = -lapf -lapf_zoltan -lcrv -ldsp -lgmi -lma -lmds -lparma -lpcu \
    -lph -lspr -lzoltan
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCORECLIB= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group
ZOLTAN_LIB=-L/p/tsc/m3dc1/lib/SCORECLib/rhel6/openmpi-1.8.4/lib -lzoltan

LIBS = 	\
	$(SCORECLIB) \
        $(ZOLTAN_LIB) \
        $(BLASLAPACKLIBS) \
        $(PETSC_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 \
	-L$(ZLIB_HOME) -lz \
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
