FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DPETSC_31 -DKSPITS $(OPTS) -DLATESTSCOREC
CCOPTS  = -c -O -DPetscDEV -DPETSC_31 -DPetscOLD #-DCJ_MATRIX_DUMP -DUSEHYBRID 

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
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


# define where you want to locate the mesh adapt libraries
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

PETSC_DIR = /p/swim/jchen/PETSC/petsc-3.5.3/
INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_HOME)/include -I$(HDF5_HOME)/lib \
	-I$(GSLHOME)/include

#	-I$(SUPERLU_DIST_HOME)/include -I$(SUPERLU_HOME)/include \
#	-I$(HYBRID_HOME)/include \
#
ifeq ($(COM), 1)
PETSC_ARCH = portalr6-intel-openmpi-1.8.4-complex
PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
	-ldmumps -lmumps_common -lpord -lcmumps -lsmumps -lzmumps \
	-lfftw3 -lfftw3_mpi \
	-lparmetis -lmetis \
	-lscalapack \
	-lsuperlu_dist_3.3 -lsuperlu_4.3 \
	-Wl,--end-group
else
PETSC_ARCH = portalr6-intel-openmpi-1.8.4
PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
	-ldmumps -lmumps_common -lpord -lcmumps -lsmumps -lzmumps \
	-lfftw3 -lfftw3_mpi \
	-lHYPRE \
	-lparmetis -lmetis \
	-lscalapack \
	-lsuperlu_dist_3.3 -lsuperlu_4.3 \
	-Wl,--end-group
endif


#	-lfblas -lflapack \
#	-lpromfei -lprometheus \
#	-lHYPRE \

#SUPERLU_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
#SUPERLU_DIST_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
#SUPERLU_LIBS = -L$(SUPERLU_HOME)/lib -lsuperlu_4.3 \
#	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_3.3 \

#	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacsCinit -lmpiblacs

#PARMETIS_LIBS = -L$(PARMETIS_HOME)/lib \
#	-Wl,-rpath,$(PARMETIS_HOME)/lib -lparmetis -lmetis

BLASLAPACKLIBS = -L$(MKLROOT)/lib/intel64 -Wl,--start-group \
	-lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_cdft_core -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core \
	-Wl,--end-group

SCORECDIR= /p/tsc/m3dc1/lib/SCORECLib/rhel6/Dec2015
ifeq ($(COM), 1)
  SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib \
             -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lzoltan -lm3dc1_scorec_complex \
             -Wl,--end-group
else
  SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib \
             -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lzoltan -lm3dc1_scorec \
             -Wl,--end-group
endif

LIBS = 	\
	$(SCORECLIB) \
        $(BLASLAPACKLIBS) \
        $(PETSC_LIBS) \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
	-L$(ZLIB_HOME) -lz \
	-L$(GSLHOME)/lib -lgsl -lgslcblas \
	-L/usr/lib -lX11

#	$(SUPERLU_LIBS) \
#	-L$(Zoltan_HOME)/lib -lzoltan \
#	$(PARMETIS_LIBS) \
#	-L$(FFTWHOME)/lib -lfftw3 \
#	-L$(ACML_HOME)/ifort64/lib -lacml \
#	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
  INCLUDE := -I$(SCORECDIR)/include  $(INCLUDE)

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
