#FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)
FOPTS = -c -r8 -implicitnone -fpp -warn all -DKSPITS -DxCJ_MATRIX_DUMP  $(OPTS) #-pg
CCOPTS  = -c -O -DKSPITS -DxCJ_MATRIX_DUMP -DxUSEHYBRID # -pg

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -vec-report=0 #-fast
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
  LOADER = mpif90 -cxxlib #-pg
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)


SCOREC_HOME=/p/tsc/m3dc1/lib/SCORECLib/rhel5/Mar2016
ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_HOME)/lib -L$(SCOREC_HOME)/lib \
    -l$(M3DC1_SCOREC_LIB) \
    -lcrv -ldsp -lph -lsize -lsam -lspr -lma -lapf_zoltan -lparma -lmds -lapf -llion -lmth -lgmi -lpcu \
    -Wl,--end-group

# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/pdslin_0.0
HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

INCLUDE = -I$(MPIHOME)/include \
	  -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
          -I$(SCOREC_HOME)/include \
	  -I$(PARMETIS_HOME)/include \
	  -I$(FFTWHOME)/include \
	  -I$(SUPERLU_DIST_HOME)/include -I$(SUPERLUHOME)/include \
	  -I$(BLACS_HOME)/include \
	  -I$(Zoltan_HOME)/include \
	  -I$(GSLHOME)/include \
	  -I$(HDF5_HOME)/include \
	  -I$(NCARG_ROOT)/include

LIBS = 	-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -ldmumps -lmumps_common -lcmumps -lzmumps -lsmumps -lpord \
        -L$(SCOREC_LIBS) \
        -L$(ZOLTAN_HOME)/lib -lzoltan \
        -L$(PARMETIS_HOME)/lib -lparmetis -lmetis \
	-L$(FFTWHOME)/lib -lfftw3 -lfftw3_mpi -lfftw3_threads \
	-L$(SUPERLU_DIST_HOME)/lib -lsuperlu_dist_3.3 -L$(SUPERLUHOME)/lib -lsuperlu_4.3 \
	-L$(SCALAPACK_HOME)/lib -lscalapack -L$(BLACS_HOME)/lib -lmpiblacs -lmpiblacsCinit -lmpiblacsF77init \
	-L$(GSLHOME)/lib -lgsl \
        -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
        -L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L/usr/X11R6/lib -lX11

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
