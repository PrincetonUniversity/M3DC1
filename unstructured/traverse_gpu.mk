  CPP = mpic++
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort

OPTS := $(OPTS) -DPETSC_VERSION=990 -DUSEBLAS

PETSC_DIR=/projects/M3DC1/petsc

ifeq ($(COM), 1)
   PETSC_ARCH=traverse-nvidia-complex
else
   PETSC_ARCH=traverse-nvidia-real
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_BASE_DIR=/projects/M3DC1/scorec/traverse-nvidia
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L/home/liuchang/zoltan/pgi/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

HDF5_DIR=/projects/M3DC1/hdf5/traverse-nvidia
NETCDFDIR=/projects/M3DC1/netcdf/traverse-nvidia
GSL_DIR=/projects/M3DC1/gsl/traverse
FFTW_DIR=/projects/M3DC1/fftw/traverse-nvidia

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lsuperlu_dist -lsuperlu -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -llapack -lblas -lstdc++

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
       -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
       -I$(NETCDFDIR)/include \
       -I$(HDF5_DIR)/include \
       -I$(FFTW_DIR)/include \
       -I$(GSL_DIR)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        -L$(NETCDFDIR)/lib -lnetcdff -lnetcdf -lzip -lcurl\
        -L$(HDF5_DIR)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz\
        -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
        -L$(GSL_DIR)/lib -lgsl -lgslcblas \
        $(PETSC_WITH_EXTERNAL_LIB)

FOPTS = -c -r8 -Mpreprocess $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -fast
  FOPTS  := $(FOPTS)  -fast
  CCOPTS := $(CCOPTS) -fast
else
  FOPTS := $(FOPTS) -Mbounds -Minfo=all -Mchkfpstk -Mchkstk -Mdalign -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
  CCOPTS := $(CCOPTS) -g
  LDOPTS := $(LDOPTS) -g
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

ifeq ($(ACC), 1)
  LDOPTS := $(LDOPTS) -acc -gpu=cuda11.3 -Minfo=accel
  FOPTS  := $(FOPTS)  -acc -gpu=cuda11.3 -Minfo=accel
  CCOPTS  := $(CCOPTS) -acc -gpu=cuda11.3 -Minfo=accel
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

