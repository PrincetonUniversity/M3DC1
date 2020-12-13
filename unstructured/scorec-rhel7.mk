FOPTS = -c -fdefault-real-8 -Wall -cpp -DPETSC_VERSION=38 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=38

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

MPICH_ROOT=/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.2.1-niuhmadq72otoeruzmhwp6f7b3rqe24y

CC = $(MPICH_ROOT)/bin/mpicc
CPP = $(MPICH_ROOT)/bin/mpicxx
F90 = $(MPICH_ROOT)/bin/mpif90
F77 = $(MPICH_ROOT)/bin/mpif90
LOADER = $(MPICH_ROOT)/bin/mpif90
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR = /lore/seol/petsc-3.8.2/
ifeq ($(COM), 1)
  PETSC_ARCH =complex-mpich3.2.1-niuhmad
else
  PETSC_ARCH =real-mpich3.2.1-niuhmad
endif

SCOREC_DIR = /lore/seol/mpich3.2.1-gcc7.3.0-install
HDF5_DIR =  $(PETSC_DIR)/$(PETSC_ARCH)
ZLIB_DIR = $(PETSC_DIR)/$(PETSC_ARCH)
FFTW_DIR = $(PETSC_DIR)/$(PETSC_ARCH)
ZOLTAN_LIB=-L$(SCOREC_DIR)/lib -lzoltan

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

PUMI_LIBS = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
SCOREC_LIBS = -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           $(PUMI_LIBS) $(M3DC1_SCOREC_LIB) -Wl,--end-group

PETSC_LIBS = -L/lore/seol/petsc-3.8.2/real-mpich3.2.1-niuhmad/lib -Wl,-rpath,/lore/seol/petsc-3.8.2/real-mpich3.2.1-niuhmad/lib -L$GCC_ROOT/lib64 -L$MPICH_ROOT/lib -L$GCC_ROOT/lib/gcc/x86_64-pc-linux-gnu/7.3.0 -L$GCC_ROOT/lib -Wl,-rpath,$GCC_ROOT/lib -Wl,-rpath,$MPICH_ROOT/lib -lpetsc -lsuperlu -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lflapack -lfblas -lparmetis -lmetis -lX11 -lpthread -lm -lmpifort -lgfortran -lquadmath -lmpicxx -lstdc++ -Wl,-rpath,$GCC_ROOT/lib -L$MPICH_ROOT/lib -L$GCC_ROOT/lib64 -L$GCC_ROOT/lib64 -L$GCC_ROOT/lib/gcc/x86_64-pc-linux-gnu/7.3.0 -L$GCC_ROOT/lib64 -L$MPICH_ROOT/lib -L$GCC_ROOT/lib64 -L$GCC_ROOT/lib -L$GCC_ROOT/lib -ldl -Wl,-rpath,$MPICH_ROOT/lib -lmpi -lgcc_s -ldl


LIBS = 	\
	$(SCOREC_LIBS) \
    $(ZOLTAN_LIB) \
    $(PETSC_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 \
    -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
	-L$(ZLIB_DIR) -lz \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(HDF5_HOME)/include \
        -I$(GSL_ROOT)/include

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
