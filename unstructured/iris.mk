CPP = mpic++
CC = mpicc
F90 = mpif90
F77 = mpif90
LOADER = mpif90

FOPTS = -c -cpp -fdefault-real-8 -fdefault-double-8 -Wall -static $(OPTS) \
        -Wno-unused-variable -ffree-line-length-0 -Wno-unused-dummy-argument \
        -Dglobalinsertval=insertval -Dglobalentdofs=entdofs

CCOPTS  = -c  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -g3 $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS) -O3
  CCOPTS := $(CCOPTS) -O3
else 
  FOPTS := $(FOPTS) -g -O0 
  CCOPTS := $(CCOPTS) -O0 
endif

 
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

SCOREC_UTIL_DIR=/fusion/projects/codes/m3dc1/scorec/tools

# define where you want to locate the mesh adapt libraries
#HYBRID_HOME =  /scratch2/scratchdirs/xyuan/Software_Hopper/pdslin_0.0
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lpdslin
MPI_DIR = /act/mpich/gcc-4.7.2
#PETSC_DIR = /fusion/usc/opt/petsc/petsc-3.5.4
PETSC_DIR=/fusion/projects/codes/m3dc1/petsc-3.5.4
NETCDF_DIR=/fusion/projects/codes/m3dc1/petsc-3.5.4
SCOREC_DIR = /fusion/projects/codes/m3dc1/scorec/Dec2015-petsc
SCOREC_CORE = -lapf -lgmi -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr

ifeq ($(COM), 1)
  SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec_complex -Wl,--end-group
  PETSC_DIR=/fusion/projects/codes/m3dc1/petsc-3.5.4-complex
else
  SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec -Wl,--end-group
  PETSC_DIR=/fusion/projects/codes/m3dc1/petsc-3.5.4-real
  HYPRE_LIB=-lHYPRE
endif
PETSC_ARCH=mpich-gcc4.7.2

PETSC_EXTERNAL_LIB_BASIC = \
       -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
       -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
-lsuperlu_dist_3.3 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_4.3 $(HYPRE_LIB) -L/act/mpich/gcc-4.7.2/lib -L/act/gcc-4.7.2/lib/gcc/x86_64-unknown-linux-gnu/4.7.2 -L/act/gcc-4.7.2/lib64 -L/act/gcc-4.7.2/lib -Wl,-rpath,/act/mpich/gcc-4.7.2/lib -lmpichcxx -lstdc++ -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lmpichf90 -lgfortran -lm -lgfortran -lm -lquadmath -lm -lmpichcxx -lstdc++ -L/act/mpich/gcc-4.7.2/lib -L/act/gcc-4.7.2/lib/gcc/x86_64-unknown-linux-gnu/4.7.2 -L/act/gcc-4.7.2/lib64 -L/act/gcc-4.7.2/lib -ldl -Wl,-rpath,/act/mpich/gcc-4.7.2/lib -lmpich -lopa -lmpl -lrt -lpthread -lgcc_s -ldl

NETCDF_LIB= -Wl,-rpath,$(NETCDF_DIR)/lib -L$(NETCDF_DIR)/lib -lnetcdf 

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	-I$(FFTW_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        -L$(SCOREC_DIR)/lib -lzoltan \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc $(PETSC_EXTERNAL_LIB_BASIC) \
        $(NETCDF_LIB) \
	-L$(FFTW_DIR)/lib -lfftw3 \
	-lgsl -lgslcblas -lhugetlbfs \
	-lstdc++

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
