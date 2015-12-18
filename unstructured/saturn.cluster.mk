CPP = mpic++
CC = mpicc
F90 = mpif90
F77 = mpif90
LOADER = mpif90

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS
FOPTS = -c -cpp -fdefault-real-8 -fdefault-double-8 -Wall -static $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS) -O3
  CCOPTS := $(CCOPTS) -O3
else
  FOPTS := $(FOPTS) -g 
  CCOPTS := $(CCOPTS)  
endif


F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)


# define where you want to locate the mesh adapt libraries
#HYBRID_HOME =  /scratch2/scratchdirs/xyuan/Software_Hopper/pdslin_0.0
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lpdslin
PETSC_DIR = /fusion/usc/opt/petsc
SCOREC_DIR = /fusion/projects/codes/m3dc1/scorec/Dec2015
SCOREC_CORE = -lapf -lgmi -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr -lzoltan

ifeq ($(COM), 1)
  SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec_complex -Wl,--end-group
  PETSC_ARCH = petsc-3.5.4-mpich-gcc-4.7.2-hdf5-netcdf-complex-nohypre
else
  SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec -Wl,--end-group
  PETSC_ARCH = petsc-3.5.4-mpich-gcc-4.7.2-hdf5-netcdf-hypre-nocomplex
endif

METIS_DIR = /fusion/usc/opt/metis/metis-5.1.0-mpich-gcc-4.7.2

PETSC_EXTERNAL_LIB_BASIC = \
       -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
       -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
       -lm -lpthread -lssl -lcrypto \
       -ldl

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	-I$(FFTW_DIR)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	-I$(SUPERLU_DIR)/include -I$(HDF5_PREFIX)/include

LIBS := $(LIBS) $(SCOREC_LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc $(PETSC_EXTERNAL_LIB_BASIC) \
	-L$(HYPRE)/lib -lHYPRE \
        -L$(HDF5_PREFIX)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
	-lgsl -lgslcblas -lhugetlbfs \
	-L$(SUPERLU_DIR)/lib -lsuperlu_4.3 -lsuperlu_dist_3.3 \
	-L$(MUMPS_DIR)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
	-L$(PARMETIS_DIR)/lib -lparmetis \
	-L$(METIS_DIR)/lib -lmetis \
	-L$(SCALAPACK_DIR)/lib -lscalapack \
	-llapack -lblas -lstdc++

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
