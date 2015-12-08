CPP = mpic++
CC = mpicc
F90 = mpif90
F77 = mpif90
LOADER = mpif90

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS
FOPTS = -c -cpp -fdefault-real-8 -Wall -static $(OPTS) \
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
SCOREC_DIR = /global/project/projectdirs/mp288/cori/scorec/Dec2015
SCOREC_CORE = -lapf -lgmi -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr

FFTW_HOME = /fusion/usc/opt/fftw/fftw-3.3.4/mpich-gcc-4.7.2

ifeq ($(COM), 1)
  SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec_complex -Wl,--end-group
  PETSC_DIR = /fusion/usc/opt/petsc/petsc-3.5.4
#  PETSC_ARCH = linux-mpich-gcc-4.7.2
  PETSC_ARCH = linux-mpich-gcc-4.7.2-hdf5-netcdf-hypre
else
      SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec -Wl,--end-group
  PETSC_DIR = /fusion/usc/opt/petsc/petsc-3.5.4
#  PETSC_ARCH = linux-mpich-gcc-4.7.2
  PETSC_ARCH = linux-mpich-gcc-4.7.2-hdf5-netcdf-hypre
endif

PETSC_EXTERNAL_LIB_BASIC = \
       -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib \
       -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
       -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord \
       -lsuperlu_4.3 -lsuperlu_dist_3.3 \
       -lparmetis -lmetis -lm -lpthread -lssl -lcrypto \
       -ldl

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	-I$(FFTW_HOME)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \

LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc $(PETSC_EXTERNAL_LIB_BASIC) \
        -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(FFTW_HOME)/lib -lfftw3 \
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
