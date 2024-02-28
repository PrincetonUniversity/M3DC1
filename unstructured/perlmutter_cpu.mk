FOPTS = -c -r8 -i4 -cpp -DPETSC_VERSION=319 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=319

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

CC = cc
CPP = CC
F90 = ftn
F77 = ftn
LOADER =  ftn
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20240207
ifeq ($(COM), 1)
  PETSC_ARCH=perlmuttercpu-nvidia-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20240207/perlmuttercpu-nvidia-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20240207/perlmuttercpu-nvidia-cplx/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lstdc++ -lquadmath
else
  PETSC_ARCH=perlmuttercpu-nvidia
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20240207/perlmuttercpu-nvidia/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20240207/perlmuttercpu-nvidia/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lstdc++ -lquadmath
endif

#SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/scorec-pmt/cpu-nvidia8.3.3-mpich8.1.25/petsc3.19.3
SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core-trunk/upgrade-nvhpc850-pcpu
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

LIBS = 	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L$(HDF5_DIR)/lib -L$(NETCDF_DIR)/lib \
	  -Wl,-rpath,$(NETCDF_DIR)/lib -Wl,-rpath,$(HDF5_DIR)/lib \
	  -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran \
	  -lnetcdf -lnetcdff \
	  -lz -Wl,--end-group
  INCLUDE += -I$(HDF5_DIR)/include -I$(NETCDF_DIR)/include
endif


%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP) $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
