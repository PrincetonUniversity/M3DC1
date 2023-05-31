FOPTS = -c -fdefault-real-8 -Wall -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=313 -DDEBUG

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

CC = $(MPICH_ROOT)/bin/mpicc
CPP = $(MPICH_ROOT)/bin/mpicxx
F90 = $(MPICH_ROOT)/bin/mpif90
F77 = $(MPICH_ROOT)/bin/mpif90
LOADER = $(MPICH_ROOT)/bin/mpif90
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/lore/seol/petsc-3.13.5
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-gcc4.8.5-v5m6xwi-mpich3.2.1-geowaxe
else
  PETSC_ARCH=real-gcc4.8.5-v5m6xwi-mpich3.2.1-geowaxe
endif

SCOREC_BASE_DIR=/lore/seol/rhel7-gcc4.8.5-v5m6xwi-mpich3.2.1-geowaxe
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

ifeq ($(ENABLE_ZOLTAN), 0)
  ZOLTAN_LIB=
  BIN_POSTFIX := $(BIN_POSTFIX)-nozoltan
  CCOPTS := $(CCOPTS) -DDISABLE_ZOLTAN
else
  ZOLTAN_LIB=-lzoltan
endif 

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-4.8.5/mpich-3.2.1-geowaxeufjfz7vgpszmfrm4l6rp4rnxo/lib -Wl,-rpath,/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-4.8.5-v5m6xwipewvkufflu3zvxs2x6jv4libr/lib:/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-4.8.5-v5m6xwipewvkufflu3zvxs2x6jv4libr/lib64 -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-4.8.5-v5m6xwipewvkufflu3zvxs2x6jv4libr/lib64 -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-4.8.5-v5m6xwipewvkufflu3zvxs2x6jv4libr/lib/gcc/x86_64-unknown-linux-gnu/4.8.5 -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-4.8.5/hdf5-1.10.4-wousvnks7b7aplenali7jg5bsh243hyd/lib -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-4.8.5/zlib-1.2.11-vhzh5cfaki5lx5sjuth5iuojq5azdkbd/lib -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-4.8.5-v5m6xwipewvkufflu3zvxs2x6jv4libr/lib -Wl,-rpath,/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-4.8.5/mpich-3.2.1-geowaxeufjfz7vgpszmfrm4l6rp4rnxo/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas $(ZOLTAN_LIB) -lpthread -lX11 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lm -lz -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl

LIBS = 	\
	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(GSL_ROOT)/lib -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
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
