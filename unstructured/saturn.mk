CPP = mpic++
CC = mpicc
F90 = mpif90
F77 = mpif90
LOADER = mpif90

FOPTS = -DPETSC_VERSION=312 -c -cpp -fdefault-real-8 -fdefault-double-8 -Wall -static $(OPTS) \
        -Wno-unused-variable -ffree-line-length-0 -Wno-unused-dummy-argument

CCOPTS = -DPETSC_VERSION=312 -c -Wall -Wwrite-strings -Wno-strict-aliasing \
         -Wno-unknown-pragmas -g3 $(OPTS)

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

MPIVER=mpich-gcc4.7.2
PETSCVER=petsc3.12.1
PETSC_VER=petsc-3.12.1
PETSC_DIR=/fusion/projects/codes/m3dc1/scorec/$(PETSC_VER)

SCOREC_BASE_DIR=/fusion/projects/codes/m3dc1/scorec/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan
SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB = m3dc1_scorec_complex
  PETSC_ARCH=cplx-$(MPIVER)
else
  M3DC1_SCOREC_LIB = m3dc1_scorec
  PETSC_ARCH=real-$(MPIVER)
endif

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/fusion/usc/opt/fftw/fftw-3.3.4-mpich-gcc-4.7.2/lib -L/fusion/usc/opt/fftw/fftw-3.3.4-mpich-gcc-4.7.2/lib -Wl,-rpath,/fusion/usc/opt/hdf5/hdf5-1.8.16-mpich-gcc-4.7.2/lib -L/fusion/usc/opt/hdf5/hdf5-1.8.16-mpich-gcc-4.7.2/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lquadmath -ldl -lstdc++

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	-I$(FFTW_DIR)/include \
        -I$(HDF5_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include

LIBS := $(LIBS) \
        -L$(SCOREC_DIR)/lib -l$(M3DC1_SCOREC_LIB) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
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
