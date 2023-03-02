FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DUSEBLAS -DPETSC_VERSION=313
CCOPTS  = -c -DPETSC_VERSION=313

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -g -check-uninit -debug all
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

CC = mpiicc
CPP = mpiicpc
F90 = mpiifort
F77 = mpiifort
LOADER = mpiifort
LDOPTS := $(LDOPTS) -cxxlib
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
MPIVER=intel2021.1.2-intelmpi2021.3.1
PETSC_VER=petsc-3.13.5
PETSCVER=petsc3.13.5
ifeq ($(COM), 1)
PETSC_DIR=/projects/M3DC1/PETSC/$(PETSC_VER)
  PETSC_ARCH=cplx-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex 
SCOREC_BASE_DIR=/projects/M3DC1/scorec/$(MPIVER)/$(PETSCVER)/202209
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib \
	-Wl,-rpath,/projects/M3DC1/PETSC/petsc-3.13.5/${PETSC_ARCH}/lib \
	-L/projects/M3DC1/PETSC/petsc-3.13.5/${PETSC_ARCH}/lib  \
	-lpetsc \
	-lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common \
	-lpord -lsuperlu -lsuperlu_dist \
	-lscalapack -lflapack -lfblas -lzoltan -lX11 \
	-lparmetis -lmetis \
	-lz -lmpifort -lmpi -lrt -lpthread \
	-lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lquadmath \
	-lstdc++ -ldl
else
PETSC_DIR=/projects/M3DC1/PETSC/petsc.20210709
  PETSC_ARCH=real-stellar-intelmpi-2021
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
#SCOREC_BASE_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
SCOREC_BASE_DIR=/projects/M3DC1/jinchen/SRC/M3DC1/m3dc1_scorec/build
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/projects/M3DC1/PETSC/petsc.20210709/real-stellar-intelmpi-2021/lib -L/projects/M3DC1/PETSC/petsc.20210709/real-stellar-intelmpi-2021/lib -Wl,-rpath,/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/mpi/2021.3.1/lib/release -L/opt/intel/oneapi/mpi/2021.3.1/lib -L/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/intel-mpi/4.7.4/lib64 -L/usr/local/gsl/2.6/x86_64/lib64 -L/usr/local/hdf5/intel-2021.1/intel-mpi/1.10.6/lib64 -L/usr/local/fftw/intel-2021.1/intel-mpi/3.3.9/lib64 -L/opt/intel/oneapi/mpi/2021.3.1/libfabric/lib -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/opt/intel/oneapi/mpi/2021.3.1/lib/release -Wl,-rpath,/opt/intel/oneapi/mpi/2021.3.1/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lz -lstdc++ -lmpifort -lmpi -ldl -lrt -lpthread -lifport -lifcoremt -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lquadmath -lstdc++
endif

SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

LIBS = 	-L$(I_MPI_ROOT)/lib -lmpicxx\
	$(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW3DIR)/lib -lfftw3 -lfftw3f_mpi -lfftw3l_mpi -lfftw3_mpi \
	-L$(HDF5DIR)/lib64 -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran \
	-L$(GSL_ROOT_DIR)/lib64 -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_DIR)/include \
        -I$(HDF5DIR)/include \
        -I$(GSL_ROOT_DIR)/include

ifeq ($(ST), 1)
  LIBS += -L$(NETCDFDIR)/lib64 -lnetcdf -lnetcdff

  INCLUDE += -I$(NETCDFDIR)/include 
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
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
