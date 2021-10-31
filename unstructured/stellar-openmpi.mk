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

CC = /usr/local/openmpi/4.1.0/intel20211/bin/mpicc
CPP = /usr/local/openmpi/4.1.0/intel20211/bin/mpicxx
F90 =/usr/local/openmpi/4.1.0/intel20211/bin/mpif90
F77 = /usr/local/openmpi/4.1.0/intel20211/bin/mpif90
LOADER =/usr/local/openmpi/4.1.0/intel20211/bin/mpif90
LDOPTS := $(LDOPTS) -cxxlib
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
MPIVER=intel2021.1.2-openmpi4.1.0
PETSC_VER=petsc-3.13.5
PETSCVER=petsc3.13.5
PETSC_DIR=/projects/M3DC1/PETSC/$(PETSC_VER)
ifeq ($(COM), 1)
  PETSC_ARCH=complex-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex 
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/projects/M3DC1/PETSC/petsc-3.13.5/complex-intel2021.1.2-openmpi4.1.0/lib -L/projects/M3DC1/PETSC/petsc-3.13.5/complex-intel2021.1.2-openmpi4.1.0/lib -Wl,-rpath,/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib -Wl,-rpath,/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/openmpi/4.1.0/intel20211/lib64 -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/4.1.0/intel20211/lib64 -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lX11 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl
else
  PETSC_ARCH=real-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib -Wl,-rpath,/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/openmpi/4.1.0/intel20211/lib64 -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/4.1.0/intel20211/lib64 -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lX11 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl
endif

SCOREC_BASE_DIR=/projects/M3DC1/scorec/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

LIBS = 	-L$(OPENMPI_HOME)/lib64 -lmpi_cxx \
	$(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW3DIR)/lib -lfftw3 -lfftw3f_mpi -lfftw3l_mpi -lfftw3_mpi \
	-L$(HDF5DIR)/lib64 -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran \
	-L$(GSL_ROOT_DIR)/lib64 -lgsl -lgslcblas 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
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
