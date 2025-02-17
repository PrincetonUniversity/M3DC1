FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=319 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=319
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -w -fallow-argument-mismatch
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g #not for gcc : noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

CC = mpicc
CPP = mpiCC
F90 = mpif90
F77 = mpif90
LOADER =  mpif90
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/prj/ntm/jin.chen/LIB/petsc
ifeq ($(COM), 1)
  PETSC_ARCH=sdumontgnu-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64 -Wl,-rpath,/scratch/app/openmpi/4.1.4_gnu+gcc-12.4/lib -L/scratch/app/openmpi/4.1.4_gnu+gcc-12.4lib -Wl,-rpath,/scratch/app/ucx/1.13+gcc-12.4/lib -L/scratch/app/ucx/1.13+gcc-12.4/lib -Wl,-rpath,/scratch/app/xpmem/2.6.5/lib -L/scratch/app/xpmem/2.6.5/lib -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib/gcc/x86_64-pc-linux-gnu/13.2.0 -L/petrobr/app_sequana/gcc/13.2/lib/gcc/x86_64-pc-linux-gnu/13.2.0 -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib/gcc -L/petrobr/app_sequana/gcc/13.2/lib/gcc -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib64 -L/petrobr/app_sequana/gcc/13.2/lib64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/clck/2017.2.019/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/clck/2017.2.019/lib/intel64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.7 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.7 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib -L/petrobr/app_sequana/gcc/13.2/lib -lpetsc -lsuperlu_dist -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lm -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
else
  PETSC_ARCH=sdumontgnu
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64 -Wl,-rpath,/scratch/app/openmpi/4.1.4_gnu+gcc-12.4/lib -L/scratch/app/openmpi/4.1.4_gnu+gcc-12.4/lib -Wl,-rpath,/scratch/app/ucx/1.13+gcc-12.4/lib -L/scratch/app/ucx/1.13+gcc-12.4/lib -Wl,-rpath,/scratch/app/xpmem/2.6.5/lib -L/scratch/app/xpmem/2.6.5/lib -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib/gcc/x86_64-pc-linux-gnu/13.2.0 -L/petrobr/app_sequana/gcc/13.2/lib/gcc/x86_64-pc-linux-gnu/13.2.0 -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib/gcc -L/petrobr/app_sequana/gcc/13.2/lib/gcc -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib64 -L/petrobr/app_sequana/gcc/13.2/lib64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/clck/2017.2.019/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/clck/2017.2.019/lib/intel64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.7 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.7 -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -Wl,-rpath,/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -L/opt/intel/parallel_studio_xe_2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -Wl,-rpath,/petrobr/app_sequana/gcc/13.2/lib -L/petrobr/app_sequana/gcc/13.2/lib -lpetsc -lsuperlu_dist -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lm -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
endif

SCOREC_BASE_DIR=/prj/ntm/jin.chen/LIB/core-240226/sdumontgnu
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

FFTW_ROOT=/scratch/app_sequana/fftw/3.3.10_openmpi-4.1.4_gnu
HDF5_DIR=/scratch/app_sequana/hdf5/1.14.5_openmpi-4.1.4
GSL_DIR=/scratch/app_sequana/gsl/2.8_gcc-9.3_gnu
LIBS = 	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW_ROOT)/lib -lfftw3 -lfftw3_mpi \
	-L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran \
	-L$(GSL_DIR)/lib -lgsl \
	-lmpi_cxx 

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include \
	-I$(FFTW_ROOT)/include \
	-I$(HDF5_DIR)/include \
	-I$(GSL_DIR)/include

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L$(HDF5_DIR)/lib -L$(NETCDF_DIR)/lib \
	  -Wl,-rpath,$(NETCDF_DIR)/lib -Wl,-rpath,$(HDF5_DIR)/lib \
	  -lnetcdf -lnetcdff \
	  -lz -Wl,--end-group
  INCLUDE += -I$(NETCDF_DIR)/include
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
