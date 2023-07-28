FOPTS = -c -r8 -i4 -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=313 -DDEBUG
R8OPTS = -r8

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

CC = mpicc
CPP = mpic++
F90 = mpifort
F77 = mpifort
LOADER =  mpifort
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

#PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609
PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4
ifeq ($(COM), 1)
  PETSC_ARCH=nvidia-hpc-sdk-6-cuda-cplx
  PETSC_CC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-cplx/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_FC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-cplx/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-cplx/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64/stubs -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lgsl -lgslcblas -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcudart -lnvToolsExt -lcuda -lmpifort_nvidia -lmpi_nvidia -lnvf -lnvomp -ldl -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
else
  ifeq ($(ST), 1)
  PETSC_ARCH=nvidia-hpc-sdk-6-cuda-st
  PETSC_CC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-st/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_FC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-st/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-st/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda-st/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64/stubs -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcudart -lnvToolsExt -lcuda -lmpifort_nvidia -lmpi_nvidia -lnvf -lnvomp -ldl -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
  else
  PETSC_ARCH=nvidia-hpc-sdk-6-cuda
  PETSC_CC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_FC_INCLUDES = -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/include -I/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.real4/nvidia-hpc-sdk-6-cuda/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64/stubs -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -lHYPRE -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lgmp -lgsl -lgslcblas -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcudart -lnvToolsExt -lcuda -lmpifort_nvidia -lmpi_nvidia -lnvf -lnvomp -ldl -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
  endif
endif

#SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/scorec-pmt/nvidia8.3.3-openmpi3.1.5/petsc3.19.1/
SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core/upgrade-nvhpc833-pgpu
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

ifeq ($(ENABLE_ZOLTAN), 0)
  ZOLTAN_LIB=
  BIN_POSTFIX := $(BIN_POSTFIX)-nozoltan
  CCOPTS := $(CCOPTS) -DDISABLE_ZOLTAN
else
  ZOLTAN_LIB=-lzoltan
endif 


LIBS = 	\
	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	/opt/cray/pe/lib64/libmpi_gtl_cuda.so.0
#-L$(HDF5_ROOT)/lib -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lhdf5_hl -lhdf5 \
#-L$(FFTW_DIR) -lfftw3f_threads -lfftw3_threads -lfftw3_mpi -lfftw3f_mpi -lfftw3f -lfftw3

INCLUDE = $(PETSC_FC_INCLUDES) #-I$(FFTW_INC) -I$(HDF5_ROOT)/include \

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L/global/homes/j/jinchen/project/NETCDF/buildnvhpc2/lib -Wl,-rpath,/global/homes/j/jinchen/project/NETCDF/buildnvhpc2/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lnetcdf -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I/global/cfs/cdirs/mp288/jinchen/NETCDF/buildnvhpc2/include
endif

ACC?=1
ifeq ($(ACC), 1)
  LDOPTS := $(LDOPTS) -acc -gpu=cuda11.7 -Minfo=accel
  FOPTS  := $(FOPTS)  -acc -gpu=cuda11.7 -Minfo=accel
  CCOPTS  := $(CCOPTS) -acc -gpu=cuda11.7 -Minfo=accel
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
