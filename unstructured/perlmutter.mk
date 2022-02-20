FOPTS = -c -r8 -i4 -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
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

CC = cc
CPP = CC
F90 = ftn
F77 = ftn
LOADER =  ftn
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927
ifeq ($(COM), 1)
  PETSC_ARCH=nvidia-hpc-sdk-6-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6-cplx/lib -L/opt/cray/pe/mpich/8.1.11/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/21.08.1.2/NVIDIA/20.7/x86_64/lib -L/global/common/software/nersc/pm-2021q4/sw/darshan/3.3.1/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.2.40-2.1_3.9__g3cf3325.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lgsl -lgslcblas -lstdc++ -ldl -ldarshan -lz -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath-avx2 -lnsnvc -lnvc -lrt -lm -lgcc_s -lquadmath -lstdc++ -ldl
else
  ifeq ($(ST), 1)
  PETSC_ARCH=nvidia-hpc-sdk-6-st
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6-st/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6-st/lib -L/opt/cray/pe/mpich/8.1.12/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/21.08.1.2/NVIDIA/20.7/x86_64/lib -L/global/common/software/nersc/pm-2021q4/sw/darshan/3.3.1/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.2.40-2.1_3.9__g3cf3325.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lstdc++ -ldl -ldarshan -lz -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath-avx2 -lnsnvc -lnvc -lrt -lm -lgcc_s -lquadmath -lstdc++ -ldl
  else
  PETSC_ARCH=nvidia-hpc-sdk-6
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20210927/nvidia-hpc-sdk-6/lib -L/opt/cray/pe/mpich/8.1.11/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/21.08.1.2/NVIDIA/20.7/x86_64/lib -L/global/common/software/nersc/pm-2021q4/sw/darshan/3.3.1/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.2.40-2.1_3.9__g3cf3325.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lgsl -lgslcblas -lstdc++ -ldl -ldarshan -lz -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath-avx2 -lnsnvc -lnvc -lrt -lm -lgcc_s -lquadmath -lstdc++ -ldl
  endif
endif

SCOREC_BASE_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
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
#-L$(HDF5_ROOT)/lib -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lhdf5_hl -lhdf5 \
#-L$(FFTW_DIR) -lfftw3f_threads -lfftw3_threads -lfftw3_mpi -lfftw3f_mpi -lfftw3f -lfftw3

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
#-I$(FFTW_INC) -I$(HDF5_ROOT)/include \


ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L/global/homes/j/jinchen/project/NETCDF/buildnvhpc/lib -Wl,-rpath,/global/homes/j/jinchen/project/NETCDF/buildnvhpc/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lnetcdf -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I/global/cfs/cdirs/mp288/jinchen/NETCDF/buildnvhpc/include
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
