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

CC = cc
CPP = CC
F90 = ftn
F77 = ftn
LOADER =  ftn
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714
ifeq ($(COM), 1)
  PETSC_ARCH=perlmuttercpu-nvidia-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia-cplx/lib -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/23.02.1.1/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.5.2-2.4_3.49__gd0f7936.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lgsl -lgslcblas -lsci_nvidia_mpi -lsci_nvidia -ldl -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
else
  ifeq ($(ST), 1)
  PETSC_ARCH=perlmuttercpu-nvidia-st
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia-st/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia-st/lib -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/23.02.1.1/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.5.2-2.4_3.49__gd0f7936.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lsci_nvidia_mpi -lsci_nvidia -ldl -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
  else
  PETSC_ARCH=perlmuttercpu-nvidia
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20230714/perlmuttercpu-nvidia/lib -L/opt/cray/pe/mpich/8.1.25/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/23.02.1.1/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.5.2-2.4_3.49__gd0f7936.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lgsl -lgslcblas -lsci_nvidia_mpi -lsci_nvidia -ldl -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
  endif
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core-trunk/upgrade-nvhpc833-pcpu
#SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/scorec-pmt/nvidia8.3.3-mpich8.1.24/petsc-3.17.4
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
  LIBS += -Wl,--start-group -L/global/homes/j/jinchen/project/NETCDF/buildnvhpc2/lib -Wl,-rpath,/global/homes/j/jinchen/project/NETCDF/buildnvhpc2/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lnetcdf -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I/global/cfs/cdirs/mp288/jinchen/NETCDF/buildnvhpc2/include
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
