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

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20250226
ifeq ($(COM), 1)
  PETSC_ARCH=perlmuttercpu-nvidia-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.30/ofi/nvidia/23.3/lib -L/opt/cray/pe/mpich/8.1.30/ofi/nvidia/23.3/lib -Wl,-rpath,/opt/cray/pe/libsci/24.07.0/NVIDIA/23.3/x86_64/lib -L/opt/cray/pe/libsci/24.07.0/NVIDIA/23.3/x86_64/lib -Wl,-rpath,/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/nvidia/23.3/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/nvidia/23.3/lib -Wl,-rpath,/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7/lib -L/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7/lib -Wl,-rpath,/opt/cray/pe/dsmml/0.3.0/dsmml/lib -L/opt/cray/pe/dsmml/0.3.0/dsmml/lib -Wl,-rpath,/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -L/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -Wl,-rpath,/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -L/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/7 -L/usr/lib64/gcc/x86_64-suse-linux/7 -lpetsc -lsuperlu_dist -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lhdf5_hl_parallel -lhdf5_parallel -ldarshan -llustreapi -lz -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lnetcdf -lnetcdff -lsci_nvidia_mpi -lsci_nvidia -ldl -lmpifort_nvidia -ldsmml -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lfftw3_mpi -lmpi_nvidia -lfftw3_threads -lfftw3 -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
else
  PETSC_ARCH=perlmuttercpu-nvidia
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.30/ofi/nvidia/23.3/lib -L/opt/cray/pe/mpich/8.1.30/ofi/nvidia/23.3/lib -Wl,-rpath,/opt/cray/pe/libsci/24.07.0/NVIDIA/23.3/x86_64/lib -L/opt/cray/pe/libsci/24.07.0/NVIDIA/23.3/x86_64/lib -Wl,-rpath,/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/nvidia/23.3/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/nvidia/23.3/lib -Wl,-rpath,/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7/lib -L/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7/lib -Wl,-rpath,/opt/cray/pe/dsmml/0.3.0/dsmml/lib -L/opt/cray/pe/dsmml/0.3.0/dsmml/lib -Wl,-rpath,/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -L/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -Wl,-rpath,/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -L/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/7 -L/usr/lib64/gcc/x86_64-suse-linux/7 -lpetsc -lsuperlu_dist -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lhdf5_hl_parallel -lhdf5_parallel -ldarshan -llustreapi -lz -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lnetcdf -lnetcdff -lsci_nvidia_mpi -lsci_nvidia -ldl -lmpifort_nvidia -ldsmml -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lfftw3_mpi -lmpi_nvidia -lfftw3_threads -lfftw3 -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lrt -lgcc_s -lm -lstdc++ -lquadmath
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/core-240226/upgrade-nvhpc850-pcpu-20250226
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
