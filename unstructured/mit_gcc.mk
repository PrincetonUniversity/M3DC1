FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=313 -DDEBUG

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -w -fallow-argument-mismatch -O2 
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

PETSC_DIR=/orcd/nese/psfc/001/jinchen/petsc/petsc20230612
ifeq ($(COM), 1)
  PETSC_ARCH=mit-gcc-openmpi-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi-cplx/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi-cplx/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
else
  ifeq ($(ST), 1)
  PETSC_ARCH=nvidia-hpc-sdk-6-st
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220915/nvidia-hpc-sdk-6-st/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220915/nvidia-hpc-sdk-6-st/lib -L/opt/cray/pe/mpich/8.1.17/ofi/nvidia/20.7/lib -L/opt/cray/pe/mpich/8.1.17/gtl/lib -L/opt/cray/pe/libsci/21.08.1.2/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/stubs -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/nvvm/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/CUPTI/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/Debugger/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/lib64 -L/opt/cray/xpmem/2.4.4-2.3_13.8__gff0e1d9.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/compilers/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lstdc++ -ldl -lcupti -lcudart -lcuda -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -lmpi_gtl_cuda -ldsmml -lxpmem -lacchost -laccdevaux -laccdevice -lcudadevice -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lquadmath -lstdc++ -ldl
  else
  PETSC_ARCH=mit-gcc-openmpi
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
  endif
endif

SCOREC_BASE_DIR=/net/eofe-netapp-nfs/home/jinch731/core/mit-gcc-openmpi
#SCOREC_BASE_DIR=/orcd/nese/psfc/001/software/scorec/gcc12.2.0-openmpi4.1.4/petsc3.19.2
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

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \

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
