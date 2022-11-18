FOPTS = -c -fallow-argument-mismatch -DMPICH_IGNORE_CXX_SEEK -target-accel=nvidia80 -Ofast -fPIC -cpp -DPETSC_VERSION=316 -DUSEBLAS $(OPTS)
CCOPTS  = -c -DMPICH_IGNORE_CXX_SEEK -target-accel=nvidia80 -O -DPETSC_VERSION=316 -DDEBUG
#FOPTS = -c -r8 -i4 -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
#CCOPTS  = -c -O -DPETSC_VERSION=313 -DDEBUG

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

MPIVER=cuda11.7-mpich8.1.17
PETSCVER=petsc3.16.3
PETSC_VER=petsc-3.16.3

PETSC_DIR=/global/cfs/cdirs/mp288/scorec-pmt/petsc/$(PETSC_VER)
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-$(MPIVER)
else
  PETSC_ARCH=real-$(MPIVER)
endif

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/scorec-pmt/petsc/petsc-3.16.3/real-cuda11.7-mpich8.1.17/lib -L/global/cfs/cdirs/mp288/scorec-pmt/petsc/petsc-3.16.3/real-cuda11.7-mpich8.1.17/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.17/ofi/gnu/9.1/lib -L/opt/cray/pe/mpich/8.1.17/ofi/gnu/9.1/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.17/gtl/lib -L/opt/cray/pe/mpich/8.1.17/gtl/lib -Wl,-rpath,/opt/cray/pe/libsci/21.08.1.2/GNU/9.1/x86_64/lib -L/opt/cray/pe/libsci/21.08.1.2/GNU/9.1/x86_64/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/stubs -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/stubs -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/nvvm/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/nvvm/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/CUPTI/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/CUPTI/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/Debugger/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/Debugger/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/lib64 -Wl,-rpath,/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -Wl,-rpath,/opt/cray/xpmem/2.4.4-2.3_13.8__gff0e1d9.shasta/lib64 -L/opt/cray/xpmem/2.4.4-2.3_13.8__gff0e1d9.shasta/lib64 -Wl,-rpath,/opt/cray/pe/gcc/11.2.0/snos/lib/gcc/x86_64-suse-linux/11.2.0 -L/opt/cray/pe/gcc/11.2.0/snos/lib/gcc/x86_64-suse-linux/11.2.0 -Wl,-rpath,/opt/cray/pe/gcc/11.2.0/snos/lib64 -L/opt/cray/pe/gcc/11.2.0/snos/lib64 -Wl,-rpath,/opt/cray/pe/gcc/11.2.0/snos/lib -L/opt/cray/pe/gcc/11.2.0/snos/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lgsl -lgslcblas -lcuda -lcudart -lcufft -lcublas -lcusparse -lcusolver -lcurand -lstdc++ -ldl -lcuda -lmpi_gtl_cuda -lxpmem -lgfortran -lm -lcupti -lcudart -lsci_gnu_82_mpi -lsci_gnu_82 -lmpifort_gnu_91 -lmpi_gnu_91 -ldsmml -lgfortran -lquadmath -lpthread -lm -lgcc_s -lquadmath -lstdc++ -ldl

ZOLTAN_DIR=/global/cfs/cdirs/mp288/scorec-pmt/$(MPIVER)/$(PETSCVER)
SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/scorec-pmt/$(MPIVER)/$(PETSCVER)
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

ZOLTAN_LIB=-lzoltan

LIBS =  \
        $(SCOREC_LIB) \
        -L$(ZOLTAN_DIR)/lib $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB)

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	 -I$(SCOREC_DIR)/include \
	 -I/opt/cray/pe/hdf5-parallel/1.12.1.5/gnu/9.1/include

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
