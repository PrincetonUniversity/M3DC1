FOPTS = -c -cpp -DPETSC_VERSION=318 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=318 -DDEBUG

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

MPIVER=nvidia8.3.3-mpich8.1.22
PETSC_VER=petsc-3.18.4
PETSCVER=petsc3.18.4
PETSC_DIR=/global/cfs/cdirs/mp288/scorec-pmt/petsc/$(PETSC_VER)

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/scorec-pmt/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ZOLTAN_DIR=$(SCOREC_BASE_DIR)
PUMI_DIR=$(SCOREC_BASE_DIR)

PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  PETSC_ARCH=cplx-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  PETSC_ARCH=real-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/opt/cray/pe/mpich/8.1.22/ofi/nvidia/20.7/lib -L/opt/cray/pe/mpich/8.1.22/gtl/lib -L/opt/cray/pe/libsci/22.11.1.2/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/fftw/3.3.10.2/x86_milan/lib -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/stubs -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/nvvm/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/CUPTI/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/extras/Debugger/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/lib64 -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/opt/cray/xpmem/2.5.2-2.4_3.20__gd0f7936.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64 -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lz -lstdc++ -ldl -lcupti -lcudart -lcuda -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lfftw3_mpi -lfftw3_threads -lfftw3 -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -lmpi_gtl_cuda -ldsmml -lxpmem -lacchost -laccdevaux -laccdevice -lcudadevice -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lquadmath -lstdc++ -ldl

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan

GSL_DIR=$(GSL_ROOT)

INCLUDE := $(INCLUDE) -I$(GSL_DIR)/include \
       -I$(SCOREC_DIR)/include \
       -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
       -I$(PETSC_DIR)/include \
       -I$(GSL_DIR)/include

LIBS = 	\
	$(SCOREC_LIB) $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas


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
