  CPP = mpicxx
  CC = mpicc
  F90 = mpifort
  F77 = mpif90
  LOADER = mpif90

OPTS := $(OPTS) -DPETSC_VERSION=313 -DUSEBLAS

MPIVER=nvhpc20.7-cuda11.0-openmpi4.0.4
PETSCVER=petsc3.13.5
PETSC_VER=petsc-3.13.5

PETSC_DIR=/projects/M3DC1/PETSC/$(PETSC_VER)
ifeq ($(COM), 1)
   PETSC_ARCH=cplx-$(MPIVER)
   M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
   PETSC_ARCH=real-$(MPIVER)
   M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_BASE_DIR=/projects/M3DC1/scorec/$(MPIVER)/$(PETSCVER)
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(SCOREC_DIR)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif
		
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/projects/M3DC1/PETSC/petsc-3.13.5/real-nvhpc20.7-cuda-11.0-openmpi4.0.4/lib -L/projects/M3DC1/PETSC/petsc-3.13.5/real-nvhpc20.7-cuda-11.0-openmpi4.0.4/lib -L/usr/local/nvhpc/lib64 -L/usr/local/nvhpc/lib64/openmpi -L/usr/local/openmpi/cuda-11.0/4.0.4/nvhpc207/ppc64le/lib64 -L/opt/nvidia/hpc_sdk/Linux_ppc64le/20.7/compilers/lib -L/usr/lib/gcc/ppc64le-redhat-linux/8 -Wl,-rpath,/usr/local/nvhpc/lib64 -Wl,-rpath,/usr/local/nvhpc/lib64/openmpi -Wl,-rpath,/usr/local/openmpi/cuda-11.0/4.0.4/nvhpc207/ppc64le/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_ppc64le/20.7/compilers/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lparmetis -lmetis -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lnvf -lnvomp -latomic -lnvhpcatm -lpthread -lnvcpumath -lnvc -lrt -lm -lgcc_s -lstdc++ -ldl

#only define them if adios-1.3 is used; otherwise use hopper default
INCLUDE := $(INCLUDE) -I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include -I/home/jinchen/LIB/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(HDF5_ROOT)/include
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
#CUDA_LIB:=-L/usr/local/cuda-11.0/lib64 -lcudart -lcusparse -lcusolver -lstdc++ -L/usr/local/cuda-11.0/lib64 -lcublas
LIBS := $(LIBS) \
	$(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-lgsl -lgslcblas

#        $(CUDA_LIB) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-lgsl -lgslcblas 

FOPTS = -c -r8 -Mpreprocess $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
# FIXME 
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -fast#-static -qopt-report
  FOPTS  := $(FOPTS)  -fast#-qopt-report
  CCOPTS := $(CCOPTS) -fast #-qopt-report
else
  FOPTS := $(FOPTS) -Mbounds -Minfo=all -Mchkfpstk -Mchkstk -Mdalign -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
  CCOPTS := $(CCOPTS) -g 
  LDOPTS := $(LDOPTS) -g
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp 
  FOPTS  := $(FOPTS)  -mp 
  CCOPTS := $(CCOPTS) -mp 
endif

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@

