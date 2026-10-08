FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp $(OPTS) -DUSEBLAS -DPETSC_VERSION=313
CCOPTS  = -c -DPETSC_VERSION=313
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -w -fallow-argument-mismatch
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g 
  CCOPTS := $(CCOPTS) -g 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

CC = mpicc
CPP = mpicxx
F90 = mpif90
F77 = mpif90
LOADER = mpif90
LDOPTS := $(LDOPTS)
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
PETSC_DIR=/projects/M3DC1/PETSC/petsc.20220201
ifeq ($(COM), 1)
  $(error COM=1 (complex build) is not ready yet on stellarai)
else
  PETSC_ARCH=real-stellarai-gcc14-cuda100-gO2
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/projects/M3DC1/PETSC/petsc.20220201/real-stellarai-gcc14-cuda100-gO2/lib -L/projects/M3DC1/PETSC/petsc.20220201/real-stellarai-gcc14-cuda100-gO2/lib -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -Wl,-rpath,/usr/local/cuda-13.2/lib64 -L/usr/local/cuda-13.2/lib64 -L/usr/local/cuda-13.2/lib64/stubs -Wl,-rpath,/usr/local/openmpi/cuda-13.2/5.0.10/gcc/lib64 -L/usr/local/openmpi/cuda-13.2/5.0.10/gcc/lib64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/14 -L/usr/lib/gcc/x86_64-redhat-linux/14 -Wl,-rpath,/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -L/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -Wl,-rpath,/usr/lib/gcc -L/usr/lib/gcc -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -lpetsc -lfftw3_mpi -lfftw3 -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lkokkoskernels -lsuperlu -lflame -lblis -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lz -lcudart -lnvtx3interop -lcufft -lcublas -lcusparse -lcusolver -lcurand -lcuda -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath
endif

SCOREC_BASE_DIR=/projects/M3DC1/scorec/stellarai-amd-gcc14-cuda100
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
SIMMETRIX_VER=2024.0-231117dev
#MESHGEN_DIR=/projects/M3DC1/scorec/stellar/$(MPIVER)/$(SIMMETRIX_VER)/bin

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

# cuDSS block-Jacobi solver (matrix_solve::solve_cudss, runtime option
# -cudsssolve <matrix_id>).  Requires m3dc1_scorec built with ENABLE_CUDSS=ON
# and libm3dc1_cudss.a from m3d-c1_preconditioner_standalone (configured with
# stellarai.sh there).  Vendored cuDSS ships only libcudss.so.0, so link it by
# full path.  No NCCL module on stellar-ai: export NCCL_DIR by hand.
# cudart/cusparse/cublas are already linked via PETSc.
#ifeq ($(CUDSS), 1)
  CUDSS_PKG_DIR = /scratch/gpfs/CSD/jinchen/petsc
  CUDSS_DIR ?= $(CUDSS_PKG_DIR)/deps/cudss
  CUDSS_LIB = -L$(SCOREC_DIR)/lib -lm3dc1_cudss \
              -Wl,-rpath,$(CUDSS_DIR)/lib $(CUDSS_DIR)/lib/libcudss.so.0 \
              -Wl,-rpath,$(CUDSS_PKG_DIR)/lib -L$(CUDSS_PKG_DIR)/lib -lnccl
#else
#  CUDSS_LIB =
#endif

HDF5DIR=$(PETSC_DIR)/$(PETSC_ARCH)
LIBS = 	$(SCOREC_LIBS) \
	$(CUDSS_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(HDF5DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran \

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_DIR)/include \
	-I$(HDF5DIR)/include \

ifeq ($(ST), 1)
  NETCDFDIR=$(PETSC_DIR)/$(PETSC_ARCH)
  LIBS += -L$(NETCDFDIR)/lib -Wl,-Bstatic -lnetcdf -Wl,-Bdynamic -lzip
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
