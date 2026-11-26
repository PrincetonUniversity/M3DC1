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
  PETSC_ARCH=cplx-stellarai-gcc14-gO2
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex 
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/projects/M3DC1/PETSC/petsc.20220201/cplx-stellarai-gcc14-gO2/lib -L/projects/M3DC1/PETSC/petsc.20220201/cplx-stellarai-gcc14-gO2/lib -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -Wl,-rpath,/usr/local/openmpi/5.0.10/gcc/lib64 -L/usr/local/openmpi/5.0.10/gcc/lib64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/14 -L/usr/lib/gcc/x86_64-redhat-linux/14 -Wl,-rpath,/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -L/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -Wl,-rpath,/usr/lib/gcc -L/usr/lib/gcc -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -lpetsc -lfftw3_mpi -lfftw3 -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lflame -lblis -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lz -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath

else
  PETSC_ARCH=real-stellarai-gcc14-gO2
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/projects/M3DC1/PETSC/petsc.20220201/real-stellarai-gcc14-gO2/lib -L/projects/M3DC1/PETSC/petsc.20220201/real-stellarai-gcc14-gO2/lib -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib -Wl,-rpath,/usr/local/openmpi/5.0.10/gcc/lib64 -L/usr/local/openmpi/5.0.10/gcc/lib64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/14 -L/usr/lib/gcc/x86_64-redhat-linux/14 -Wl,-rpath,/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -L/usr/local/hdf5/gcc/openmpi-5.0.10/1.14.6/lib64 -Wl,-rpath,/usr/lib/gcc -L/usr/lib/gcc -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -L/opt/AMD/aocl/aocl-linux-gcc-5.3.0/gcc/ST/lib_LP64 -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_spectrum-x_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/nccl_rdma_sharp_plugin/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/sharp/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/hcoll/lib -Wl,-rpath,/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -L/opt/hpcx-v2.51-gcc-doca_ofed-redhat10-cuda13-x86_64/ucc/lib -lpetsc -lfftw3_mpi -lfftw3 -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lflame -lblis -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lz -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -lquadmath
endif

SCOREC_BASE_DIR=/projects/M3DC1/scorec/stellarai-amd-gcc14
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

LIBS = 	$(SCOREC_LIBS) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(HDF5DIR)/lib64 -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran \

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
