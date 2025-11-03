FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=319 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=319
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -w -fallow-argument-mismatch
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g #not for gcc : noarg_temp_created 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
CC = cc
CPP = CC
F90 = ftn
F77 = ftn
LOADER =  ftn
endif

FOPTS := $(FOPTS)
F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606
ifeq ($(COM), 1)
  PETSC_ARCH=perlmuttercpu-gnu-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttercpu-gnu-cplx/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttercpu-gnu-cplx/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.30/ofi/gnu/12.3/lib -L/opt/cray/pe/mpich/8.1.30/ofi/gnu/12.3/lib -Wl,-rpath,/opt/cray/pe/libsci/24.07.0/GNU/12.3/x86_64/lib -L/opt/cray/pe/libsci/24.07.0/GNU/12.3/x86_64/lib -Wl,-rpath,/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -L/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -Wl,-rpath,/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/gnu/12.3/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/gnu/12.3/lib -Wl,-rpath,/opt/cray/pe/hdf5-parallel/1.12.2.9/gnu/12.3/lib -L/opt/cray/pe/hdf5-parallel/1.12.2.9/gnu/12.3/lib -Wl,-rpath,/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -L/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -Wl,-rpath,/opt/cray/pe/dsmml/0.3.0/dsmml/lib -L/opt/cray/pe/dsmml/0.3.0/dsmml/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/13 -L/usr/lib64/gcc/x86_64-suse-linux/13 -Wl,-rpath,/usr/x86_64-suse-linux/lib -L/usr/x86_64-suse-linux/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -ldarshan -llustreapi -lz -ldl -lxpmem -lgfortran -lm -lhdf5_hl_parallel -lhdf5_parallel -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lfftw3_mpi -lfftw3_threads -lfftw3 -lnetcdf -lnetcdff -lmpifort_gnu_123 -lsci_gnu_mpi -lmpi_gnu_123 -lsci_gnu -ldsmml -lgfortran -lquadmath -lpthread -lm -lgcc_s -lstdc++ -lquadmath
else
  PETSC_ARCH=perlmuttercpu-gnu
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttercpu-gnu/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606/perlmuttercpu-gnu/lib -Wl,-rpath,/opt/cray/pe/mpich/8.1.30/ofi/gnu/12.3/lib -L/opt/cray/pe/mpich/8.1.30/ofi/gnu/12.3/lib -Wl,-rpath,/opt/cray/pe/libsci/24.07.0/GNU/12.3/x86_64/lib -L/opt/cray/pe/libsci/24.07.0/GNU/12.3/x86_64/lib -Wl,-rpath,/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -L/opt/cray/pe/fftw/3.3.10.6/x86_milan/lib -Wl,-rpath,/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/gnu/12.3/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/gnu/12.3/lib -Wl,-rpath,/opt/cray/pe/hdf5-parallel/1.12.2.9/gnu/12.3/lib -L/opt/cray/pe/hdf5-parallel/1.12.2.9/gnu/12.3/lib -Wl,-rpath,/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -L/global/common/software/nersc9/darshan/3.4.6-gcc-13.2.1/lib -Wl,-rpath,/opt/cray/pe/dsmml/0.3.0/dsmml/lib -L/opt/cray/pe/dsmml/0.3.0/dsmml/lib -Wl,-rpath,/usr/lib64/gcc/x86_64-suse-linux/13 -L/usr/lib64/gcc/x86_64-suse-linux/13 -Wl,-rpath,/usr/x86_64-suse-linux/lib -L/usr/x86_64-suse-linux/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lkokkoskernels -lkokkoscontainers -lkokkoscore -lkokkossimd -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -ldarshan -llustreapi -lz -ldl -lxpmem -lgfortran -lm -lhdf5_hl_parallel -lhdf5_parallel -lhdf5hl_fortran_parallel -lhdf5_fortran_parallel -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lfftw3_mpi -lfftw3_threads -lfftw3 -lnetcdf -lnetcdff -lmpifort_gnu_123 -lsci_gnu_mpi -lmpi_gnu_123 -lsci_gnu -ldsmml -lgfortran -lquadmath -lpthread -lm -lgcc_s -lstdc++ -lquadmath
endif

SCOREC_BASE_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/core-240527/upgrade-gnu850-pcpu
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
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW_ROOT)/lib -lfftw3 -lfftw3_mpi

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include \
	-I$(FFTW_ROOT)/include

ifeq ($(ST), 1)
  LIBS += -Wl,--start-group -L$(HDF5_DIR)/lib -L$(NETCDF_DIR)/lib \
	  -Wl,-rpath,$(HDF5_ROOT)/lib -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran \
	  -Wl,-rpath,$(NETCDF_DIR)/lib -lnetcdf -lnetcdff \
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
