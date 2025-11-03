FOPTS = -c -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS)
CCOPTS  = -c -O -DPETSC_VERSION=313 -DDEBUG
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -w -O2 #-qopt-report=5 -qopt-report-phase=vec,loop
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g 
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
CPP = mpicxx
F90 = mpif90
F77 = mpif77
LOADER =  mpif90
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/orcd/pool/008/jinch731/petsc20250515
ifeq ($(COM), 1)
  PETSC_ARCH=mit-gcc-ompi-cplx
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/orcd/pool/008/jinch731/petsc20250515/mit-gcc-ompi-cplx/lib -L/orcd/pool/008/jinch731/petsc20250515/mit-gcc-ompi-cplx/lib -Wl,-rpath,/orcd/software/core/001/pkg/openmpi/5.0.8/lib -L/orcd/software/core/001/pkg/openmpi/5.0.8/lib -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib64 -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib64 -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib -lpetsc -lfftw3_mpi -lfftw3 -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lflapack -lfblas -lzoltan -lnetcdf -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl_f90cstub -lhdf5_f90cstub -lhdf5_hl -lhdf5 -lparmetis -lmetis -lgsl -lgslcblas -lz -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
else
  PETSC_ARCH=mit-gcc-ompi
  PETSC_WITH_EXTERNAL_LIB = -Wl,-rpath,/orcd/pool/008/jinch731/petsc20250515/mit-gcc-ompi/lib -L/orcd/pool/008/jinch731/petsc20250515/mit-gcc-ompi/lib -Wl,-rpath,/orcd/software/core/001/pkg/openmpi/5.0.8/lib -L/orcd/software/core/001/pkg/openmpi/5.0.8/lib -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib64 -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib64 -Wl,-rpath,/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib -L/orcd/software/core/001/spack/pkg/gcc/12.2.0/yt6vabm/lib -lpetsc -lfftw3_mpi -lfftw3 -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lsuperlu -lflapack -lfblas -lzoltan -lnetcdf -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl_f90cstub -lhdf5_f90cstub -lhdf5_hl -lhdf5 -lparmetis -lmetis -lgsl -lgslcblas -lz -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
endif

SCOREC_BASE_DIR=/orcd/pool/008/jinch731/pumi/gcc-ompi
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
M3DC1_SCOREC_DIR=$(SCOREC_BASE_DIR)

ifdef SCORECVER
  SCOREC_DIR=$(M3DC1_SCOREC_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(M3DC1_SCOREC_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

LIBS = 	\
	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \

ifeq ($(ST), 1)
  NETCDF_CDIR=$(PETSC_DIR)/$(PETSC_ARCH)
  NETCDF_FDIR=$(PETSC_DIR)/$(PETSC_ARCH)
  LIBS += -Wl,--start-group -L$(NETCDF_CDIR)/lib -Wl,-rpath,$(NETCDF_CDIR)/lib -lnetcdff -Wl,--end-group \
          -Wl,--start-group -L$(NETCDF_FDIR)/lib -Wl,-rpath,$(NETCDF_FDIR)/lib -lnetcdf -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl_f90cstub -lhdf5_f90cstub -lhdf5_hl -lhdf5 -lz -Wl,--end-group
  INCLUDE += -I$(NETCDF_CDIR)/include \
             -I$(NETCDF_FDIR)/include
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
