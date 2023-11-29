FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=313 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=313
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

PETSC_DIR=/orcd/nese/psfc/001/jinchen/petsc/petsc20230612
ifeq ($(COM), 1)
  PETSC_ARCH=mit-intel-impi-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi-cplx/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi-cplx/lib -Wl,-rpath,/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -L/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lfftw3_mpi -lfftw3 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -lgfortran -lstdc++ -lquadmath
else
  ifeq ($(ST), 1)
  PETSC_ARCH=mit-intel-impi-st
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi/lib -Wl,-rpath,/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -L/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lfftw3_mpi -lfftw3 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -lgfortran -lstdc++ -lquadmath
  else
  PETSC_ARCH=mit-intel-impi
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-intel-impi/lib -Wl,-rpath,/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -L/orcd/nese/psfc/001/software/spack/2023-07-01-physics-rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/intel-oneapi-mkl-2023.1.0-seow5nciajy2dbojmltvs7abacdaudve/mkl/2023.1.0/lib/intel64 -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lfftw3_mpi -lfftw3 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -lgfortran -lstdc++ -lquadmath
  endif
endif

SCOREC_BASE_DIR=$(PETSC_DIR)/$(PETSC_ARCH)
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
  LIBS += -Wl,--start-group -L$(NETCDF_FDIR)/lib -Wl,-rpath,$(NETCDF_FDIR)/lib -lnetcdff -lz -Wl,--end-group -Wl,--start-group -L$(NETCDF_CDIR)/lib -Wl,-rpath,$(NETCDF_CDIR)/lib -lnetcdf -Wl,--end-group \
          -Wl,--start-group -L$(NETCDF_FDIR)/lib -Wl,-rpath,$(NETCDF_FDIR)/lib -lnetcdff -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -Wl,--end-group
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
