OPTS := $(OPTS) 
FOPTS =  $(OPTS) -c -fdefault-real-8 -Wall -cpp -DPETSC_VERSION=37 #-DRESTART_FACTOR
CCOPTS  = -c -O  -DPETSC_VERSION=37

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -O0
endif

MPI_DIR = /usr/local/openmpi/latest
GCC_DIR = /usr/lib/gcc/x86_64-linux-gnu/4.4.5

CC = $(MPI_DIR)/bin/mpicc
CPP = $(MPI_DIR)/bin/mpicxx
F90 = $(MPI_DIR)/bin/mpif90
F77 = $(MPI_DIR)/bin/mpif90
LOADER = $(MPI_DIR)/bin/mpif90 -gcc

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

MPI_DIR = /usr/local/openmpi/latest
HDF5_DIR =  $(PETSC_DIR)/$(PETSC_ARCH)
GSL_DIR = /usr/lib64
ZLIB_DIR = /usr/lib64

TRILINOS_DIR = /fasttmp/seol/openmpi-gcc4.4.5-install
STDCPP_DIR = /usr/lib/gcc/x86_64-linux-gnu/4.4.5

PETSC_DIR = /lore/seol/petsc-3.7.6
ifeq ($(COM), 1)
  PETSC_ARCH =cplx-openmpi
else
  PETSC_ARCH =real-openmpi
endif

INCLUDE = -I$(MPI_DIR)/include \
	  -I$(HDF5_DIR)/include \
	  -I$(PETSC_DIR)/include \
          -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
          -I$(SCOREC_DIR)/include

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(MPI_DIR)/lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lpetsc -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lsuperlu -lfftw3_mpi -lfftw3 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lX11 -lgfortran -lmpi_f90 -lmpi_f77 -lmpi_cxx -lstdc++ -L$(MPI_DIR)lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lm -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl

BLASLAPACK_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
	-lflapack -lfblas  \
        -Wl,--end-group

ifeq ($(COM), 1)
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
  endif
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif 
endif

PUMI_DIR=/lore/seol/openmpi-petsc3.7.6-install/debug
SCOREC_BASE_DIR=/lore/seol/openmpi-petsc3.7.6-install
ZOLTAN_DIR = /lore/seol/openmpi-petsc3.7.6-install
ifdef SCORECVER
  SCOREC_DIR=/lore/seol/openmpi-petsc3.7.6-install/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

SCOREC_LIBS= -L$(SCOREC_DIR)/lib -l$(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifeq ($(TRILINOS), 1)
TRILINOS_DIR = /fasttmp/seol/openmpi-gcc4.4.5-install
TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_DIR)/lib -L$(TRILINOS_DIR)/lib \
                -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
endif

LIBS = 	\
	$(SCOREC_LIBS) \
        $(TRILINOS_LIBS) \
        -L$(STDCPP_DIR) -lstdc++ \
        -L$(ZOLTAN_DIR)/lib -lzoltan \
        $(PETSC_WITH_EXTERNAL_LIB) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 \
	-L$(ZLIB_DIR) -lz \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas \
	-L/usr/lib64 -lX11

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
