OPTS := $(OPTS) -DUSEBLAS
FOPTS = -c -fdefault-real-8 -Wall -cpp -DPETSC_VERSION=37 -DUSESCOREC $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=37

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g noarg_temp_created 
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

PETSC_DIR = /lore/seol/petsc-3.7.6
ifeq ($(COM), 1)
  PETSC_ARCH =complex-openmpi1.6.5
else
  PETSC_ARCH =real-openmpi
endif


SCOREC_DIR = /lore/seol/openmpi-petsc3.7.6-install
ZOLTAN_DIR = /lore/seol/openmpi-petsc3.7.6-install

HDF5_DIR =  $(PETSC_DIR)/$(PETSC_ARCH)
GSL_DIR = /usr/lib64
ZLIB_DIR = /usr/lib64
FFTW_DIR = /usr/lib64

TRILINOS_DIR = /fasttmp/seol/openmpi-gcc4.4.5-install

INCLUDE = -I$(MPI_DIR)/include \
	  -I$(HDF5_DIR)/include \
	  -I$(PETSC_DIR)/include \
          -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
          -I$(SCOREC_DIR)/include

ifeq ($(COM), 1)
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib  -lpetsc -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist_3.3 -lsuperlu_4.3 -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L$(MPI_DIR)/lib -L/usr/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -ldl -lmpi -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl
else
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(MPI_DIR)/lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lpetsc -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lsuperlu -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lX11 -lgfortran -lmpi_f90 -lmpi_f77 -lmpi_cxx -lstdc++ -L$(MPI_DIR)lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lm -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl
endif

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

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

ifeq ($(TRILINOS), 1)
TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_DIR)/lib -L$(TRILINOS_DIR)/lib \
                -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib  -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L$(MPI_DIR)/lib -L/usr/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_DIR)/lib -L$(GCC_DIR) -ldl -lmpi -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl
endif

LIBS = 	\
	$(SCOREC_LIBS) \
        $(TRILINOS_LIBS) \
        -L$(STDCPP_DIR) -lstdc++ \
        -L$(ZOLTAN_DIR)/lib -lzoltan \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 \
        -L$(FFTW_DIR)/lib -lfftw3 -lfftw3_threads \
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
