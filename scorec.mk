FOPTS = -c  -fdefault-real-8 -Wall -cpp -DPetscDEV -DPETSC_31 -DKSPITS $(OPTS) -DLATESTSCOREC
CCOPTS  = -c -O -DPetscDEV -DPETSC_31 -DPetscOLD #-DCJ_MATRIX_DUMP -DUSEHYBRID 

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g noarg_temp_created 
endif

  CC = /usr/local/openmpi/latest/bin/mpicc
  CPP = /usr/local/openmpi/latest/bin/mpicxx
  F90 = /usr/local/openmpi/latest/bin/mpif90
  F77 = /usr/local/openmpi/latest/bin/mpif90
  LOADER = /usr/local/openmpi/latest/bin/mpif90 -gcc
  FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)


ifeq ($(COM), 1)
  PETSC_DIR = /fasttmp/seol/petsc-3.5.4-complex
  else
  PETSC_DIR = /fasttmp/seol/petsc-3.5.4-real
endif
PETSC_ARCH = openmpi1.6.5

MPI_HOME = /usr/local/openmpi/latest
ZOLTAN_HOME = $(PETSC_DIR)/$(PETSC_ARCH)
HDF5_HOME =  $(PETSC_DIR)/$(PETSC_ARCH)
GSL_HOME = /usr/lib64
ZLIB_HOME = /usr/lib64
FFTW_HOME = /usr/lib64
SCOREC_HOME = /users/seol/public
TRILINOS_HOME = /usr/local/trilinos/latest
STDCPP_HOME = /users/granzb/lib64

INCLUDE = -I$(MPI_HOME)/include \
	  -I$(HDF5_HOME)/include \
	  -I$(PETSC_DIR)/include \
          -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
          -I$(SCOREC_HOME)/include

ifeq ($(COM), 1)
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib  -lpetsc -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist_3.3 -lsuperlu_4.3 -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lmpi_cxx -lstdc++ -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl
else
  PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib  -lpetsc -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist_3.3 -lsuperlu_4.3 -lHYPRE -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lmpi_cxx -lstdc++ -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -L$(MPI_HOME)/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lnuma -lrt -lnsl -lutil -lgcc_s -lpthread -ldl
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

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_HOME)/lib -L$(SCOREC_HOME)/lib \
               -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -l$(M3DC1_SCOREC_LIB) \
               -Wl,--end-group
ifeq ($(TRILINOS), 1)
TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_HOME)/lib -L$(TRILINOS_HOME)/lib \
                -lamesos2 -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
endif

LIBS = 	\
	$(SCOREC_LIBS) \
        $(TRILINOS_LIBS) \
        -L$(STDCPP_HOME) -lstdc++ \
        -L$(ZOLTAN_HOME)/lib -lzoltan \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
        -L$(FFTW_HOME)/lib -lfftw3 -lfftw3_threads \
	-L$(ZLIB_HOME) -lz \
	-L$(GSL_HOME)/lib -lgsl -lgslcblas \
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
