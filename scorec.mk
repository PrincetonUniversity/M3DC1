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

MPIHOME = /usr/local/openmpi/latest
ZOLTAN_HOME = /usr/local/zoltan/latest/openmpi1.6.5_gcc4.4.5
PETSC_DIR = /lore/fzhang/petsc-real/petsc-3.5.1
PETSC_ARCH = test
HDF5_HOME = /lore/fzhang/petsc-real/petsc-3.5.1/test
GSL_HOME = /usr/lib64
ZLIB_HOME = /usr/lib64
FFTW_HOME = /usr/lib64
SCOREC_HOME = /users/seol/develop
TRILINOS_HOME = /usr/local/trilinos/latest
STDCPP_HOME = /users/granzb/lib64

INCLUDE = -I$(MPIHOME)/include \
	  -I$(PETSC_DIR)/include \
          -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	  -I$(HDF5_HOME)/include \
          -I$(SCOREC_HOME)/include

ifeq ($(COM), 1)
PETSC_ARCH = portalr6-intel-openmpi-1.8.4-complex
PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
	-ldmumps -lmumps_common -lpord -lcmumps -lsmumps -lzmumps \
	-lparmetis -lmetis \
	-lscalapack \
	-lsuperlu_dist_3.3 -lsuperlu_4.3 \
	-Wl,--end-group
else
PETSC_ARCH = portalr6-intel-openmpi-1.8.4
PETSC_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
        -lpetsc \
	-ldmumps -lmumps_common -lpord -lcmumps -lsmumps -lzmumps \
	-lpord \
	-lparmetis -lmetis \
	-lscalapack \
	-lsuperlu_dist_3.3 -lsuperlu_4.3 \
	-Wl,--end-group
endif

BLASLAPACK_LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,--start-group \
	-lfblas -lflapack \
        -Wl,--end-group

ifeq ($(COM), 1)
  SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_HOME)/lib -L$(SCOREC_HOME)/lib \
               -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lm3dc1_scorec_complex \
               -Wl,--end-group
else
  SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_HOME)/lib -L$(SCOREC_HOME)/lib \
               -lapf -lgmi -lma -lparma -lph -lmds -lpcu -lspr -lapf_zoltan -lm3dc1_scorec_trilinos \
               -Wl,--end-group
endif

TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_HOME)/lib -L$(TRILINOS_HOME)/lib \
                -lamesos2 -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group

LIBS = 	\
	$(SCOREC_LIBS) \
        $(TRILINOS_LIBS) \
        -L$(STDCPP_HOME) -lstdc++ \
        -L$(ZOLTAN_HOME)/lib -lzoltan \
        $(PETSC_LIBS) \
        $(BLASLAPACK_LIBS) \
	-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 \
	-Wl,-rpath -Wl,$(HDF5_HOME)/lib \
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
