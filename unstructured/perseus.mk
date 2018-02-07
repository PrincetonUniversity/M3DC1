FOPTS = -c -r8 -implicitnone -fpp -warn all -DNEXTPetscDEV -DPETSC_VERSION=38 $(OPTS)
CCOPTS  = -c -O -DPETSC_VERSION=38 -DNEXTPetscDEV

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -vec-report0 # -fast
#  FOPTS  := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
endif

ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=../select.tau
  CC     = tau_cc.sh $(TAU_OPTIONS)
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CC = mpicc
  CPP = mpicxx
  F90 = mpif90
  F77 = mpif90
  LOADER = mpif90 -cxxlib
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define petsc and other libs
PETSC_DIR = /home/jinchen/LIB/petsc-3.8.2
ifeq ($(COM), 1)
  PETSC_ARCH = complex-intel17.0-openmpi1.10.2
else
  PETSC_ARCH = real-intel17.0-openmpi1.10.2
endif

LIB64_DIR=/usr/local/intel/lib64
MPI_DIR=/usr/local/openmpi/1.10.2/intel170/x86_64/lib64
MKL_DIR=/opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/lib/intel64
COMPILER_DIR=/opt/intel/compilers_and_libraries_2017.5.239/linux/compiler/lib/intel64
COMPILER_LIN_DIR=/opt/intel/compilers_and_libraries_2017.5.239/linux/compiler/lib/intel64_lin
ifeq ($(COM), 1)
  HYPRE_LIB=
else
  HYPRE_LIB=-lHYPRE
endif

PETSC_WITH_EXTERNAL_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(LIB64_DIR) -L$(LIB64_DIR)/openmpi -L$(MPI_DIR) -L$(COMPILER_DIR) -L$(MKL_DIR) -L$(COMPILER_LIN_DIR) -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -Wl,-rpath,$(LIB64_DIR) -Wl,-rpath,$(LIB64_DIR)/openmpi -Wl,-rpath,$(MPI_DIR) -lpetsc -lsuperlu -lsuperlu_dist $(HYPRE_LIB) -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lfftw3_mpi -lfftw3 -lflapack -lfblas -lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lX11 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lifport -lifcoremt_pic -lmpi_cxx -lintlc -lrt -lm -lpthread -lz -ldl -lstdc++ -lmpi -limf -lsvml -lirng -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lpthread -lirc_s 

SCOREC_DIR = /home/jinchen/LIB/scorec/intel17.0-openmpi1.10.2/Jan2018
ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif
SCOREC_LIBS=-Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
            $(M3DC1_SCOREC_LIB) -lpumi -lcrv -lph -lsam -lspr -lma \
            -lapf_zoltan -lparma -lmds -lapf -llion -lmth -lgmi -lpcu \
            -Wl,--end-group

ZOLTAN_DIR=/home/jinchen/LIB/scorec/intel17.0-openmpi1.10.2
ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan

INCLUDE = -I$(MPIHOME)/include  -I$(SCOREC_DIR)/include \
          -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

LIBS = $(SCOREC_LIBS) \
       $(ZOLTAN_LIB) \
       $(PETSC_WITH_EXTERNAL_LIB) \
       -L/usr/local/gsl/2.4/x86_64/lib64  -lssl -lgsl -lgslcblas

LDOPTS := $(LDOPTS) -qopenmp -lmpi_cxx
FOPTS  := $(FOPTS)  -qopenmp
CCOPTS := $(CCOPTS) -qopenmp

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
