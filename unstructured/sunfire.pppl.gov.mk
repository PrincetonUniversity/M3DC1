FOPTS = $(OPTS) -DPETSC_VERSION=37 -c -r8 -implicitnone -fpp -warn all -DLATESTSCOREC -DUSEBLAS
# FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DLATESTSCOREC -DUSEPARTICLES
CCOPTS  = -c -DPETSC_VERSION=37

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
  CCOPTS := $(CCOPTS) -g -check-uninit -debug all
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
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
  LOADER = mpif90
  LDOPTS := $(LDOPTS) -cxxlib
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver
PETSC_VER=petsc-3.7.6
PETSCVER=petsc3.7.6

PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/$(PETSC_VER)
ifeq ($(COM), 1)
PETSC_ARCH=complex-openmpi-3.0.0
HYPRE_LIB=
else
PETSC_ARCH=real-intel2018-openmpi3.0.0-gcc6.1.0
HYPRE_LIB=-lHYPRE
endif

SCOREC_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/intel2018-openmpi3.0.0-gcc6.1.0/$(PETSCVER)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

SCOREC_UTIL_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/openmpi-3.0.0/bin

ifeq ($(TRILINOS), 1)
  TRILINOS_DIR=/usr/pppl/intel/2015-pkgs/openmpi-1.10.3-pkgs/trilinos-11.12.1
  ZOLTAN_LIB=-L$(TRILINOS_DIR)/lib -lzoltan
  TRILINOS_LIBS = -Wl,--start-group,-rpath,$(TRILINOS_DIR)/lib -L$(TRILINOS_DIR)/lib \
        -lstdc++  -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra -lepetraext \
        -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics -lteuchosremainder
  PETSC_LIBS=
else
  ZOLTAN_LIB=-L/p/tsc/m3dc1/lib/SCORECLib/rhel6/intel2018-openmpi3.0.0-gcc6.1.0/lib -lzoltan
  TRILINOS_LIBS=
  PETSC_LIBS =-L/p/tsc/m3dc1/lib/SCORECLib/rhel6/petsc-3.7.6/real-intel2018-openmpi3.0.0-gcc6.1.0/lib -Wl,-rpath,/p/tsc/m3dc1/lib/SCORECLib/rhel6/petsc-3.7.6/real-intel2018-openmpi3.0.0-gcc6.1.0/lib -L/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/pppl/gcc/6.1.0/lib/gcc/x86_64-pc-linux-gnu/6.1.0 -L/usr/pppl/gcc/6.1.0/lib64 -L/usr/pppl/gcc/6.1.0/lib -Wl,-rpath,/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -lpetsc -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lsuperlu -lHYPRE -lscalapack -lfftw3_mpi -lfftw3 -lflapack -lfblas -lhwloc -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lX11 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lifport -lifcoremt_pic -lintlc -ldl -L/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -lmpi -L/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/gcc/6.1.0/lib/gcc/x86_64-pc-linux-gnu/6.1.0 -L/usr/pppl/gcc/6.1.0/lib64 -L/usr/pppl/gcc/6.1.0/lib64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/pppl/gcc/6.1.0/lib -Wl,-rpath,/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -limf -lsvml -lirng -lm -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lpthread -lirc_s -L/usr/pppl/intel/2018-pkgs/openmpi-3.0.0/lib -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/gcc/6.1.0/lib/gcc/x86_64-pc-linux-gnu/6.1.0 -L/usr/pppl/gcc/6.1.0/lib64 -L/usr/pppl/gcc/6.1.0/lib64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.4 -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/usr/pppl/intel/2018.u1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64_lin/gcc4.4 -L/usr/pppl/gcc/6.1.0/lib -ldl 
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB=-lm3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB=-lm3dc1_scorec
  endif
endif

SCORECLIB= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
           $(PUMI_LIB) $(M3DC1_SCOREC_LIB) -Wl,--end-group

LIBS = 	\
	$(SCORECLIB) \
        $(TRILINOS_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_LIBS) \
	-L$(GSL_HOME)/lib -lgsl -lgslcblas \
	-lX11

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(GSL_HOME)/include

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
