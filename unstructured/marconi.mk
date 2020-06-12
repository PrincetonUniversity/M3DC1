ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpiicpc
  CC = mpiicc
  F90 = mpiifort
  F77 = mpiifort
  LOADER = mpiifort
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif
 
OPTS := $(OPTS) -DPETSC_VERSION=39 -DUSEBLAS

PUMI_DIR=/marconi/home/userexternal/jchen000/project/PUMI/core/build
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
SCOREC_UTIL_DIR=$(PUMI_LIB)/bin

ZOLTAN_DIR=/marconi_work/FUA34_WPJET1/jchen000/ZOLTAB/Zoltan_v3.8/build
ZOLTAN_LIB=-L$(ZOLTAN_DIR)/lib -lzoltan

ifdef SCORECVER
  SCOREC_DIR=$(PUMI_LIB)/$(SCORECVER)
else
  SCOREC_DIR=/marconi/home/userexternal/jchen000/project/SRC/M3DC1/m3dc1_scorec/build-real
endif

ifeq ($(COM), 1)
  SCOREC_DIR=/marconi/home/userexternal/jchen000/project/SRC/M3DC1/m3dc1_scorec/build-cplx
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  SCOREC_DIR=/marconi/home/userexternal/jchen000/project/SRC/M3DC1/m3dc1_scorec/build-real
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

PETSC_DIR=/marconi/home/userexternal/jchen000/project/PETSC/master
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-intel-intelmpi
else
  PETSC_ARCH=real-intel-intelmpi
endif

MKL_LIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
#PETSC_WITH_EXTERNAL_LIB = -Wl,--start-group -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist  -lparmetis -lmetis -lfblas -lflapack -lscalapack -Wl,--end-group -lstdc++
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -L/cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.4.274/intel64/lib/release_mt -L/cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.4.274/intel64/lib -L/marconi/prod/opt/compilers/intel/pe-xe-2018/binary/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64_lin -L/cineca/prod/opt/compilers/gnu/7.4.0/none/lib64/ -Wl,-rpath,/cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.4.274/intel64/lib/release_mt -Wl,-rpath,/cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.4.274/intel64/lib -Wl,-rpath,/opt/intel/mpi-rt/2017.0.0/intel64/lib/release_mt -Wl,-rpath,/opt/intel/mpi-rt/2017.0.0/intel64/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lflapack  -lparmetis -lmetis -ldl -lstdc++ -lmpifort -lmpi -lmpigi -lrt -lpthread -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lquadmath -ldl -lstdc++

#GSL_DIR=/marconi/home/userexternal/jchen000/project/GSL/build

INCLUDE := $(INCLUDE) -I$(PUMI_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_HOME)/include \
	   -I$(HDF5_HOME)/include \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
        $(SCOREC_LIB) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	-L$(FFTW_HOME)/lib -lfftw3_mpi -lfftw3 \
        -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 \
	-L$(GSL_HOME)/lib -lgsl \
	$(MKL_LIB)
#	-L/usr/lib64/ -ldl -lrt -lpthread -lm -lc
#        $(HYBRID_LIBS) \

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -assume no2underscores -ftz -fPIC -O
  FOPTS  := $(FOPTS)  
  CCOPTS := $(CCOPTS) -std=c++14
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -fopenmp 
  FOPTS  := $(FOPTS)  -fopenmp 
  CCOPTS := $(CCOPTS) -fopenmp 
endif

F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
