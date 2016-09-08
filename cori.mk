ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = CC
  CC = cc
  F90 = ftn
  F77 = ftn
  LOADER = ftn
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

# define where you want to locate the mesh adapt libraries
#HYBRID_HOME =  /scratch2/scratchdirs/xyuan/Software_Hopper/pdslin_0.0
#HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lpdslin

ifeq ($(COM), 1)
      SCOREC_DIR = /global/project/projectdirs/mp288/cori/scorec/Dec2015
      SCOREC_CORE = -lcrv -ldsp -lph -lsam -lspr -lma -lparma -lapf_zoltan -lmds -lapf -llion -lmth -lgmi -lpcu
      # SCOREC_CORE = -lcrv -ldsp -lph -lsize -lsam -lspr -lma -lparma -lapf_zoltan -lmds -lapf -llion -lmth -lgmi -lpcu
      PETSC_DIR = /global/project/projectdirs/mp288/cori/petsc-3.5.4
      SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec_complex -Wl,--end-group
      PETSC_ARCH = complex-intel-mpich7.3
      HYPRE_LIB = 
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(HYPRE_LIB) \
       -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_4.3 -lsuperlu_dist_3.3 \
       -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -lnetcdf -ldl -lstdc++
else
      SCOREC_DIR = /global/project/projectdirs/mp288/cori/scorec/Dec2015
      SCOREC_LIBS=-L$(SCOREC_DIR)/lib -Wl,--start-group $(SCOREC_CORE) -lm3dc1_scorec -Wl,--end-group
      SCOREC_CORE = -lcrv -ldsp -lph -lsam -lspr -lma -lparma -lapf_zoltan -lmds -lapf -llion -lmth -lgmi -lpcu
      # SCOREC_CORE = -lcrv -ldsp -lph -lsize -lsam -lspr -lma -lparma -lapf_zoltan -lmds -lapf -llion -lmth -lgmi -lpcu
#o    PETSC_DIR = /global/project/projectdirs/mp288/cori/petsc-3.5.4
#o    PETSC_ARCH = real-intel-mpich7.3
      PETSC_DIR = /global/homes/j/jinchen/project/PETSC/master.noomp-nostrumpack
      PETSC_ARCH = next-noomp-nostrumpack
#n2      PETSC_DIR = /global/homes/j/jinchen/project/PETSC/master.omp-nostrumpack
#n2      PETSC_ARCH = next-omp-nostrumpack
#n3      PETSC_DIR = /global/homes/j/jinchen/project/PETSC/master.omp-strumpack
#n3      PETSC_ARCH = next-omp-strumpack
      HYPRE_LIB = -lHYPRE
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(HYPRE_LIB) \
       -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu -lsuperlu_dist -lstrumpack_sparse \
       -lflapack -lfblas -lparmetis -lmetis -lpthread -lssl -lcrypto -ldl -lstdc++ \
       -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis -lscotch -lscotcherr -lscotcherrexit

      OPTS := $(OPTS) -DNEXTPetscDEV
endif

ifeq ($(USEADIOS), 1)
  OPTS := $(OPTS) -DUSEADIOS

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
#ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
ADIOS_DIR=/usr/common/usg/adios/1.4.1
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L/usr/common/usg/minixml/2.7/lib -lm -lmxml \
             -L/usr/lib64/ -llustreapi
endif

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS #-DUSEHYBRID -DCJ_MATRIX_DUMP

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
           $(FFTW_INCLUDE_OPTS) \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include
#
LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc $(PETSC_EXTERNAL_LIB_BASIC) \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	$(FFTW_POST_LINK_OPTS) -lfftw3 \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB)
#        $(HYBRID_LIBS) \

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
        -Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic
  FOPTS  := $(FOPTS)  -g -dynamic
  CCOPTS := $(CCOPTS) -g -dynamic
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -dynamic
  FOPTS  := $(FOPTS)  -O3
  CCOPTS := $(CCOPTS) -O3
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
