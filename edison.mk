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
      SCORECDIR = /global/project/projectdirs/mp288/edison/scorec/complex
      PETSC_DIR = /global/project/projectdirs/mp288/edison/petsc-complex/petsc-3.5.2
      PETSC_ARCH = arch-xc30-opt
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,/global/project/projectdirs/mp288/edison/petsc-complex/petsc-3.5.2/arch-xc30-opt/lib -L/global/project/projectdirs/mp288/edison/petsc-complex/petsc-3.5.2/arch-xc30-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_4.3 -lsuperlu_dist_3.3 -lparmetis -lmetis -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
else
      SCORECDIR = /global/project/projectdirs/mp288/edison/scorec/real
      PETSC_DIR =/global/homes/h/hzhang/petsc/
      PETSC_ARCH =arch-xc30-opt
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,/global/u2/h/hzhang/petsc/arch-xc30-opt/lib -L/global/u2/h/hzhang/petsc/arch-xc30-opt/lib -lHYPRE -lsuperlu_4.3 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lml -lsuperlu_dist_3.3 -lparmetis -lmetis -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
endif

SCORECLIB="-lapf -lgmi -lm3dc1_scorec -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr"
SCOREC_LIBS =-L$(SCORECDIR)/lib -Wl,--start-group -lm3dc1_scorec $(SCORECLIB) -Wl,--end-group -lzoltan
INCLUDE := $(INCLUDE) -I$(SCORECDIR)/include
LIBS := $(LIBS) \
        -L$(SCORECDIR)/lib $(SCOREC_LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc $(PETSC_EXTERNAL_LIB_BASIC)  \
        -lstdc++

ifeq ($(USEADIOS), 1)
  OPTS := $(OPTS) -DUSEADIOS
endif

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS #-DUSEHYBRID -DCJ_MATRIX_DUMP

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
#ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs
ADIOS_DIR=/usr/common/usg/adios/1.4.1
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L/usr/common/usg/minixml/2.7/lib -lm -lmxml \
             -L/usr/lib64/ -llustreapi

INCLUDE := $(INCLUDE) -I$(HDF5_DIR)/include $(FFTW_INCLUDE_OPTS) \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	-I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include

LIBS := $(LIBS) -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz \
	$(FFTW_POST_LINK_OPTS) -lfftw3 \
	$(HYPRE) $(MUMPS) $(PARMETIS) -ldl \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB)
#        $(HYBRID_LIBS) \


FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c -O $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)  -O0
  CCOPTS := $(CCOPTS)
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)  
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
