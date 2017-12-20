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

SCORECDIR = /global/project/projectdirs/mp288/hopper/scorec/May2015
ifeq ($(COM),1)
  SCORECLIB=-lapf -lgmi -lm3dc1_scorec_complex -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr
else
  SCORECLIB=-lapf -lgmi -lm3dc1_scorec -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr
endif
SCOREC_LIBS =-L$(SCORECDIR)/lib -Wl,--start-group $(SCORECLIB) -Wl,--end-group
INCLUDE := $(INCLUDE) -I$(SCORECDIR)/include
LIBS := $(LIBS) $(SCOREC_LIBS) #-lC -lstd

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DxUSEADIOS -DKSPITS -DNO_STOP_MESSAGE=1 -DPetscOLD #-DUSEHYBRID -DCJ_MATRIX_DUMP
#PETSC_DIR = /project/projectdirs/mp288/lib/hopper2/petsc/petsc-dev-SUPERLU-HYPRE-MUMPS/petsc-dev-060711/petsc-dev
#PETSC_ARCH = arch-cray-xt6-pkgs-opt
#SUPERLU_DIST = -lsuperlu_dist_2.5
#HYPRE = -lHYPRE
#MUMPS = -ldmumps -lmumps_common -lpord

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/global/homes/p/pnorbert/adios/hopper
ADIOS_DIR=/global/homes/p/pnorbert/adios/1.3.1/hopper/pgi/
ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf -L/global/homes/p/pnorbert/mxml/mxml.hopper/lib -lm -lmxml -llustreapi -pgcpplibs

INCLUDE := $(INCLUDE) -I$(HDF5_INCLUDE_OPTS) $(FFTW_INCLUDE_OPTS) \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
        -I$(HYBRID_HOME)/include -I$(GSL_DIR)/include

LIBS := $(LIBS) $(HDF5_POST_LINK_OPTS) -lhdf5_fortran -lhdf5 \
	$(FFTW_POST_LINK_OPTS) -lfftw3 \
	$(HYPRE) $(MUMPS) $(PARMETIS) -ldl \
        $(HYBRID_LIBS) \
        $(ADIOS_FLIB) \
	-L$(GSL_DIR)/lib -lgsl

#FOPTS = -c -Mr8 -Mpreprocess -Minform=warn $(OPTS) \

FOPTS = -c -s real64 -e Fm -rm $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c -O $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
#  LDOPTS := $(LDOPTS) -fastsse -Mipa=fast,inline
#  FOPTS  := $(FOPTS) -fastsse -Mipa=fast,inline
  LDOPTS := $(LDOPTS) 
  FOPTS  := $(FOPTS) 
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -Mbounds -Mchkptr -traceback
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
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
