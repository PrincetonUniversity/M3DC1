ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CC     = tau_cc.sh $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
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
HYBRID_HOME = /p/swim/jchen/hybrid.test

ifeq ($(USESCOREC), 1)
    ifndef SCORECDIR
      SCORECDIR = /project/projectdirs/mp288/lib/hopper2/scorec/install-Opt
    endif

SCOREC_LIBS =  \
        $(SCORECDIR)/lib/libFUSIONAPP.a \
        $(SCORECDIR)/lib/libSOLVER.a \
        $(SCORECDIR)/lib/libMESHADAPTMAP.a \
        $(SCORECDIR)/lib/libSOLTRANSFER.a \
        $(SCORECDIR)/lib/libSOLVER.a \
        $(SCORECDIR)/lib/libFEMANALYSIS.a \
        $(SCORECDIR)/lib/libASSEMBLER.a \
        $(SCORECDIR)/lib/libMeshAdapt.a \
        $(SCORECDIR)/lib/libDISCERRORESTIM.a \
        $(SCORECDIR)/lib/libASF.a \
        $(SCORECDIR)/lib/libSCORECModel.a \
        $(SCORECDIR)/lib/libmeshModel.a \
        $(SCORECDIR)/lib/libFMDB.a \
        $(SCORECDIR)/lib/libSCORECUtil.a \
        $(SCORECDIR)/lib/libipcomman.a \
        $(SCORECDIR)/lib/libzoltan.a \
        $(SCORECDIR)/lib/libSPARSKIT.a \
        $(SCORECDIR)/lib/libSCORECModel.a \
        $(SCORECDIR)/lib/libmeshModel.a \
        $(SCORECDIR)/lib/libSCORECUtil.a

  INCLUDE := $(INCLUDE) -I$(SCORECDIR)/include
  LIBS := $(LIBS) $(SCOREC_LIBS) -lC -lstd

  PARMETIS = -lparmetis -lmetis

else
#  OPTS := $(OPTS) -DPetscDEV
endif   # on USESCOREC

OPTS := $(OPTS) -DPetscDEV
PETSC_DIR = /project/projectdirs/mp288/lib/hopper2/petsc/petsc-dev-SUPERLU-HYPRE-MUMPS/petsc-dev-060711/petsc-dev
PETSC_ARCH = arch-cray-xt6-pkgs-opt
SUPERLU_DIST = -lsuperlu_dist_2.5
HYPRE = -lHYPRE
MUMPS = -ldmumps -lmumps_common -lpord

INCLUDE := $(INCLUDE) $(HDF5_INCLUDE_OPTS) \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include
LIBS := $(LIBS) $(HDF5_POST_LINK_OPTS) -lhdf5_fortran -lhdf5 \
	-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
	$(SUPERLU_DIST) $(HYPRE) $(MUMPS) $(PARMETIS) -ldl

FOPTS = -c -Mr8 -Mpreprocess -Minform=warn $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c -O $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) -fastsse -Mipa=fast,inline
  FOPTS  := $(FOPTS) -fastsse -Mipa=fast,inline
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -Mbounds
  CCOPTS := $(CCOPTS)
endif


F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)


%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
