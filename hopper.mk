FOPTS = -c -Mr8 -Mpreprocess -fastsse -Minform=warn -Mipa=fast,inline $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs \
	-DPetscDEV
CCOPTS  = -c -O -DPetscDEV

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
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS)
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/hybrid.test

ifeq ($(USESCOREC), 1)
    ifndef SCORECDIR
      SCORECDIR = /project/projectdirs/mp288/lib/hopper2/install/03032011
    endif

    SCOREC_LIBS = \
        -L$(SCORECDIR)/lib \
	-Wl,-rpath,$(SCORECDIR)/lib \
        -lsuperlu_dist_2.5 \
	-lFUSIONAPP$(SCORECOPT) \
        -lSOLVER$(SCORECOPT) \
        -lMESHADAPTMAP$(SCORECOPT) \
        -lSOLTRANSFER$(SCORECOPT) \
        -lSOLVER$(SCORECOPT) \
        -lFEMANALYSIS$(SCORECOPT) \
        -lASSEMBLER$(SCORECOPT) \
        -lMeshAdapt$(SCORECOPT) \
        -lDISCERRORESTIM$(SCORECOPT) \
        -lASF$(SCORECOPT) \
        -lSCORECModel$(SCORECOPT) \
        -lmeshModel$(SCORECOPT) \
        -lFMDB$(SCORECOPT) \
        -lSCORECUtil$(SCORECOPT) \
	-lipcomman$(SCORECOPT) \
	-lzoltan \
        -lSPARSKIT$(SCORECOPT) \
        -lSCORECModel$(SCORECOPT) \
        -lmeshModel$(SCORECOPT) \
        -lSCORECUtil$(SCORECOPT) \
	-lparmetis -lmetis \
	-lsuperlu_dist_2.5 \
	-lC -lstd

  PETSC_DIR = /project/projectdirs/mp288/lib/hopper2/petsc/petsc-dev
  PETSC_ARCH = arch-cray-xt5-opt

  INCLUDE := $(INCLUDE) -I$(SCORECDIR)/include
  LIBS := $(LIBS) $(SCOREC_LIBS)

endif   # on USESCOREC

INCLUDE := $(INCLUDE) $(HDF5_INCLUDE_OPTS) \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIBS := $(LIBS) $(HDF5_POST_LINK_OPTS) -lhdf5_fortran -lhdf5 \
	-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc


%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
