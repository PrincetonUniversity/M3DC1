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
CORE_LIBS = -lapf -lgmi -lma -lparma -lph -lapf_zoltan -lmds -lpcu -lspr

ifeq ($(TRILINOS), 1)
  SCORECDIR = /global/project/projectdirs/mp288/seol/trilinos_dev
  SCORECLIB = $(CORE_LIBS) -lm3dc1_scorec_trilinos
  TRILINOS_LIBS = -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra -lsacado -lteuchosparameterlist -lteuchoscomm -lteuchoscore -lteuchosnumerics -lteuchosremainder
else
SCORECDIR = /global/project/projectdirs/mp288/edison/scorec/Jun2015
endif

ifeq ($(COM), 1)
  ifeq ($(TRILINOS), 1)
    SCORECLIB = $(CORE_LIBS) -lm3dc1_scorec_trilinos_complex
  else
    SCORECLIB = $(CORE_LIBS) -lm3dc1_scorec_complex
  endif
  PETSC_DIR =/opt/cray/petsc/3.5.2.1/complex/INTEL/140/sandybridge
  PETSC_ARCH =
  PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/opt/cray/tpsl/1.4.3/INTEL/140/sandybridge/lib -lsuperlu -lcmumps -ldmumps -lesmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu_dist -lparmetis -lmetis -lptscotch -lscotch -lptscotcherr -lscotcherr -lsci_intel_mpi_mp -lsci_intel_mp -liomp5 -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
else
  ifeq ($(TRILINOS), 1)
    SCORECLIB = $(CORE_LIBS) -lm3dc1_scorec_trilinos
  else
    SCORECLIB = $(CORE_LIBS) -lm3dc1_scorec
  endif
  PETSC_DIR =/opt/cray/petsc/3.5.2.1/real/INTEL/140/sandybridge
  PETSC_ARCH =
  PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/opt/cray/tpsl/1.4.3/INTEL/140/sandybridge/lib -lHYPRE -lsuperlu -lcmumps -ldmumps -lesmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu_dist -lparmetis -lmetis -lptscotch -lscotch -lptscotcherr -lscotcherr -lsci_intel_mpi_mp -lsci_intel_mp -liomp5 -lsundials_cvode -lsundials_cvodes -lsundials_ida -lsundials_idas -lsundials_kinsol -lsundials_nvecparallel -lsundials_nvecserial -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
endif

SCOREC_LIBS =-L$(SCORECDIR)/lib -Wl,--start-group $(SCORECLIB) -Wl,--end-group 
INCLUDE := $(INCLUDE) -I$(SCORECDIR)/include -I/opt/cray/tpsl/1.4.3/INTEL/140/sandybridge/include
LIBS := $(LIBS) -L$(SCORECDIR)/lib $(SCOREC_LIBS)

ifeq ($(COM), 1)
LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcraypetsc_intel_complex $(PETSC_EXTERNAL_LIB_BASIC)  
else
LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcraypetsc_intel_real $(PETSC_EXTERNAL_LIB_BASIC)
endif

ifeq ($(TRILINOS), 1)
INCLUDE := $(INCLUDE) -I$(CRAY_TRILINOS_PREFIX_DIR)/include
LIBS := $(LIBS) \
        -L$(CRAY_TRILINOS_PREFIX_DIR)/lib $(TRILINOS_LIBS) \
        -L$(PARALLEL_NETSDF_DIR/lib) -lpnetcdf \
        -L$(MKL_LIBDIR) -lmkl_lapack95_lp64 -lmkl_blas95_lp64 \
        -L$(BOOST_LIB) -lboost_mpi -lboost_exception
endif

LIBS := $(LIBS) -lstdc++

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
        -L$(CRAY_TRILINOS_PREFIX)/lib -lzoltan \
	$(FFTW_POST_LINK_OPTS) -lfftw3 \
	$(HYPRE) $(MUMPS) $(PARMETIS) -ldl \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
	$(ADIOS_FLIB)
#        $(HYBRID_LIBS) \


FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)  -O3
  CCOPTS := $(CCOPTS) -O3
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
