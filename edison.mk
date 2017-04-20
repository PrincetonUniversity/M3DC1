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

ifeq ($(COM), 1)
  SCOREC_DIR = /global/project/projectdirs/mp288/edison/scorec/Apr2017-mpich7.2.5
  M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
  SCOREC_DIR = /global/project/projectdirs/mp288/edison/scorec/Mar2017-mpich7.4.1
  ifeq ($(TRILINOS), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_trilinos
  else
    M3DC1_SCOREC_LIB = m3dc1_scorec
  endif
endif

  SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
               -lpumi -lcrv -lph -lsam -lspr -lma \
               -lapf_zoltan -lparma -lmds -lapf -llion -lmth -lgmi -lpcu -l$(M3DC1_SCOREC_LIB) \
               -Wl,--end-group

ifeq ($(TRILINOS),1)
TRILINOS_LIBS = -Wl,--start-group,-rpath,$(CRAY_TRILINOS_PREFIX_DIR)/lib -L$(CRAY_TRILINOS_PREFIX_DIR)/lib \
                -lamesos -ltpetra -lkokkosnodeapi -ltpi -laztecoo -lepetra \
                -lsacado -lteuchosparameterlist -lteuchoscomm_intel -lteuchoscore -lteuchosnumerics \
                -lteuchosremainder -Wl,--end-group
else
TRILINOS_LIBS =
endif

ifeq ($(COM), 1)
      PETSC_DIR =/global/project/projectdirs/mp288/edison/petsc-3.5.4-complex
      PETSC_ARCH =cray-mpich-7.2
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_4.3 -lsuperlu_dist_3.3 -lzoltan -lparmetis -lmetis -lsci_intel_mpi_mp -lsci_intel_mp -liomp5 -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
      PETSC_LIB =  -lpetsc
      ZOLTAN_LIB =
else
      PETSC_DIR =/opt/cray/petsc/3.7.4.0/real/INTEL/15.0/sandybridge
      PETSC_ARCH =
ifeq ($(TRILINOS), 1)
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L/opt/cray/tpsl/1.4.3/INTEL/140/sandybridge/lib -lparmetis -lmetis -lpthread -lssl -lcrypto -Wl,-rpath,/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -L/opt/cray/hdf5-parallel/1.8.11/intel/130/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
      PETSC_LIB =
else
      PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(CRAY_TPSL_DIR)/INTEL/150/sandybridge/lib -lHYPRE -lsuperlu -lcmumps -ldmumps -lesmumps -lsmumps -lzmumps -lmumps_common -lptesmumps -lpord -lsuperlu_dist -lparmetis -lmetis -lptscotch -lscotch -lptscotcherr -lscotcherr -lsci_intel_mpi_mp -lsci_intel_mp -liomp5 -lsundials_cvode -lsundials_cvodes -lsundials_ida -lsundials_idas -lsundials_kinsol -lsundials_nvecparallel -lsundials_nvecserial -lpthread -lssl -lcrypto -Wl,-rpath,$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++
      PETSC_LIB = -lcraypetsc_intel_real
endif
      ZOLTAN_LIB = -L$(CRAY_TRILINOS_PREFIX_DIR)/lib -lzoltan
      OPTS := $(OPTS) -DNEXTPetscDEV
endif

ifeq ($(USEADIOS), 1)
  OPTS := $(OPTS) -DUSEADIOS
  ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
             -L/usr/common/usg/minixml/2.7/lib -lm -lmxml \
             -L/usr/lib64/ -llustreapi
else
  ADIOS_FLIB =
endif

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS #-DUSEHYBRID -DCJ_MATRIX_DUMP

INCLUDE := $(INCLUDE) -I$(HDF5_DIR)/include $(FFTW_INCLUDE_OPTS) \
        -I$(SCOREC_DIR)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include \
        -I$(PETSC_DIR)/include \
	-I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include

LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(TRILINOS_LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(PETSC_LIB) \
        $(PETSC_EXTERNAL_LIB_BASIC) -lstdc++ \
        -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz \
        $(FFTW_POST_LINK_OPTS) -lfftw3 \
        -L$(GSL_DIR)/lib -lgsl -lhugetlbfs \
        $(ADIOS_FLIB)

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c $(OPTS)

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)  -O3
  CCOPTS := $(CCOPTS) -O3
else
  FOPTS := $(FOPTS) -g -Mbounds -check all -check noarg-temp-created -fpe0 -warn -traceback -debug extended
  CCOPTS := $(CCOPTS)  
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -openmp -Wl,-z,muldefs 
  FOPTS  := $(FOPTS)  -openmp 
  CCOPTS := $(CCOPTS) -openmp 
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
