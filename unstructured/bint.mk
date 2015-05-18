ifeq ($(TAU), 1)
  TAU_OPTIONS = -optCPPOpts=-DUSETAU -optVerbose -optPreProcess -optMpi -optTauSelectFile=select.tau
  CPP    = tau_cxx.sh $(TAU_OPTIONS)
  CC     = tau_cc.sh  $(TAU_OPTIONS)
  F90    = tau_f90.sh $(TAU_OPTIONS)
  F77    = tau_f90.sh $(TAU_OPTIONS)
  LOADER = tau_f90.sh $(TAU_OPTIONS)
else
  CPP = mpiicpc -DMPICH_IGNORE_CXX_SEEK -mmic
  CC = mpiicc -mmic
  F90 = mpiifort -mmic -align array64byte
  F77 = mpiifort -mmic -align array64byte
  LOADER = mpiifort -mmic -align array64byte
endif

ifeq ($(HPCTK), 1)
  OPTS := $(OPTS) -gopt
  LOADER := hpclink $(LOADER)
endif

    SCORECDIR = /chos/global/project/projectdirs/mp288/babbage/scorec/May2015
ifeq ($(COM), 1)
    SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib -lm3dc1_scorec_complex -lpcu -lgmi -lapf -lmds -lspr -lapf_zoltan -lparma -lma -lph -Wl,--end-group 
    PETSC_DIR= /usr/common/usg/petsc/3.5.2-fee0b69/complex
    PETSC_ARCH= knc_201501
    PETSCLIB= -Wl,--start-group,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lmumps_common -lpetsc -lpord -lscalapack -lsmumps -lsuperlu_4.3 -lsuperlu_dist_3.3 -lzmumps -lparmetis -lmetis -Wl,--end-group

else
    SCORECLIB= -Wl,--start-group,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib -lm3dc1_scorec -lpcu -lgmi -lapf -lmds -lspr -lapf_zoltan -lparma -lma -lph -Wl,--end-group
    PETSC_DIR= /usr/common/usg/petsc/3.5.2-fee0b69/real
    PETSC_ARCH= knc_201501
    PETSCLIB = -Wl,--start-group,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lcmumps -ldmumps -lmumps_common -lpetsc -lpord -lscalapack -lsmumps -lsuperlu_4.3 -lsuperlu_dist_3.3 -lzmumps -lparmetis -lmetis -Wl,--end-group
endif

  INCLUDE := -I$(SCORECDIR)/include
  LIBS := $(LIBS) \
	$(SCORECLIB) \
        -L$(ZOLTAN_DIR)/lib -lzoltan \
	$(PETSCLIB) \
        -L$(MKLROOT)/lib/mic/ -lmkl_blacs_intelmpi_lp64 -lmkl_blas95_lp64 -lmkl_cdft_core -lmkl_core -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_sequential \
        -lstdc++

AUX = d1mach.o i1mach.o r1mach.o fdump.o dbesj0.o dbesj1.o

OPTS := $(OPTS) -DPetscDEV -DKSPITS -DxUSEADIOS #-DUSEHYBRID -DCJ_MATRIX_DUMP

#only define them if adios-1.3 is used; otherwise use hopper default
#ADIOS_DIR=/usr/common/usg/adios/1.4.1
#ADIOS_FLIB = -L${ADIOS_DIR}/lib -ladiosf_v1 -ladiosreadf_v1 \
#             -L/usr/common/usg/minixml/2.7/lib -lm -lmxml \
#             -L/usr/lib64/ -llustreapi

INCLUDE := $(INCLUDE) -I$(HDF5_PAR_DIR)/include -I$(FFTW_DIR)/include \
	-I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	-I$(GSL_DIR)/include

LIBS := $(LIBS) \
	-Wl,--start-group,-rpath,$(HDF5_PAR_DIR)/lib -L$(HDF5_PAR_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lstdc++ -ldl -Wl,--end-group \
	-L$(FFTW_DIR)/lib -lfftw3 \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) \
	-Dglobalinsertval=insertval -Dglobalentdofs=entdofs
CCOPTS  = -c -O $(OPTS)

# Optimization flags
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS)
  FOPTS  := $(FOPTS)  -O0
  CCOPTS := $(CCOPTS)
else
  FOPTS := $(FOPTS) -g -Mbounds
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
