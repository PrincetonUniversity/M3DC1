#FOPTS = -c -r8 -implicitnone -fpp -warn all -DxPetscDEV -DPETSC_31 $(OPTS)
FOPTS = -c -r8 -implicitnone -fpp -warn all -DPetscDEV -DKSPITS -DxCJ_MATRIX_DUMP  $(OPTS) #-pg
CCOPTS  = -c -O -DPetscDEV -DKSPITS -DxCJ_MATRIX_DUMP -DxUSEHYBRID #-DPETSC_31 -pg

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -fast
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv
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
  LOADER = mpif90 -cxxlib #-pg
  FOPTS := $(FOPTS)
endif
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)


# define where you want to locate the mesh adapt libraries
HYBRID_HOME = /p/swim/jchen/pdslin_0.0
#HYBRID_HOME = /p/swim/jchen/hybrid.test
#HYBRID_HOME = /u/iyamazak/release/v2/hybrid.test
HYBRID_LIBS = -L$(HYBRID_HOME)/lib -lhsolver

INCLUDE = -I$(MPIHOME)/include \
	-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(GSLHOME)/include \
	-I$(HYBRID_HOME)/include

#	-L$(MUMPS_HOME)/lib -ldmumps -lmumps_common -lpord \
#	-L$(SCALAPACK_HOME) -lscalapack \
#	-L$(BLACS_HOME)/lib -lmpiblacsF77init -lmpiblacs -lmpiblacsCinit -lmpiblacs

LIBS = 	$(PETSC_LIBS) \
	-L$(FFTWHOME)/lib -lfftw3 -lfftw3_mpi -lfftw3_threads \
        -L$(CCHOME)/mkl/lib/em64t -lmkl -lmkl_lapack \
	-L$(CCHOME)/lib/intel64 -lguide \
	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c \
	-L$(GSLHOME)/lib -lgsl \
	-L/usr/X11R6/lib -lX11

ifeq ($(USESCOREC), 1)
  #SCORECDIR = /p/tsc/m3dc1/lib/develop.petsc3.Fan/build_greene/cmake
  SCORECDIR=/p/tsc/m3dc1/lib/develop.petsc3.Fan/build_greene/newpetsc
  SCORECLIB=-lapf -lapf_pumi -lpumi_util -lpumi_geom -lpcu -lpumi_geom_meshmodel -lpumi_mesh -lmeshadapt
  SCOREC_LIBS = -Wl,-rpath,$(SCORECDIR)/lib -L$(SCORECDIR)/lib -Wl,--start-group -lPPPLFusion $(SCORECLIB) -Wl,--end-group -lzoltan
 PETSC_EXTERNAL_LIB_BASIC = -Wl,-rpath,/p/tsc/m3dc1/lib/SCORECLib/petsc-3.5.1-greene/petsc-3.5.1/intel-ompi/lib -L/p/tsc/m3dc1/lib/SCORECLib/petsc-3.5.1-greene/petsc-3.5.1/intel-ompi/lib -lsuperlu_dist_3.3 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -Wl,-rpath,/usr/pppl/intel/11-pkgs/acml-4.3.0/ifort64/lib -L/usr/pppl/intel/11-pkgs/acml-4.3.0/ifort64/lib -lacml -lparmetis -lmetis -lX11 -lpthread -lssl -lcrypto -lhdf5_fortran -lhdf5_hl -lhdf5 -Wl,-rpath,/usr/pppl/intel/11-pkgs/openmpi-1.4.1/lib -L/usr/pppl/intel/11-pkgs/openmpi-1.4.1/lib -Wl,-rpath,/usr/pppl/intel/11.0/081/lib/intel64 -L/usr/pppl/intel/11.0/081/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -lmpi_f90 -lmpi_f77 -lifport -lifcore -lm -lmpi_cxx -lstdc++ -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -limf -lsvml -lipgo -ldecimal -lirc -lgcc_s -lpthread -lirc_s -ldl
  LIBS := $(SCOREC_LIBS) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib $(PETSC_EXTERNAL_LIB_BASIC) -lpetsc $(LIBS)
  INCLUDE := -I$(SCORECDIR)/include \
        $(INCLUDE)
endif   # on USESCOREC

%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
