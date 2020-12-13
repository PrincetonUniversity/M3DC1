  CPP = mpic++
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort

OPTS := $(OPTS) -DUSEADIOS -DPETSC_VERSION=37 -DUSEBLAS

SCOREC_BASE_DIR=/gpfs/wolf/gen127/proj-shared/scorec/pgi18.10-cuda10.1-mpi10.3
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

SCOREC_DIR=$(SCOREC_BASE_DIR)

ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

ZOLTAN_LIB=

SCOREC_LIBS= -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -l$(M3DC1_SCOREC_LIB) -Wl,--end-group

PETSC_DIR=/gpfs/wolf/gen127/proj-shared/petsc/petsc-3.7.6
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-pgi18.10-cuda-mpi10.3.0
else
  PETSC_ARCH=real-pgi18.10-cuda-mpi10.3.0
endif

PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath, $(PETSC_DIR)/$(PETSC_ARCH)/lib $(OLCF_PGI_ROOT)/linuxpower/18.10/lib/pgi.ld –L$(MPI_ROOT)/lib -L$(OLCF_PGI_ROOT)/linuxpower/18.10/lib -L/usr/lib/gcc/ppc64le-redhat-linux/4.8.5 -Wl,-rpath, $(OLCF_PGI_ROOT)/linuxpower/18.10/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_dist -lparmetis -lmetis -lsuperlu -lscalapack -lflapack -lfblas -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lmpi_ibm_usempi -lmpi_ibm_mpifh -lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -lrt -lpgatm -lstdc++ -lrt -lm -lpthread -lz -L$(MPI_ROOT)/lib -L$(OLCF_PGI_ROOT)/linuxpower/18.10/lib -L/usr/lib/gcc/ppc64le-redhat-linux/4.8.5 -ldl -lpthread -lmpiprofilesupport -lmpi_ibm -Wl,-rpath, $(OLCF_PGI_ROOT)/linuxpower/18.10/lib -latomic -lpgkomp -lpgompstub -lomptarget -lpgmath -lpgc -lmass_simdp9 -lmassvp9 -lmassp9 -lm -lgcc_s -ldl

#only define them if adios-1.3 is used; otherwise use hopper default
INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(GSL_DIR)/include # \
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
LIBS := $(LIBS) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	-L$(GSL_DIR)/lib -lgsl -lhugetlbfs 

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
# FIXME 
ifeq ($(OPT), 1)
  LDOPTS := $(LDOPTS) #-static -qopt-report
  FOPTS  := $(FOPTS)  #-qopt-report
  CCOPTS := $(CCOPTS) #-qopt-report
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
