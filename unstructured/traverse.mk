  CPP = mpic++
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort

OPTS := $(OPTS) -DPETSC_VERSION=990 -DUSEBLAS

PETSCVER=petsc3.14.4
PETSC_VER=petsc-3.14.4

PETSC_DIR=/home/jinchen/project/PETSC/petsc
#PETSC_DIR=/home/jinchen/project/PETSC/petsc.rh7
ifeq ($(COM), 1)
   PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-cplx-master
  #PETSC_ARCH=traverse-pgi-openmpi-199-cplx-master
else
   PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-master
  #PETSC_ARCH=traverse-pgi-openmpi-199-master
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

#SCOREC_BASE_DIR=/projects/M3DC1/scorec/pgi19.9-openmpi4.0.2/$(PETSCVER)
SCOREC_BASE_DIR=/projects/M3DC1/scorec/pgi20.4-openmpi4.0.4/petsc3.13.4/202209
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

#zoltan is not available
ZOLTAN_LIB=

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif
		
ifeq ($(COM), 1)
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/cuda-10.2/lib64 -L/usr/local/cuda-10.2/lib64 /opt/pgi/20.4/linuxpower/20.4/lib/pgi.ld -L/usr/local/pgi/lib64 -L/usr/local/pgi/lib64/openmpi -L/usr/local/openmpi/4.0.4/pgi204/ppc64le/lib64 -L/opt/pgi/20.4/linuxpower/20.4/lib -L/usr/lib/gcc/ppc64le-redhat-linux/8 -Wl,-rpath,/usr/local/pgi/lib64 -Wl,-rpath,/usr/local/pgi/lib64/openmpi -Wl,-rpath,/usr/local/openmpi/4.0.4/pgi204/ppc64le/lib64 -Wl,-rpath,/opt/pgi/20.4/linuxpower/20.4/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lparmetis -lmetis -lcufft -lcublas -lcudart -lcusparse -lcusolver -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -latomic -lnvomp -lpthread -lpgmath -lnvc -lrt -lmass_simdp9 -lmassvp9 -lmassp9 -lm -lgcc_s -lstdc++ -ldl
else
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/cuda-10.2/lib64 -L/usr/local/cuda-10.2/lib64 /opt/pgi/20.4/linuxpower/20.4/lib/pgi.ld -L/usr/local/pgi/lib64 -L/usr/local/pgi/lib64/openmpi -L/usr/local/openmpi/4.0.4/pgi204/ppc64le/lib64 -L/opt/pgi/20.4/linuxpower/20.4/lib -L/usr/lib/gcc/ppc64le-redhat-linux/8 -Wl,-rpath,/usr/local/pgi/lib64 -Wl,-rpath,/usr/local/pgi/lib64/openmpi -Wl,-rpath,/usr/local/openmpi/4.0.4/pgi204/ppc64le/lib64 -Wl,-rpath,/opt/pgi/20.4/linuxpower/20.4/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lml -lfftw3_mpi -lfftw3 -lflapack -lfblas -lparmetis -lmetis -lcufft -lcublas -lcudart -lcusparse -lcusolver -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -latomic -lnvomp -lpthread -lpgmath -lnvc -lrt -lmass_simdp9 -lmassvp9 -lmassp9 -lm -lgcc_s -lstdc++ -ldl
endif

#only define them if adios-1.3 is used; otherwise use hopper default
INCLUDE := $(INCLUDE) -I$(SCOREC_BASE_DIR)/include -I$(SCOREC_DIR)/include -I/home/jinchen/LIB/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
	   -I$(HDF5_ROOT)/include
#        -I$(HYBRID_HOME)/include
#           -I$(CRAY_TPSL_DIR)/INTEL/150/haswell/include \
#
CUDA_LIB:=-L/usr/local/cuda-10.1/lib64 -lcudart -lcusparse -lcusolver -lstdc++ -L/usr/local/cuda-10.0/lib64 -lcublas
LIBS := $(LIBS) \
        $(CUDA_LIB) \
        $(SCOREC_LIBS) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(HDF5_ROOT)/lib64 -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
	    -lgsl -lgslcblas -lfftw3_mpi -lfftw3

FOPTS = -c -r8 -Mpreprocess $(OPTS)

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
  LDOPTS := $(LDOPTS) -fast#-static -qopt-report
  FOPTS  := $(FOPTS)  -fast#-qopt-report
  CCOPTS := $(CCOPTS) -fast #-qopt-report
else
  FOPTS := $(FOPTS) -Mbounds -Minfo=all -Mchkfpstk -Mchkstk -Mdalign -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
  CCOPTS := $(CCOPTS) -g 
  LDOPTS := $(LDOPTS) -g
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp 
  FOPTS  := $(FOPTS)  -mp 
  CCOPTS := $(CCOPTS) -mp 
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

