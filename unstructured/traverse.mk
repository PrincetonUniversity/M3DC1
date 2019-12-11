  CPP = mpic++
  CC = mpicc
  F90 = mpifort
  F77 = mpifort
  LOADER = mpifort

OPTS := $(OPTS) -DPETSC_VERSION=990 -DUSEBLAS #-DNEWSOLVERDEVELOPMENT

PETSCVER=petsc
PETSC_VER=petsc

SCOREC_BASE_DIR=/projects/M3DC1/PETSC/petsc/traverse-pgi-openmpi-199-gpu-cuda-master/scorec199/

SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ifdef SCORECVER
    SCOREC_DIR=/projects/M3DC1/PETSC/petsc/traverse-pgi-openmpi-199-gpu-cuda-master/pumi199/
else	  
    SCOREC_DIR=/projects/M3DC1/PETSC/petsc/traverse-pgi-openmpi-199-gpu-cuda-master/pumi199/
endif
		
#zoltan is not available		
ZOLTAN_LIB=
		
ifeq ($(COM), 1)
    M3DC1_SCOREC_LIB = m3dc1_scorec_complex
else
    M3DC1_SCOREC_LIB = m3dc1_scorec
endif

SCOREC_LIBS= -L$(SCOREC_BASE_DIR)/lib -l$(M3DC1_SCOREC_LIB) -Wl,--start-group,-rpath,$(SCOREC_DIR)/lib -L$(SCOREC_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

PETSC_DIR=/home/jinchen/project/PETSC/petsc
ifeq ($(COM), 1)
  #PETSC_ARCH=cplx-pgi-cuda-mpi10.3.0
  PETSC_ARCH=traverse-pgi-openmpi-2019-2cplx
else
  ##PETSC_ARCH=real-pgi-cuda-mpi10.3.0
  #PETSC_ARCH=traverse-pgi-openmpi-2019-2
  #PETSC_ARCH=traverse-pgi-openmpi-2019-2-hdf
  #PETSC_ARCH=traverse-pgi-openmpi-2019-gpu
  #PETSC_ARCH=traverse-pgi-openmpi-2019-gpu-cuda
  #PETSC_ARCH=traverse-pgi-openmpi-del-gpu-cuda
   PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-master
  #PETSC_ARCH=traverse-pgi-openmpi-195-gpu-cuda-master
  #PETSC_ARCH=traverse-pgi-openmpi-195-master
  #PETSC_ARCH=traverse-pgi-openmpi-199-master
  #PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda
  #PETSC_ARCH=traverse-pgi-openmpi-195-gpu-cuda
endif

ifeq ($(PAR), 1)
  OPTS := $(OPTS) -DUSEPARTICLES
endif
		
#PETSC_WITH_EXTERNAL_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath, $(PETSC_DIR)/$(PETSC_ARCH)/lib $(PGI)/linuxpower/19.9/lib/pgi.ld -L$(PGI)/linuxpower/19.9/lib -L/usr/lib/gcc/ppc64le-redhat-linux/4.8.5 -Wl,-rpath, $(PGI)/linuxpower/19.9/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lsuperlu_dist -lparmetis -lmetis -lsuperlu -lscalapack -lflapack -lfblas -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lmpi_usempif08 -lmpi_mpifh -lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -lrt -lpgatm -lstdc++ -lrt -lm -lpthread -lz -L$(PGI)/linuxpower/19.9/lib -L/usr/lib/gcc/ppc64le-redhat-linux/4.8.5 -ldl -lpthread -lmpi -Wl,-rpath, $(PGI)/linuxpower/19.9/lib -latomic -lpgkomp -lomptarget -lpgmath -lpgc -lmass_simdp9 -lmassvp9 -lmassp9 -lm -lgcc_s -ldl

PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib /opt/pgi/19.5/linuxpower/19.5/lib/pgi.ld -L/usr/local/pgi/lib64 -L/usr/local/pgi/lib64/openmpi -L/usr/local/openmpi/4.0.2rc1/pgi195/ppc64le/lib64 -L/opt/pgi/19.5/linuxpower/19.5/lib -L/usr/lib/gcc/ppc64le-redhat-linux/4.8.5 -Wl,-rpath,/usr/local/pgi/lib64 -Wl,-rpath,/usr/local/pgi/lib64/openmpi -Wl,-rpath,/usr/local/openmpi/4.0.2rc1/pgi195/ppc64le/lib64 -Wl,-rpath,/opt/pgi/19.5/linuxpower/19.5/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lparmetis -lmetis -lml -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lpgf90rtl -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -latomic -lpgkomp -lomp -lomptarget -lpthread -lpgmath -lpgc -lrt -lmass_simdp9 -lmassvp9 -lmassp9 -lm -lgcc_s -lstdc++ -ldl

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

