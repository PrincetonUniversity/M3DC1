FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DUSEBLAS -DPETSC_VERSION=322
CCOPTS  = -c -DPETSC_VERSION=322
R8OPTS = -r8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O3
  CCOPTS := $(CCOPTS) -O3
else
  FOPTS := $(FOPTS) -g -check all -check noarg_temp_created -debug all -ftrapuv -traceback -fpe=all
  CCOPTS := $(CCOPTS) -g -check=uninit -debug all
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

CC = mpicc
CPP = mpicxx
F90 = mpif90
F77 = mpif77
LOADER = mpifort
LDOPTS := $(LDOPTS) -cxxlib
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
MPIVER=intel2020
PETSC_VER=petsc-3.22.4
PETSCVER=petsc3.22.4
PETSC_DIR=/fusion/projects/codes/m3dc1/omega_install/$(PETSC_VER)
ifeq ($(COM), 1)
  PETSC_ARCH=cplx-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/fusion/projects/codes/m3dc1/omega_install/petsc-3.22.4/cplx-intel2020/lib -L/fusion/projects/codes/m3dc1/omega_install/petsc-3.22.4/cplx-intel2020/lib -Wl,-rpath,/fusion/usc/opt/fftw/fftw-3.3.4-mpich-3.2_intel_2020/lib -L/fusion/usc/opt/fftw/fftw-3.3.4-mpich-3.2_intel_2020/lib -Wl,-rpath,/fusion/usc/opt/hdf5/hdf5-1.8.19_mpich-3.2_intel-2020/lib -L/fusion/usc/opt/hdf5/hdf5-1.8.19_mpich-3.2_intel-2020/lib -Wl,-rpath,/fusion/usc/opt/mpich/mpich-3.2/intel2020/lib -L/fusion/usc/opt/mpich/mpich-3.2/intel2020/lib -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/ipp/lib/intel64 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/ipp/lib/intel64 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/compiler/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/compiler/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/mkl/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/mkl/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64/gcc4.7 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64/gcc4.7 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/daal/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/daal/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64_lin/gcc4.4 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64_lin/gcc4.4 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/8 -L/usr/lib/gcc/x86_64-redhat-linux/8 -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lpthread -lparmetis -lmetis -lhdf5_hl -lhdf5 -lX11 -lmpifort -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lmpicxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lquadmath -lmpicxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
else
  PETSC_ARCH=real-$(MPIVER)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/fusion/projects/codes/m3dc1/omega_install/petsc-3.22.4/real-intel2020/lib -L/fusion/projects/codes/m3dc1/omega_install/petsc-3.22.4/real-intel2020/lib -Wl,-rpath,/fusion/usc/opt/fftw/fftw-3.3.4-mpich-3.2_intel_2020/lib -L/fusion/usc/opt/fftw/fftw-3.3.4-mpich-3.2_intel_2020/lib -Wl,-rpath,/fusion/usc/opt/hdf5/hdf5-1.8.19_mpich-3.2_intel-2020/lib -L/fusion/usc/opt/hdf5/hdf5-1.8.19_mpich-3.2_intel-2020/lib -Wl,-rpath,/fusion/usc/opt/mpich/mpich-3.2/intel2020/lib -L/fusion/usc/opt/mpich/mpich-3.2/intel2020/lib -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/ipp/lib/intel64 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/ipp/lib/intel64 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/compiler/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/compiler/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/mkl/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/mkl/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64/gcc4.7 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64/gcc4.7 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/daal/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/daal/lib/intel64_lin -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64_lin/gcc4.4 -L/fusion/usc/opt/intel2020/compilers_and_libraries/linux/tbb/lib/intel64_lin/gcc4.4 -Wl,-rpath,/fusion/usc/opt/intel2020/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -L/fusion/usc/opt/intel2020/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/8 -L/usr/lib/gcc/x86_64-redhat-linux/8 -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lpthread -lparmetis -lmetis -lhdf5_hl -lhdf5 -lX11 -lmpifort -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lmpicxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lquadmath -lmpicxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
endif

SCOREC_BASE_DIR=/fusion/projects/codes/m3dc1/omega_install/scorec/omega/$(MPIVER)/$(PETSCVER)/release
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ZOLTAN_LIB=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lzoltan

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

INCLUDE := $(INCLUDE) -I$(SCOREC_DIR)/include \
        -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(HDF5_DIR)/include \
	-I$(FFTW_DIR)/include

LIBS := $(LIBS) \
        -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
        $(SCOREC_LIBS) \
        $(PETSC_WITH_EXTERNAL_LIB) \
        -L$(HDF5_DIR)/lib -lhdf5_hl_fortran -lhdf5_fortran -lz\
        -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
        -lgsl -lgslcblas -lhugetlbfs \
        -lstdc++


ifeq ($(ST), 1)
  NETCDF_DIR=/fusion/projects/codes/m3dc1/scorec/netcdf/intel2020-mpich3.2
  INCLUDE += -I$(NETCDF_DIR)/include
  LIBS += -L$(NETCDF_DIR)/lib -lnetcdff  -L$(NETCDF_DIR)/lib64 -lnetcdf
endif

ifneq ($(USEADAS), 1)
  USEADAS = 1
  OPTS := $(OPTS) -DUSEADAS
endif
ADAS_LIB = -L$(ADASHOME)/lib/ifort -ladaslib


%.o : %.c
	$(CC)  $(CCOPTS) $(INCLUDE) $< -o $@

%.o : %.cpp
	$(CPP) $(CCOPTS) $(INCLUDE) $< -o $@

%.o: %.f
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.F
	$(F77) $(F77OPTS) $(INCLUDE) $< -o $@

%.o: %.f90
	$(F90) $(F90OPTS) $(INCLUDE) -fpic $< -o $@
