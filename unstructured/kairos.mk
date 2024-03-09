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
 
OPTS := $(OPTS) -DPETSC_VERSION=39 -DUSEBLAS

SCOREC_BASE_DIR=/opt/nfri/lib/m3dc1/LIB/core-trunk/KAIROS
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

PETSC_DIR=/home/jchen/LIB/petsc-3.19.5
ifeq ($(COM), 1)
  PETSC_ARCH=kairos-intel-craympich-cplx
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/home/jchen/LIB/petsc-3.19.5/kairos-intel-craympich-cplx/lib -L/home/jchen/LIB/petsc-3.19.5/kairos-intel-craympich-cplx/lib -L/opt/cray/pe/fftw/3.3.8.4/x86_64/lib -L/opt/cray/pe/libsci/19.06.1/INTEL/16.0/x86_64/lib -L/opt/cray/pe/hdf5-parallel/1.10.5.2/intel/19.0/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.6.3.2/intel/19.0/lib -L/opt/cray/dmapp/default/lib64 -L/opt/cray/pe/mpt/7.7.11/gni/mpich-intel/16.0/lib -L/opt/nfri/lib/gsl/2.6/GNU/74/lib -L/opt/cray/rca/2.2.20-7.0.1.1_4.28__g8e3fb5b.ari/lib64 -L/opt/cray/pe/atp/2.1.3/libApp -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/libfabric/lib -L/opt/intel/compilers_and_libraries_2020.0.166/linux/ipp/lib/intel64 -L/opt/intel/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.0.166/linux/tbb/lib/intel64/gcc4.8 -L/opt/intel/compilers_and_libraries_2020.0.166/linux/daal/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -L/usr/lib64/gcc/x86_64-suse-linux/7 -L/usr/x86_64-suse-linux/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lrca -lfftw3_mpi -lfftw3_threads -lfftw3 -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lhugetlbfs -limf -lAtpSigHandler -lAtpSigHCommData -lhdf5_hl_parallel -lhdf5_parallel -lnetcdf_c++4 -lnetcdf -lsci_intel_mpi -lsci_intel -lmpich_intel -lmpichcxx_intel -lstdc++ -lifcoremt -lifport -lpthread -lsvml -lirng -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lgsl -lgslcblas -lm -lrca -lfftw3_mpi -lfftw3_threads -lfftw3 -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lhugetlbfs -limf -lAtpSigHandler -lAtpSigHCommData -lhdf5_hl_parallel -lhdf5_parallel -lnetcdf_c++4 -lnetcdf -lsci_intel_mpi -lsci_intel -lmpich_intel -lmpichcxx_intel -lstdc++ -lifcoremt -lifport -lpthread -lsvml -lirng -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
else
  PETSC_ARCH=kairos-intel-craympich-real
PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/home/jchen/LIB/petsc-3.19.5/kairos-intel-craympich-real/lib -L/home/jchen/LIB/petsc-3.19.5/kairos-intel-craympich-real/lib -L/opt/cray/pe/fftw/3.3.8.4/x86_64/lib -L/opt/cray/pe/libsci/19.06.1/INTEL/16.0/x86_64/lib -L/opt/cray/pe/hdf5-parallel/1.10.5.2/intel/19.0/lib -L/opt/cray/pe/netcdf-hdf5parallel/4.6.3.2/intel/19.0/lib -L/opt/cray/dmapp/default/lib64 -L/opt/cray/pe/mpt/7.7.11/gni/mpich-intel/16.0/lib -L/opt/nfri/lib/gsl/2.6/GNU/74/lib -L/opt/cray/rca/2.2.20-7.0.1.1_4.28__g8e3fb5b.ari/lib64 -L/opt/cray/pe/atp/2.1.3/libApp -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/libfabric/lib -L/opt/intel/compilers_and_libraries_2020.0.166/linux/ipp/lib/intel64 -L/opt/intel/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.0.166/linux/tbb/lib/intel64/gcc4.8 -L/opt/intel/compilers_and_libraries_2020.0.166/linux/daal/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -L/usr/lib64/gcc/x86_64-suse-linux/7 -L/usr/x86_64-suse-linux/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lzoltan -lparmetis -lmetis -lgsl -lgslcblas -lm -lrca -lfftw3_mpi -lfftw3_threads -lfftw3 -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lhugetlbfs -limf -lAtpSigHandler -lAtpSigHCommData -lhdf5_hl_parallel -lhdf5_parallel -lnetcdf_c++4 -lnetcdf -lsci_intel_mpi -lsci_intel -lmpich_intel -lmpichcxx_intel -lstdc++ -lifcoremt -lifport -lpthread -lsvml -lirng -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lgsl -lgslcblas -lm -lrca -lfftw3_mpi -lfftw3_threads -lfftw3 -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lhugetlbfs -limf -lAtpSigHandler -lAtpSigHCommData -lhdf5_hl_parallel -lhdf5_parallel -lnetcdf_c++4 -lnetcdf -lsci_intel_mpi -lsci_intel -lmpich_intel -lmpichcxx_intel -lstdc++ -lifcoremt -lifport -lpthread -lsvml -lirng -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
endif

MKL_LIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

ZOLTAN_LIB=-L$(SCOREC_BASE_DIR)/lib -lzoltan
GSL_DIR=$(PACKAGE_DIR)

INCLUDE := $(INCLUDE) -I$(SCOREC_BASE_DIR)/include \
	   -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include \
#   -I$(HDF5_DIR)/include \
#   -I$(GSL_DIR)/include 

LIBS := $(LIBS) \
        -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
        $(SCOREC_LIB) \
        $(ZOLTAN_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \
	$(MKL_LIB)\

#-L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 \
#-L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz \
#-L$(GSL_DIR)/lib -lgsl -lgslcblas \

ifeq ($(ST), 1)
LIBS += -Wl,--start-group -L$(NETCDF_DIR)/lib -Wl,-rpath,$(NETCDF_DIR)/lib -lnetcdf -lnetcdff -lz -Wl,--end-group
INCLUDE += -I$(NETCDF_DIR)/include
endif

FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS)

CCOPTS  = -c $(OPTS)

# Optimization flags
ifeq ($(VTUNE), 1)
  LDOPTS := $(LDOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
  FOPTS  := $(FOPTS)  -g -dynamic -debug inline-debug-info -parallel-source-info=2
  CCOPTS := $(CCOPTS) -g -dynamic -debug inline-debug-info -parallel-source-info=2
endif

# Optimization flags
ifeq ($(OPT), 1)
# LDOPTS := $(LDOPTS) -static -qopt-report
  LDOPTS := $(LDOPTS) -qopt-report
  FOPTS  := $(FOPTS)  -qopt-report
  CCOPTS := $(CCOPTS) -qopt-report
else
  FOPTS := $(FOPTS) -g -Mbounds -check noarg_temp_created -fpe0 -warn -traceback -debug extended
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
