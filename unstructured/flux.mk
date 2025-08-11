FOPTS = -c -r8 -implicitnone -fpp -warn all $(OPTS) -DUSEBLAS -DPETSC_VERSION=313
CCOPTS  = -c -DPETSC_VERSION=313
R8OPTS = -r8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -O2 -qopt-report=0 -qopt-report-phase=vec
  CCOPTS := $(CCOPTS) -O
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
F77 = mpif90
LOADER = mpif90
LDOPTS := $(LDOPTS) -cxxlib
F90OPTS = $(F90FLAGS) $(FOPTS) -gen-interfaces
F77OPTS = $(F77FLAGS) $(FOPTS)

# define where you want to locate the mesh adapt libraries
MPIVER=intel2021.1.2-intelmpi2021.3.1
PETSC_VER=petsc-3.13.5
PETSCVER=petsc3.13.5
PETSC_DIR=/p/swim/jchen/PETSC/petsc.20231214
ifeq ($(COM), 1)
	PETSC_ARCH=cplx-spark-m3dc1
	M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
	PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/p/swim/jchen/PETSC/petsc.20231214/cplx-spark-m3dc1/lib -L/p/swim/jchen/PETSC/petsc.20231214/cplx-spark-m3dc1/lib -Wl,-rpath,/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -L/opt/pppl/intel/2023.2.0/ipp/2021.9.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ippcp/2021.8.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/lib -L/opt/pppl/intel/2023.2.0/dal/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ccl/2021.10.0/lib/cpu_gpu_dpcpp -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/lib -L/opt/pppl/intel/2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8 -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/11 -Wl,-rpath,/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -lX11 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lquadmath -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
else
	PETSC_ARCH=real-spark-m3dc1
	M3DC1_SCOREC_LIB=-lm3dc1_scorec
	PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/p/swim/jchen/PETSC/petsc.20231214/real-spark-m3dc1/lib -L/p/swim/jchen/PETSC/petsc.20231214/real-spark-m3dc1/lib -Wl,-rpath,/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -L/opt/pppl/intel/2023.2.0/ipp/2021.9.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ippcp/2021.8.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/lib -L/opt/pppl/intel/2023.2.0/dal/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ccl/2021.10.0/lib/cpu_gpu_dpcpp -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/lib -L/opt/pppl/intel/2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8 -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/11 -Wl,-rpath,/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -lX11 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lquadmath -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
endif

SCOREC_BASE_DIR=/p/swim/jchen/PETSC/core-240527/spark-20231214
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin

SIM_VER=2023.1-240113
MESHGEN_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel9/intel2023.2.0-openmpi4.1.6/2025.0-250217-dev/bin

ifdef SCORECVER
  SCOREC_DIR=$(SCOREC_BASE_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(SCOREC_BASE_DIR)
endif

SCOREC_LIBS= -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
             -Wl,--start-group,-rpath,$(SCOREC_BASE_DIR)/lib -L$(SCOREC_BASE_DIR)/lib \
             -lpumi -lapf -lapf_zoltan -lgmi -llion -lma -lmds -lmth -lparma \
             -lpcu -lph -lsam -lspr -lcrv -Wl,--end-group

LIBS = 	$(SCOREC_LIBS) \
        $(PETSC_WITH_EXTERNAL_LIB)

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
	-I$(SCOREC_DIR)/include \

ifeq ($(ST), 1)
	LIBS += -L$(NETCDF_C_HOME)/lib64 -lnetcdf -L$(NETCDF_F_HOME)/lib -lnetcdff
	INCLUDE += -I$(NETCDF_C_HOME)/include -I$(NETCDF_F_HOME)/include 
endif

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
