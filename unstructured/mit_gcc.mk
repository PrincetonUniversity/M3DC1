FOPTS = -c -fdefault-real-8 -fdefault-double-8 -cpp -DPETSC_VERSION=319 -DUSEBLAS $(OPTS) 
CCOPTS  = -c -O -DPETSC_VERSION=319
R8OPTS = -fdefault-real-8 -fdefault-double-8

ifeq ($(OPT), 1)
  FOPTS  := $(FOPTS) -w -fallow-argument-mismatch -O2
  CCOPTS := $(CCOPTS) -O
else
  FOPTS := $(FOPTS) -g 
endif

ifeq ($(PAR), 1)
  FOPTS := $(FOPTS) -DUSEPARTICLES
endif

ifeq ($(OMP), 1)
  LDOPTS := $(LDOPTS) -mp
  FOPTS  := $(FOPTS)  -mp
  CCOPTS := $(CCOPTS) -mp
endif

CC = mpicc
CPP = mpicxx
F90 = mpif90
F77 = mpif90
LOADER =  mpif90
FOPTS := $(FOPTS)

F90OPTS = $(F90FLAGS) $(FOPTS) 
F77OPTS = $(F77FLAGS) $(FOPTS)

PETSC_DIR=/orcd/nese/psfc/001/jinchen/petsc/petsc20230612
ifeq ($(COM), 1)
  PETSC_ARCH=mit-gcc-openmpi-cplx
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi-cplx/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi-cplx/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -lpetsc -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
else
  ifeq ($(ST), 1)
  PETSC_ARCH=mit-gcc-openmpi
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
  else
  PETSC_ARCH=mit-gcc-openmpi
  PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/orcd/nese/psfc/001/jinchen/petsc/petsc20230612/mit-gcc-openmpi/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib64 -L/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/gcc-12.2.0-or6pfydukwucqlbwbijl5pgpgknm4jc5/lib/gcc/x86_64-pc-linux-gnu/12.2.0 -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/openmpi-4.1.4-3r4zaihkaqj2gmfvtzk4adiu3qxlzgj5/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/hwloc-2.8.0-a56oj35bkhqi7rpsxyrzv2cvjhk6f4nl/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/ucx-1.13.1-xxapazpv4rciyeuajxom5vfll4djhakq/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libevent-2.1.12-dnasb7atyzwlagnyyrplzk5if6efrfbe/lib -Wl,-rpath,/nfs/software001/home/software-r8-x86_64/spack-20230328/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/pmix-4.1.2-d7f7fwzoomvmc6hwotduhlzmdoc6oz7o/lib -lpetsc -ldmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu -lsuperlu_dist -lfftw3_mpi -lfftw3 -lflapack -lfblas -lzoltan -lparmetis -lmetis -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgsl -lgslcblas -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl
  endif
endif

SCOREC_BASE_DIR=/orcd/nese/psfc/001/scorec/gcc12.2.0-openmpi4.1.4/petsc3.19.2
#SCOREC_BASE_DIR=/orcd/nese/psfc/001/software/scorec/gcc12.2.0-openmpi4.1.4/petsc-3.19.2
SCOREC_UTIL_DIR=$(SCOREC_BASE_DIR)/bin
PUMI_DIR=$(SCOREC_BASE_DIR)
PUMI_LIB = -lpumi -lapf -lapf_zoltan -lcrv -lsam -lspr -lmth -lgmi -lma -lmds -lparma -lpcu -lph -llion
M3DC1_SCOREC_DIR=$(SCOREC_BASE_DIR)

ifdef SCORECVER
  SCOREC_DIR=$(M3DC1_SCOREC_DIR)/$(SCORECVER)
else
  SCOREC_DIR=$(M3DC1_SCOREC_DIR)
endif

ifeq ($(COM), 1)
  M3DC1_SCOREC_LIB=-lm3dc1_scorec_complex
else
  M3DC1_SCOREC_LIB=-lm3dc1_scorec
endif

SCOREC_LIB = -L$(SCOREC_DIR)/lib $(M3DC1_SCOREC_LIB) \
            -Wl,--start-group,-rpath,$(PUMI_DIR)/lib -L$(PUMI_DIR)/lib \
           $(PUMI_LIB) -Wl,--end-group

LIBS = 	\
	$(SCOREC_LIB) \
        $(PETSC_WITH_EXTERNAL_LIB) \

INCLUDE = -I$(PETSC_DIR)/include \
        -I$(PETSC_DIR)/$(PETSC_ARCH)/include \

ifeq ($(ST), 1)
  NETCDF_CDIR=/orcd/nese/psfc/001/software/spack/2023-05-01-physics_rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/netcdf-c-4.9.2-kogljzthr4j5xi2bffxdxp6smzie7b62
  NETCDF_FDIR=/orcd/nese/psfc/001/software/spack/2023-05-01-physics_rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/netcdf-fortran-4.6.0-cylyrw3rlsvkzgfwy7pcx344xbftebgw/
  LIBS += -Wl,--start-group -L$(NETCDF_CDIR)/lib -Wl,-rpath,$(NETCDF_CDIR)/lib -lnetcdf -Wl,--end-group \
          -Wl,--start-group -L$(NETCDF_FDIR)/lib -Wl,-rpath,$(NETCDF_FDIR)/lib -lnetcdff -lz -Wl,--end-group
  INCLUDE += -I$(NETCDF_CDIR)/include \
             -I$(NETCDF_FDIR)/include
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
	$(F90) $(F90OPTS) $(INCLUDE) $< -o $@
