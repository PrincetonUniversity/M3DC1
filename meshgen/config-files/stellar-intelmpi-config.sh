MPIVER=intel2021.1.2-intelmpi2021.3.1
SIM_VER=17.0-220903
SIM_DIR=/home/PPPL/simmetrix/simmodsuite/$SIM_VER
SIM_ARCHOS=x64_rhel7_gcc48
PETSC_VER=petsc-3.13.5
PETSCVER=petsc3.13.5
PETSC_DIR=/projects/M3DC1/PETSC/$PETSC_VER
PETSC_ARCH=real-$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
LAPACK_DIR=$PETSC_DIR/$PETSC_ARCH
PUMI_DIR=/projects/M3DC1/scorec/$MPIVER/$SIM_VER
PREFIX=/projects/M3DC1/scorec/$MPIVER/$SIM_VER
#module purge
#module load intel/2021.1.2 intel-mpi/intel/2021.3.1 cmake/3.19.7
#module load fftw/intel-2021.1/intel-mpi/3.3.9
#module load hdf5/intel-2021.1/intel-mpi/1.10.6
#module load rlm/pppl simmodsuite/pppl/17.0-220903
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64:/home/PPPL/simmetrix/simmodsuite/17.0-220903/lib/x64_rhel7_gcc48
#unset I_MPI_HYDRA_BOOTSTRAP
#unset I_MPI_PMI_LIBRARY
#OPTFLAGS=""
cmake3 .. \
  -DCMAKE_C_COMPILER="mpiicc" \
  -DCMAKE_CXX_COMPILER="mpiicpc" \
  -DCMAKE_Fortran_COMPILER="mpiifort" \
  -DCMAKE_C_FLAGS="-O2 -g -DOLDMA -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-O2 -g -DOLDMA -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR=$PUMI_DIR/include \
  -DSCOREC_LIB_DIR=$PUMI_DIR/lib \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_SIMMETRIX=ON \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DLAPACK_LIB_DIR=$LAPACK_DIR/lib \
  -DENABLE_TESTING=OFF \
  -DENABLE_PPPL=ON \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Debug
