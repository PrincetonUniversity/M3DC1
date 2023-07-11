MPIVER=intel2021.1.2-openmpi4.1.0
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
PREFIX=YOUR_DIR
#module load intel/2021.1.2 openmpi/intel-2021.1/4.1.0
#module load gsl/2.6 fftw/intel-2021.1/openmpi-4.1.0/3.3.9 hdf5/intel-2021.1/openmpi-4.1.0/1.10.6
#
#module load rlm/pppl simmodsuite/pppl/17.0-220903
#export SimModSuite_licenseFile="5053@license2.pppl.gov"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/gsl/2.6/x86_64/lib64:/home/PPPL/simmetrix/simmodsuite/17.0-220903/lib/x64_rhel7_gcc48
#cp module load rlm/pppl simmodsuite/pppl/17.0-220903
#cp meshgen/config-files/stellar-FindSimModSuite.cmake meshgen/cmake/FindSimModSuite.cmake
cmake .. \
  -DCMAKE_C_COMPILER=/usr/local/openmpi/4.1.0/intel20211/bin/mpicc \
  -DCMAKE_CXX_COMPILER=/usr/local/openmpi/4.1.0/intel20211/bin/mpicxx \
  -DCMAKE_Fortran_COMPILER=/usr/local/openmpi/4.1.0/intel20211/bin/mpif90 \
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
  -DENABLE_SIMLICENSE=ON \
  -DENABLE_PU=ON \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Debug
