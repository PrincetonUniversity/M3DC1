MPIVER=intel2019u3-openmpi4.0.3
SIM_VER=16.0-220226
SIM_DIR=/usr/pppl/Simmetrix/simmodsuite/$SIM_VER
SIM_ARCHOS=x64_rhel7_gcc48
PETSC_VER=petsc-3.13.5
PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/PETSC/$PETSC_VER
PETSC_ARCH=real-rhel7-$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
LAPACK_DIR=$LAPACK_HOME
PUMI_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel7/intel2019u3-openmpi4.0.3/16.0-220226/
PREFIX=$PUMI_DIR
#cd ~/develop/core-m3dc1/
#ssh sunfire10
#module load simmodeler/10.0-220226
#module load simmodsuite/16.0-220226
#module load lapack/3.9.0
#module load intel/2019.u3 openmpi/4.0.3 cmake
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
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
