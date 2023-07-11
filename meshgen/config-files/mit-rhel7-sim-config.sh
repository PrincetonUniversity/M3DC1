MPIVER=rhel7-gcc6.2.0-openmpi4.0.4
PETSC_DIR=/orcd/nese/psfc/001/software/scorec/petsc-3.18.2
PETSC_ARCH=real-rhel7-gcc6.2.0-openmpi4.0.4
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
PETSCVER=3.19.2
SIM_VER=18.0-230521
SIM_DIR=/orcd/nese/psfc/001/software/simmetrix/SimModSuite$SIM_VER
SIM_ARCHOS=x64_rhel7_gcc48
LAPACK_DIR=$PETSC_DIR/$PETSC_ARCH
PUMI_DIR=/orcd/nese/psfc/001/software/scorec/$MPIVER/sim$SIM_VER
PREFIX=$PUMI_DIR
# module load module load gcc/6.2.0 openmpi/4.0.4
$PETSC_DIR/$PETSC_ARCH/bin/cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS="-O2 -g -Wall -DMIT" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall -DMIT -std=c++11" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR=$PUMI_DIR/include \
  -DSCOREC_LIB_DIR=$PUMI_DIR/lib \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_SIMMETRIX=ON \
  -DENABLE_SIMLICENSE=ON \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DLAPACK_LIB_DIR=$LAPACK_DIR/lib \
  -DTIRPC_LIB_DIR=/usr/lib64 \
  -DENABLE_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Debug
