MPIVER=intel2019u3-openmpi4.0.3
SIM_VER=16.0-220226
SIM_DIR=/usr/pppl/Simmetrix/simmodsuite/$SIM_VER
SIM_ARCHOS=x64_rhel7_gcc48
PREFIX=/p/tsc/m3dc1/lib/SCORECLib/rhel7/$MPIVER/$SIM_VER
PETSC_VER=petsc-3.13.5
PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/PETSC/$PETSC_VER
PETSC_ARCH=real-rhel7-$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
#ssh sunfire10
#cd /u/sseol/develop/core-m3dc1
#module load simmodeler/10.0-220226
#module load simmodsuite/16.0-220226
#module load intel/2019.u3 openmpi/4.0.3 cmake
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DENABLE_SIMMETRIX=ON \
  -DSIM_DISCRETE=OFF \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DENABLE_FIELDSIM=OFF \
  -DBUILD_EXES=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Debug
