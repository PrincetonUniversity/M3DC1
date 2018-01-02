PETSC_DIR=/fasttmp/seol/petsc-3.5.4-real
PETSC_ARCH=openmpi1.6.5
ZOLTAN_DIR=/fasttmp/seol/openmpi-gcc4.4.5-install
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=/fasttmp/seol/openmpi-gcc4.4.5-sim-install

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O2" \
  -DCMAKE_CXX_FLAGS=" -g -O2" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_SIMMETRIX=ON \
  -DSIMMETRIX_INCLUDE_DIR=/net/common/meshSim/latest/include \
  -DSIMMETRIX_LIB_DIR=/net/common/meshSim/latest/lib/x64_rhel5_gcc41 \
  -DSIM_MPI="openmpi18" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
