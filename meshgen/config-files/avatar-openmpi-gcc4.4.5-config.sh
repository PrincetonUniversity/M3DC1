PETSC_DIR=/lore/seol/petsc-3.5.4
PETSC_ARCH=real-openmpi1.6.5
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/lore/seol/openmpi-gcc4.4.5-install
SIM_ARCHOS=x64_rhel6_gcc44
SIM_VER=14.0-190513dev
PREFIX=/lore/seol/openmpi-gcc4.4.5-$SIM_VER-install
cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/usr/local/openmpi/latest/bin/mpif90" \
  -DCMAKE_C_FLAGS=" -g -O2 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PREFIX/include \
  -DSCOREC_LIB_DIR=$PREFIX/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DENABLE_SIMMETRIX=ON \
  -DSIM_MPI=openmpi165-ib \
  -DSIMMETRIX_INCLUDE_DIR=/net/common/meshSim/$SIM_VER/include \
  -DSIMMETRIX_LIB_DIR=/net/common/meshSim/$SIM_VER/lib/$SIM_ARCHOS \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=/lore/seol/openmpi-gcc4.4.5-$SIM_VER-install \
  -DENABLE_PETSC=ON \
  -DENABLE_TRILINOS=OFF \
  -DLAPACK_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib"



