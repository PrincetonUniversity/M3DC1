# edit PREFIX to check installation directory
# change SCOREC_DIR to local install of SCOREC::core
# load a petsc module to get PETSC_DIR and PETSC_ARCH in the environment
# ENABLE_COMPLEX=ON for complex build

PREFIX=$DEVROOT/install/m3dc1_upstream
SCOREC_DIR=$DEVROOT/install/core

CC=`which mpicc`
CXX=`which mpicxx`
FORT=`which mpif90`

cmake .. \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_Fortran_COMPILER=$FORT \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_C_FLAGS=" -g -DDEBUG" \
  -DCMAKE_CXX_FLAGS=" -g -DDEBUG" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$SCOREC_DIR/include \
  -DSCOREC_LIB_DIR=$SCOREC_DIR/lib \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
