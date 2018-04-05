# edit PREFIX to check installation directory
# change SCOREC_DIR to local install of SCOREC::core
# load a petsc module to get PETSC_DIR and PETSC_ARCH in the environment
# ENABLE_COMPLEX=ON for complex build

PREFIX=/lore/seol/mpich3-gcc4.9.2-install
SCOREC_DIR=/lore/seol/mpich3-gcc4.9.2-install
PETSC_DIR=/lore/seol/petsc-3.7.6/real-mpich3
PETSC_ARCH=

CC=`which mpicc`
CXX=`which mpicxx`
FORT=`which mpif90`

cmake .. \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_Fortran_COMPILER=$FORT \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_C_FLAGS=" -g -DDEBUG -I/lore/seol/petsc-3.7.6/include" \
  -DCMAKE_CXX_FLAGS=" -g -DDEBUG -I/lore/seol/petsc-3.7.6/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$SCOREC_DIR/include \
  -DSCOREC_LIB_DIR=$SCOREC_DIR/lib \
  -DPETSC_INCLUDE_DIRS=$PETSC_DIR/$PETSC_ARCH/include \
  -DPETSC_LIBRARIES=$PETSC_DIR/$PETSC_ARCH/lib \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
