HOST=cori
MPIVER=mpich7.6.2
ARCH=knl
DATE=Jul2019
#SWTYPE=release
CMAKETYPE=Release

PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.9.3
PETSC_ARCH=cori-knl-mpich773-real-nomkl-510

PARMETIS_DIR=${PETSC_DIR}/${PETSC_ARCH}/
ZOLTAN_DIR=${PETSC_DIR}/${PETSC_ARCH}/
PREFIX=/global/project/projectdirs/mp288/cori/scorec/mpich7.7.3/knl-petsc3.9.3/reordered

cmake .. \
  -DCMAKE_C_COMPILER="cc" \
  -DCMAKE_CXX_COMPILER="CC" \
  -DCMAKE_Fortran_COMPILER="ftn" \
  -DCMAKE_C_FLAGS=" -g -O2 -DPETSCMASTER" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -DPETSCMASTER" \
  -DCMAKE_Fortran_FLAGS="-fpic -g -O2" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR="${PREFIX}/include " \
  -DSCOREC_LIB_DIR="${PREFIX}/lib" \
  -DPETSC_DIR="$PETSC_DIR" \
  -DPETSC_ARCH="$PETSC_ARCH" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
