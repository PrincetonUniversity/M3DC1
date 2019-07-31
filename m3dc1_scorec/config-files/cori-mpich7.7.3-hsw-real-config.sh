HOST=cori
MPIVER=mpich7.7.6
ARCH=hsw
CMAKETYPE=Release
PETSC_VER=petsc-3.9.3
PETSCVER=petsc3.9.3
#PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/$PETSC_VER
#PETSC_ARCH=cori-hsw-mpich773-real-nomkl-510
PETSC_DIR=/global/homes/j/jinchen/project/PETSC/$PETSC_VER
PETSC_ARCH=cori-hsw-mpich776-cplx-nomkl-510
#/global/project/projectdirs/mp288/cori/petsc/$PETSC_VER
#PETSC_ARCH=real-intel-mpi7.7.3-$ARCH
#load module cray-hdf5-parallel
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/global/project/projectdirs/mp288/cori/scorec/$MPIVER/$ARCH-$PETSCVER
PREFIX=$ZOLTAN_DIR
cmake .. \
  -DCMAKE_C_COMPILER="cc" \
  -DCMAKE_CXX_COMPILER="CC" \
  -DCMAKE_Fortran_COMPILER="ftn" \
  -DCMAKE_C_FLAGS=" -g -O2 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic -g -O2" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DSCOREC_LIB_DIR="$ZOLTAN_DIR/lib" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=ON \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
