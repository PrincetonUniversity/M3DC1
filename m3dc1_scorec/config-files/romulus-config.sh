#module load gcc/7.3.0-bt47fwr  
#module load cmake/3.13.1-hwuwodx
#module load mpich/3.2.1-niuhmad
#module load gsl/2.5-hnrtuln
PETSC_DIR=/lore/seol/petsc-3.8.2
PETSC_ARCH=real-mpich3.2.1-gcc7.3.0
PREFIX=/lore/seol/mpich3.2.1-gcc7.3.0-install
ZOLTAN_DIR=$PREFIX
cmake .. \
  -DCMAKE_C_COMPILER="$MPICH_ROOT/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$MPICH_ROOT/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/$MPICH_ROOT/bin/mpif90" \
  -DCMAKE_C_FLAGS=" -g -O0  -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PREFIX/include \
  -DSCOREC_LIB_DIR=$PREFIX/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DENABLE_PETSC=ON \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DHDF5_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=$PREFIX

