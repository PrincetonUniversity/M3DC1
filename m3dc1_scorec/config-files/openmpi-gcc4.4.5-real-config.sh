PETSC_DIR=/lore/seol/petsc-3.7.6
PETSC_ARCH=real-openmpi1.6.5
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
#module load openmpi/1.6.5-ib cmake git-2.21.0-gcc-4.9.2-5sxljib
PUMI_DIR=/lore/seol/openmpi1.6.5-petsc3.7.6-install
ZOLTAN_DIR=$PUMI_DIR
PREFIX=$PUMI_DIR
cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/usr/local/openmpi/latest/bin/mpif90" \
  -DCMAKE_C_FLAGS="-g -O0 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-g -O0 -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PUMI_DIR/include \
  -DSCOREC_LIB_DIR=$PUMI_DIR/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DENABLE_PETSC=ON \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=$PREFIX
