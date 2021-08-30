PETSC_DIR=/lore/seol/petsc-3.13.5
PETSC_ARCH=real-gcc4.8.5-v5m6xwi-mpich3.2.1-geowaxe
PARMETIS_INSTALL_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_INSTALL_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=/lore/seol/rhel7-gcc4.8.5-v5m6xwi-mpich3.2.1-geowaxe
#module load gcc/4.8.5-v5m6xwi mpich/3.2.1-geowaxe cmake gsl/2.5-px4dg7h
#module unload zlib/1.2.11-vhzh5cf
#add -DPETSCMASTER for petsc 3.8.3 or higher
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PREFIX/include \
  -DSCOREC_LIB_DIR=$PREFIX/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_INSTALL_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_INSTALL_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_INSTALL_DIR/lib/libmetis.a" \
  -DENABLE_PETSC=ON \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=$PREFIX
