#module load gcc/7.3.0-bt47fwr  
#module load cmake/3.13.1-hwuwodx
#module load mpich/3.2.1-niuhmad
#module load gsl/2.5-hnrtuln
MPI_DIR=/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.2.1-niuhmadq72otoeruzmhwp6f7b3rqe24y
PREFIX=/lore/seol/mpich3.2.1-niuhmad-install
ZOLTAN_DIR=$PREFIX
PETSC_DIR=/lore/seol/petsc-3.8.2
PETSC_ARCH=real-mpich3.2.1
cmake .. \
  -DCMAKE_C_COMPILER="$MPI_DIR/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$MPI_DIR/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/$MPI_DIR/bin/mpif90" \
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

