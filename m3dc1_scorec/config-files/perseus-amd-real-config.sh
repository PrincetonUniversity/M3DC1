MPIVER=intel18.0-mpi2018.3-amd
CMAKETYPE=Release
PETSC_VER=petsc-3.12.4
PETSCVER=petsc3.12.4
PETSC_DIR=/home/jinchen/LIB/petsc
PETSC_ARCH=perseusamd-intelmpi2018-master
#PETSC_ARCH=perseusamd-intelmpi2018-master-cplx
# module load intel/18.0/64/18.0.3.222 intel-mpi/intel/2018.3/64
# module load hdf5/intel-17.0/intel-mpi/1.10.0
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/projects/M3DC1/scorec/$MPIVER/$PETSCVER
PREFIX=$ZOLTAN_DIR
cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_Fortran_COMPILER="mpif90" \
  -DCMAKE_C_FLAGS="-g -O0 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-g -O0 -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR="$PREFIX/include" \
  -DSCOREC_LIB_DIR="$PREFIX/lib" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
