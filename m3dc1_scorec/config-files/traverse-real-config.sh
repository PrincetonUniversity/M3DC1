MPIVER=pgi19.9-openmpi4.0.2
CMAKETYPE=Release
PETSC_VER=petsc-3.12.4
PETSCVER=petsc3.12.4
PETSC_DIR=/home/jinchen/project/PETSC/petsc
#   PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-cplx-master
PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-master
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
SCOREC_DIR=/projects/M3DC1/scorec/$MPIVER/$PETSCVER
ZOLTAN_DIR=$SCOREC_DIR
PREFIX=$SCOREC_DIR
#load module load pgi openmpi/pgi-19.5/4.0.2rc1/64
#add -DPETSCMASTER for petsc 3.8.3 or higher
cmake3 .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR="$SCOREC_DIR/include" \
  -DSCOREC_LIB_DIR="$SCOREC_DIR/lib" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_PETSC=ON \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
