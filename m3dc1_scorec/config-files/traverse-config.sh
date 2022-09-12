HOST=traverse
#/projects/M3DC1/scorec/traverse/pgi20.4-openmpi4.0.4/petsc3.16.2
MPIVER=pgi20.4-openmpi4.0.4
CMAKETYPE=Release
PETSC_VER=petsc-3.13.4
PETSCVER=petsc3.13.4
PETSC_DIR=/home/jinchen/project/PETSC/petsc
PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-master
#PETSC_ARCH=traverse-pgi-openmpi-199-gpu-cuda-cplx-master
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/projects/M3DC1/scorec/$MPIVER/$PETSCVER
SCOREC_DIR=/projects/M3DC1/scorec/$MPIVER/$PETSCVER
PREFIX=$SCOREC_DIR/202209
# module load pgi/20.4/64 openmpi/pgi-20.4/4.0.4/64
#add -DPETSCMASTER for petsc 3.8.3 or higher
cmake3 .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O0 -DOLDMA -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -DOLDMA -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR="$SCOREC_DIR/include" \
  -DSCOREC_LIB_DIR="$SCOREC_DIR/lib" \
  -DENABLE_ZOLTAN=ON \
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
