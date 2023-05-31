MPIVER=gnu8.3.3-cuda11.7-mpich8.1.17
CMAKETYPE=Release
PETSC_VER=petsc-3.18.2
PETSCVER=petsc3.18.2
PETSC_DIR=/global/cfs/cdirs/mp288/scorec-pmt/petsc/$PETSC_VER
PETSC_ARCH=real-$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
METIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/global/cfs/cdirs/mp288/scorec-pmt/$MPIVER/$PETSCVER
PREFIX=$ZOLTAN_DIR #enter your M3DC1_SCOREC install dir
#module load PrgEnv-gnu
#module load cudatoolkit/11.7 craype-accel-nvidia80 cmake/3.22.0
#module unload darshan
cmake .. \
  -DCMAKE_C_COMPILER="cc" \
  -DCMAKE_CXX_COMPILER="CC" \
  -DSCOREC_CXX_WARNINGS=OFF \
  -DSCOREC_CXX_OPTIMIZE=OFF \
  -DENABLE_ZOLTAN=ON \
  -DBUILD_EXES=ON \
  -DIS_TESTING=OFF \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_INCLUDE_DIR="$METIS_DIR/include" \
  -DMETIS_LIBRARY="$METIS_DIR/lib/libmetis.so" \
  -DBUILD_EXES=ON \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
