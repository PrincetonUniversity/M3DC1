# These should be consistent with unstructued/stellar.mk,
# but values there will overwrite these when compiled from M3D-C1 "make scorec"
MPIVER=${MPIVER:-intel2021.1.2-intelmpi2021.3.1}
PETSC_VER=${PETSC_VER:-petsc-3.15.5}
PETSCVER=${PETSCVER:-petsc3.15.5}
PETSC_DIR=${PETSC_DIR:-/projects/M3DC1/PETSC/$PETSC_VER}
SCOREC_COMPLEX=${SCOREC_COMPLEX:-OFF}
if [ $SCOREC_COMPLEX == ON ]
then
  PETSC_ARCH=${PETSC_ARCH:-cplx-$MPIVER}
else
  PETSC_ARCH=${PETSC_ARCH:-real-$MPIVER}
fi
SCOREC_BASE_DIR=${SCOREC_BASE_DIR:-/projects/M3DC1/scorec/$MPIVER/$PETSCVER/202209}

# This is different from when compiling with M3D-C1 "make scorec"
SCOREC_DIR=${SCOREC_DIR:-$SCORE_BASE_DIR}

# Always defined here
CMAKETYPE=Release
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=$SCOREC_DIR

cmake3 .. \
  -DCMAKE_C_COMPILER="mpiicc" \
  -DCMAKE_CXX_COMPILER="mpiicpc" \
  -DCMAKE_Fortran_COMPILER="mpiifort" \
  -DCMAKE_C_FLAGS="-g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR="$SCOREC_BASE_DIR/include" \
  -DSCOREC_LIB_DIR="$SCOREC_BASE_DIR/lib" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=$SCOREC_COMPLEX \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
