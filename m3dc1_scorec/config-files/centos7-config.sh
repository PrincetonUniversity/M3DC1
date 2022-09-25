# These should be consistent with unstructued/centos7.mk,
# but values there will overwrite these when compiled from M3D-C1 "make scorec"
MPIVER=${MPIVER:-intel2019u3-openmpi4.0.3}
PETSCVER=${PETSC_VER:-petsc3.13.5}
PETSC_VER=${PETSCVER:-petsc-3.13.5}
PETSC_DIR=${PETSC_DIR:-/p/tsc/m3dc1/lib/SCORECLib/PETSC/$PETSC_VER}
SCOREC_COMPLEX=${SCOREC_COMPLEX:-OFF}
if [ $SCOREC_COMPLEX == ON ]
then
  PETSC_ARCH=${PETSC_ARCH:-cplx-rhel7-$MPIVER}
else
  PETSC_ARCH=${PETSC_ARCH:-real-rhel7-$MPIVER}
fi

# This is different from when compiling with M3D-C1 "make scorec"
SCOREC_DIR=${SCOREC_DIR:-/p/tsc/m3dc1/lib/SCORECLib/rhel7/$MPIVER/$PETSCVER}

# Always defined here
CMAKETYPE=Release
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=$SCOREC_DIR

cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR="$SCOREC_BASE_DIR/include" \
  -DSCOREC_LIB_DIR="$SCOREC_BASE_DIR/lib" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_PETSC=ON \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_COMPLEX=$SCOREC_COMPLEX \
  -DENABLE_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
