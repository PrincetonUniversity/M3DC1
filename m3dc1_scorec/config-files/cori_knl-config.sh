# These should be consistent with unstructued/cori_knl.mk,
# but values there will overwrite these when compiled from M3D-C1 "make scorec"
PETSC_DIR=${PETSC_DIR:-/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220107}
SCOREC_COMPLEX=${SCOREC_COMPLEX:-OFF}
if [ $SCOREC_COMPLEX == ON ]
then
  PETSC_ARCH=${PETSC_ARCH:-coriknl-PrgEnvintel6010-craympich7719-master-real}
else
  PETSC_ARCH=${PETSC_ARCH:-coriknl-PrgEnvintel6010-craympich7719-master-cplx}
fi
SCOREC_BASE_DIR=${SCOREC_BASE_DIR:-/global/cfs/cdirs/mp288/jinchen/PETSC/core/upgrade-intel6610-craympich7719-knl}

# This is different from when compiling with M3D-C1 "make scorec"
SCOREC_DIR=${SCOREC_DIR:-$SCOREC_BASE_DIR}

# Always defined here
CMAKETYPE=Release
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=$SCOREC_DIR

#add -DPETSCMASTER for petsc 3.8.3 or higher
cmake .. \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_C_FLAGS=" -g -fPIC -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -fPIC -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fPIC "\
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
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
  -DCMAKE_BUILD_TYPE=$CMAKETYPE \
