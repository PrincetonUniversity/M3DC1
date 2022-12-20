MPIVER=gnu8.3.3-cuda11.7-mpich8.1.17
CMAKETYPE=Release
PETSC_VER=petsc-3.18.2
PETSCVER=petsc3.18.2
PETSC_DIR=/global/cfs/cdirs/mp288/scorec-pmt/petsc/$PETSC_VER
PETSC_ARCH=real-$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
METIS_DIR=$PETSC_DIR/$PETSC_ARCH
SCOREC_DIR=/global/cfs/cdirs/mp288/scorec-pmt/$MPIVER/$PETSCVER
ZOLTAN_DIR=$SCOREC_DIR
#module load PrgEnv-gnu
#module load cudatoolkit/11.7 craype-accel-nvidia80 cmake/3.22.0
#module unload darshan
PUMI_DIR=$SCOREC_DIR #enter your PUMI install dir
PREFIX=$PUMI_DIR #enter your M3DC1_SCOREC install dir
cmake .. \
  -DCMAKE_C_COMPILER="cc" \
  -DCMAKE_CXX_COMPILER="CC" \
  -DCMAKE_Fortran_COMPILER="ftn" \
  -DCMAKE_C_FLAGS="-g -target-accel=nvidia80 -Ofast -fPIC -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-g -target-accel=nvidia80 -Ofast -fPIC -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-g -target-accel=nvidia80 -Ofast -fPIC" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.so" \
  -DSCOREC_INCLUDE_DIR="$PUMI_DIR/include" \
  -DSCOREC_LIB_DIR="$PUMI_DIR/lib" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
