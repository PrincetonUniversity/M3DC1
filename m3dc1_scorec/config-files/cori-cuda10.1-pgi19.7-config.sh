HOST=cori
MPIVER=mvapich2.3.1
PETSCVER=petsc3.11.3
PETSC_VER=petsc-3.11.3
#module purge
#module load cuda pgi mvapich2
#salloc -C gpu -t 60 -N 1 -c 10 --gres=gpu:1 -A m1759
ARCH=cuda10.1-pgi19.7
CMAKETYPE=Release
PETSC_DIR=/global/homes/j/jinchen/project/PETSC/$PETSC_VER
PETSC_ARCH=cori-gpu-mvapich2-231-real-nomkl-611
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=/global/project/projectdirs/mp288/$HOST/scorec/$MPIVER/$ARCH-$PETSCVER
ZOLTAN_DIR=$PREFIX
cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpic++" \
  -DCMAKE_Fortran_COMPILER="mpif90" \
  -DCMAKE_C_FLAGS="-DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR="$PREFIX/include" \
  -DSCOREC_LIB_DIR="$PREFIX/lib" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
