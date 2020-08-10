#export RLM_LICENSE=sdumont12@2800
#export LD_LIBRARY_PATH=/scratch/ntm/software/Simmetrix/extra-libs:$LD_LIBRARY_PATH
#source /scratch/app/modulos/intel-psxe-2019.sh
#module load openmpi/icc/4.0.4 cmake git intel_psxe/2019
SWTYPE=debug
CMAKETYPE=Debug
PETSCVER=3.9.4
MPIVER=intel-psxe2019-openmpiicc4.0.4
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$PETSCVER
PETSC_ARCH=cplx-$MPIVER
PREFIX=/scratch/ntm/software/scorec/$MPIVER/petsc$PETSCVER
PARMETIS_DIR=$PETSC_ARCH/$PETSC_VER
ZOLTAN_DIR=$PETSC_ARCH/$PETSC_VER

HDF5_DIR=$PETSC_DIR/$PETSC_ARCH
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -ftz -fPIC -O -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -shared-intel -ftz -fPIC -O -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-assume no2underscores -ftz -fPIC -O" \
  -DSCOREC_INCLUDE_DIR="$PREFIX/include" \
  -DSCOREC_LIB_DIR="$PREFIX/lib" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_PETSC=ON \
  -DENABLE_TRILINOS=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DENABLE_COMPLEX=ON \
  -DENABLE_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
