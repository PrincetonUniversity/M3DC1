CMAKETYPE=Release
PETSCVER=3.18.2
MPIVER=intel-psxe2020-openmpiicc4.1.4
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$PETSCVER
PETSC_ARCH=cplx-$MPIVER
SCOREC_DIR=/scratch/ntm/software/scorec/$MPIVER
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
SIMVER=12.0-181027
SIM_DIR=/scratch/ntm/software/Simmetrix/$SIMVER
SIM_ARCHOS=x64_rhel6_gcc44
PREFIX=$SCOREC_DIR/petsc-$PETSCVER
#source /scratch/app/modulos/intel-psxe-2020.sh
#module load openmpi/icc/4.1.4 intel_psxe/2020 git cmake/3.23.2 
cmake .. \
  -DCMAKE_C_COMPILER="/scratch/app/openmpi/icc/4.1.4/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/scratch/app/openmpi/icc/4.1.4/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/scratch/app/openmpi/icc/4.1.4/bin/mpif90" \
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
