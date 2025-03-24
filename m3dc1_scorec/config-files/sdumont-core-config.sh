CMAKETYPE=Release
PETSCVER=3.18.2
MPIVER=intel-psxe2020-openmpiicc4.1.4
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$PETSCVER
PETSC_ARCH=real-$MPIVER
SCOREC_DIR=/scratch/ntm/software/scorec/$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
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
  -DCMAKE_C_FLAGS="-ftz -fPIC -O" \
  -DCMAKE_CXX_FLAGS="-shared-intel -ftz -fPIC -O" \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DBUILD_EXES=ON \
  -DIS_TESTING=OFF \
  -DENABLE_ZOLTAN=ON \
  -DENABLE_SIMMETRIX=OFF \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
