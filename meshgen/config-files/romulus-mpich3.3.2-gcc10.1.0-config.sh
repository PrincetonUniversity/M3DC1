MPIVER=mpich3.3.2
PETSC_DIR=/lore/seol/petsc-3.13.5 
PETSC_ARCH=real-$MPIVER-gcc10.1.0
PARMETIS_DIR=$PARMETIS_ROOT
METIS_DIR=$METIS_ROOT
ZOLTAN_DIR=$ZOLTAN_ROOT
SIM_ARCHOS=x64_rhel7_gcc48
SIM_VER=17.0-210808dev
#17.0-210808dev
SIM_DIR=$SIMMETRIX_SIMMODSUITE_ROOT
PREFIX=/lore/seol/$MPIVER-gcc10.1.0-$SIM_VER-install
#module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
#module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
#module load gcc/10.1.0 mpich/3.3.2 zoltan/3.83-int32 cmake/3.20.0
#module load simmetrix-simmodsuite/16.0-210606dev
# simmetrix-simmodsuite/17.0-210808dev simmetrix/simModeler/11.0-210812-dev
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O2 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PREFIX/include \
  -DSCOREC_LIB_DIR=$PREFIX/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.so" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_LIBRARY="$METIS_DIR/lib/libmetis.so" \
  -DENABLE_SIMMETRIX=ON \
  -DSIM_MPI=$MPIVER \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DENABLE_PETSC=ON \
  -DLAPACK_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib"
