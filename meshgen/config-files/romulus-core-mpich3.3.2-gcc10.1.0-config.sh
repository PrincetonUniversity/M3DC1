MPIVER=mpich3.3.2
PETSC_DIR=/lore/seol/petsc-3.13.5 
PETSC_ARCH=real-$MPIVER-gcc10.1.0
PARMETIS_DIR=$PARMETIS_ROOT
METIS_DIR=$METIS_ROOT
ZOLTAN_DIR=$ZOLTAN_ROOT
SIM_ARCHOS=x64_rhel7_gcc48
SIM_VER=17.0-210808dev
SIM_DIR=$SIMMETRIX_SIMMODSUITE_ROOT
PREFIX=/lore/seol/$MPIVER-gcc10.1.0-$SIM_VER-install
#module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
#module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
#module load gcc/10.1.0 mpich/3.3.2 zoltan/3.83-int32 cmake/3.20.0 simmetrix-simmodsuite/17.0-210808dev simmetrix/simModeler/11.0-210812-dev
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.so" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_INCLUDE_DIR="$METIS_DIR/include" \
  -DMETIS_LIBRARY="$METIS_DIR/lib/libmetis.so" \
  -DENABLE_SIMMETRIX=ON \
  -DSIM_MPI=$MPIVER \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_EXES=OFF \
  -DCMAKE_INSTALL_PREFIX=$PREFIX
