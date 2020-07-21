MPIVER=intel2019u3-openmpi4.0.1
ZOLTAN_DIR=/p/tsc/m3dc1/lib/SCORECLib/rhel6/$MPIVER/petsc3.9.4
PETSC_DIR=/p/tsc/m3dc1/lib/SCORECLib/PETSC/petsc-3.9.4
PETSC_ARCH=real-intel2019u3-openmpi4.0.1
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
SIM_VER=14.0-190402dev
SIM_DIR=/usr/pppl/Simmetrix/simmodsuite/$SIM_VER
SIM_ARCHOS=x64_rhel6_gcc44
PREFIX=/p/tsc/m3dc1/lib/SCORECLib/rhel6/$MPIVER
#cd /u/sseol/develop/core-sim-m3dc1
#module load intel/2019.u3 openmpi/4.0.1 
#module load lapack simmodeler/7.0-190402dev simmodsuite/14.0-190402dev
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DENABLE_SIMMETRIX=ON \
  -DENABLE_FIELDSIM=OFF \
  -DBUILD_EXES=OFF \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
