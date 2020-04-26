#export RLM_LICENSE=sdumont12@2800
#export LD_LIBRARY_PATH=/scratch/ntm/software/Simmetrix/extra-libs:$LD_LIBRARY_PATH
#source /scratch/app/modulos/intel-psxe-2019.sh
#module load cmake git intel_psxe/2019
SWTYPE=debug
CMAKETYPE=Debug
PETSCVER=3.9.3
MPIVER=intel-psxe2019
PETSC_DIR=/scratch/ntm/software/petsc/petsc-$PETSCVER
PETSC_ARCH=real-intel-psxe2019
SCOREC_DIR=/scratch/ntm/software/scorec/$MPIVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$SCOREC_DIR/petsc$PETSCVER
SIM_VER=12.0-181027
SIM_DIR=/scratch/ntm/software/Simmetrix/$SIMVER
#/scratch/ntm/software/Simmetrix/12.0-181027
SIM_ARCHOS=x64_rhel6_gcc44
PREFIX=$ZOLTAN_DIR
#PREFIX=$SCOREC_DIR/sim$SIM_VER
cmake .. \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_CXX_COMPILER=mpiicpc \
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
