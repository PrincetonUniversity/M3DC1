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
PREFIX=/scratch/ntm/software/scorec/$MPIVER/petsc$PETSCVER
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PREFIX
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
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBUILD_EXES=ON \
  -DIS_TESTING=OFF \
  -DENABLE_ZOLTAN=ON \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE
