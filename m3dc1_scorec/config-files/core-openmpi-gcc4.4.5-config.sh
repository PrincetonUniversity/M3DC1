PETSC_DIR=/lore/seol/petsc-3.7.6
PETSC_ARCH=real-openmpi1.6.5
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
#module load openmpi/1.6.5-ib cmake git-2.21.0-gcc-4.9.2-5sxljib
PUMI_DIR=/lore/seol/openmpi1.6.5-petsc3.7.6-install
ZOLTAN_DIR=$PUMI_DIR
PREFIX=$PUMI_DIR
cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DSCOREC_CXX_WARNINGS=OFF \
  -DSCOREC_CXX_OPTIMIZE=OFF \
  -DBUILD_EXES=ON \
  -DIS_TESTING=OFF \
  -DMESHES=/lore/seol/meshes \
  -DMPIRUN=/usr/local/openmpi/latest/bin/mpirun \
  -DCMAKE_BUILD_TYPE=Debug

