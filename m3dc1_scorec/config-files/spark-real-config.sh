MPIVER=
CMAKETYPE=Release
PETSCVER=petsc3.19.6
PETSC_VER=petsc-3.19.6
PETSC_DIR=/p/swim/jchen/PETSC/petsc.20231214
PETSC_ARCH=real-spark-m3dc1
PARMETIS_DIR=/p/swim/jchen/PETSC/petsc.20231214/real-spark-m3dc1
METIS_DIR=/p/swim/jchen/PETSC/petsc.20231214/real-spark-m3dc1
ZOLTAN_INSTALL_DIR=/p/swim/jchen/PETSC/petsc.20231214/real-spark-m3dc1
NETCDF_F_HOME=/opt/pppl/spack-pkgs/linux-rocky9-zen3/intel-2023.2.0/netcdf-fortran-4.6.0-hz3b2mq6ld2u2u6ktdodbx72k37opcdf
NETCDF_C_HOME=/opt/pppl/spack-pkgs/linux-rocky9-zen3/intel-2023.2.0/netcdf-c-4.9.2-3wxlqzxkv3tfm4rpo5utrn2lc3yejtnu
SCOREC_DIR=/p/swim/jchen/PETSC/core-240527/spark-20231214
PREFIX=$SCOREC_DIR

cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS=" -g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -std=c++11 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_INCLUDE_DIR="$METIS_DIR/include" \
  -DMETIS_LIBRARY="$METIS_DIR/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_PETSC=ON \
  -DSCOREC_CXX_WARNINGS=OFF \
  -DSCOREC_CXX_OPTIMIZE=OFF \
  -DBUILD_EXES=ON \
  -DIS_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE \
  -DCMAKE_INSTALL_PREFIX=$PREFIX
