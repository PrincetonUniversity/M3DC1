MPIVER=gcc12.2.0-openmpi4.1.4
#PETSC_DIR=/orcd/nese/psfc/001/software/scorec/petsc-3.19.2
#PETSC_ARCH=real-$MPIVER
PETSC_DIR=/orcd/nese/psfc/001/jinchen/petsc/petsc20230612
PETSC_ARCH=mit-gcc-openmpi
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=$PETSC_DIR/$PETSC_ARCH
PETSCVER=3.19.2
SIM_VER=18.0-230521
SIM_DIR=/orcd/nese/psfc/001/software/simmetrix/SimModSuite$SIM_VER
SIM_ARCHOS=x64_rhel8_gcc83
PUMI_DIR=/orcd/nese/psfc/001/software/scorec/gcc12.2.0-openmpi4.1.4/sim18.0-230521
#/orcd/nese/psfc/001/software/scorec/$MPIVER/sim$SIM_VER
LAPACK_DIR=$PETSC_DIR/$PETSC_ARCH
TIRPC_DIR=/orcd/nese/psfc/001/software/spack/2023-05-01-physics_rpp/spack/opt/spack/linux-rocky8-x86_64/gcc-12.2.0/libtirpc-1.2.6-rbpmdt36beybzvwcrrtxk4n2unp7l6cm
PREFIX=$PUMI_DIR
#module load cmake gcc/12.2.0-x86_64
#module use /orcd/nese/psfc/001/software/modulefiles/spack-modules
#module load openmpi-4.1.4-gcc-12.2.0-3r4zaih
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_C_FLAGS="-O2 -g -Wall -DMIT" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall -DMIT" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DSCOREC_INCLUDE_DIR=$PUMI_DIR/include \
  -DSCOREC_LIB_DIR=$PUMI_DIR/lib \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_SIMMETRIX=ON \
  -DENABLE_SIMLICENSE=ON \
  -DSIMMETRIX_INCLUDE_DIR=$SIM_DIR/include \
  -DSIMMETRIX_LIB_DIR=$SIM_DIR/lib/$SIM_ARCHOS \
  -DLAPACK_LIB_DIR=$LAPACK_DIR/lib \
  -DTIRPC_LIB_DIR=$TIRPC_DIR/lib \
  -DENABLE_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Debug
