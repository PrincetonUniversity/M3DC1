#EDIT PETSC_DIR, PETSC_ARCH, PREFIX
#FOR COMPLEX VERSION, SET 
PETSC_DIR=/fasttmp/seol/petsc-3.5.4
PETSC_ARCH=real-mpich3-gcc4.7.0
PREFIX=/fasttmp/seol/mpich3-gcc4.7.0-install

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/mpich3/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/mpich3/latest/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/usr/local/mpich3/latest/bin/mpif90" \
  -DCMAKE_C_FLAGS=" -g -O2 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR="$PREFIX/include" \
  -DSCOREC_LIB_DIR="$PREFIX/lib" \
  -DZOLTAN_LIBRARY="$PREFIX/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DHDF5_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_MESHGEN=OFF \
  -DSIMMETRIX_INCLUDE_DIR=/net/common/meshSim/latest/include \
  -DSIMMETRIX_LIB_DIR=/net/common/meshSim/latest/lib/x64_rhel5_gcc41 \
  -DENABLE_TRILINOS=OFF \
  -DTRILINOS_INCLUDE_DIR="/fasttmp/seol/openmpi-gcc4.7.0-install/include" \
  -DTRILINOS_LIB_DIR="/fasttmp/seol/openmpi-gcc4.7.0-install/lib" \
  -DLAPACK_LIB_DIR="/fasttmp/seol/petsc-3.5.4-real/openmpi1.6.5/lib" \
  -DBOOST_LIB_DIR="/fasttmp/seol/openmpi-gcc4.7.0-install/lib" \
  -DSTDCPP_LIBRARY="/usr/local/gcc/4.7.0/lib64/libstdc++.a" \
  -DNETCDF_LIBRARY="/fasttmp/seol/petsc-3.5.4-real/openmpi1.6.5/lib/libnetcdf.a" \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"

