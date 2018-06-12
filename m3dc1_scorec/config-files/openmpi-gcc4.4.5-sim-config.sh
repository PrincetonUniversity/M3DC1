PETSC_DIR=/lore/seol/petsc-3.5.4
PETSC_ARCH=real-openmpi1.6.5
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
ZOLTAN_DIR=/lore/seol/openmpi-gcc4.4.5-install
PREFIX=/lore/seol/openmpi-gcc4.4.5-sim-install
SIM_ARCHOS=x64_rhel6_gcc44
cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/usr/local/openmpi/latest/bin/mpif90" \
  -DCMAKE_C_FLAGS=" -g -O2 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DSCOREC_INCLUDE_DIR=$PREFIX/include \
  -DSCOREC_LIB_DIR=$PREFIX/lib \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DENABLE_SIMMETRIX=ON \
  -DSIM_MPI=openmpi165-ib \
  -DSIMMETRIX_INCLUDE_DIR=/net/common/meshSim/11.0-170826/include \
  -DSIMMETRIX_LIB_DIR=/net/common/meshSim/11.0-170826/lib/$SIM_ARCHOS \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=/lore/seol/openmpi-gcc4.4.5-sim-install \
  -DENABLE_PETSC=ON \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DHDF5_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DENABLE_TRILINOS=OFF \
  -DTRILINOS_INCLUDE_DIR="/fasttmp/seol/openmpi-gcc4.4.5-install/include" \
  -DTRILINOS_LIB_DIR="/fasttmp/seol/openmpi-gcc4.4.5-install/lib" \
  -DLAPACK_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DBOOST_LIB_DIR="/fasttmp/seol/openmpi-gcc4.4.5-install/lib" \
  -DSTDCPP_LIBRARY="/usr/lib/gcc/x86_64-linux-gnu/4.4.5/libstdc++.a" \
  -DNETCDF_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libnetcdf.a" 


