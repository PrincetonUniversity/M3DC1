#!/bin/bash

CC=$(which mpicc)
CXX=$(which mpicxx)
FORT=$(which mpif90)

cmake \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_Fortran_COMPILER=$FORT \
    -DCMAKE_INSTALL_PREFIX=$DEVROOT/install/m3dc1 \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
    -DUSE3D=1 \
    -DUSESCOREC=1 \
    -DM3DC1_ARCH=scorec \
    -DSCOREC_LIB_DIR=$DEVROOT/install/core/lib \
    -DSCOREC_INCLUDE_DIR=$DEVROOT/install/core/include \
    -DPETSC_DIR=$PETSC_DIR \
    -DPETSC_ARCH=$PETSC_ARCH \
    -DENABLE_TESTING=1 \
    ..
