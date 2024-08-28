# Modify these paths for your system.
setenv MPIHOME /usr/pppl/intel/2015-pkgs/openmpi-1.8.4
setenv PETSC_DIR /p/swim/jchen/PETSC/petsc-3.5.3
setenv PETSC_ARCH portalr6-intel-openmpi-1.8.4
setenv BOOST_DIR /usr/pppl/boost/1.52.0
setenv HDF5_HOME /usr/pppl/intel/2015-pkgs/openmpi-1.8-pkgs/hdf5-1.8.14-parallel
setenv NETCDF_DIR /p/tsc/m3dc1/lib/SCORECLib/rhel6/trilinos
setenv MKL_DIR /usr/pppl/intel/2015.u1/composerxe/mkl
setenv INTEL_LICENSE_FILE /usr/pppl/intel/licenses/server.lic

/p/tsc/m3dc1/lib/SCORECLib/rhel6/cmake-3.3.2/bin/cmake \
\
 -D Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON \
 -D CMAKE_INSTALL_PREFIX:PATH="/p/tsc/m3dc1/lib/SCORECLib/rhel6/trilinos/real" \
 -D CMAKE_BUILD_TYPE:STRING=DEBUG \
 -D TPL_ENABLE_MPI:BOOL=ON \
 -D MPI_BASE_DIR:PATH=$MPIHOME \
 -D CMAKE_C_COMPILER:STRING="/usr/pppl/intel/2015.u1/bin/icc" \
 -D CMAKE_CXX_COMPILER:STRING="/usr/pppl/intel/2015.u1/bin/icpc" \
 -D CMAKE_Fortran_COMPILER:STRING="/usr/pppl/intel/2015.u1/bin/ifort" \
# -D CMAKE_C_COMPILER:STRING="/usr/pppl/intel/2015-pkgs/openmpi-1.8.4/bin/mpicc" \
# -D CMAKE_CXX_COMPILER:STRING="/usr/pppl/intel/2015-pkgs/openmpi-1.8.4/bin/mpicxx" \
# -D CMAKE_Fortran_COMPILER:STRING="/usr/pppl/intel/2015-pkgs/openmpi-1.8.4/bin/mpif90" \
 -D CMAKE_C_FLAGS:STRING="-O2 -g " \
 -D CMAKE_CXX_FLAGS:STRING="-O2 -std=c++11 -ggdb -Wno-sign-compare" \
 -D Trilinos_ENABLE_CXX11:BOOL=ON \
 -D Trilinos_CXX11_FLAGS:STRING="-std=c++11" \
 -D CMAKE_ENABLE_Fortran:BOOL=ON \
 -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
 -D BUILD_SHARED_LIBS:BOOL=ON \
# -D TPL_FIND_SHARED_LIBS:BOOL=OFF \
# -D Trilinos_LINK_SEARCH_START_STATIC:BOOL=ON \
 -D Trilinos_EXTRA_LINK_FLAGS:STRING="-ldl" \
\
 -D Trilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=OFF \
\
 -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -D Trilinos_ENABLE_TESTS:BOOL=OFF \
 -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
 -D Trilinos_ENABLE_Teuchos:BOOL=ON \
 -D Trilinos_ENABLE_Shards:BOOL=ON \
 -D Trilinos_ENABLE_Sacado:BOOL=ON \
 -D Trilinos_ENABLE_Epetra:BOOL=ON \
 -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
 -D Trilinos_ENABLE_Ifpack:BOOL=ON \
 -D Trilinos_ENABLE_AztecOO:BOOL=ON \
 -D Trilinos_ENABLE_Amesos:BOOL=ON \
 -D Trilinos_ENABLE_Anasazi:BOOL=ON \
 -D Trilinos_ENABLE_Belos:BOOL=ON \
 -D Trilinos_ENABLE_ML:BOOL=ON \
 -D Trilinos_ENABLE_Phalanx:BOOL=ON \
 -D Trilinos_ENABLE_Intrepid:BOOL=ON \
 -D Trilinos_ENABLE_NOX:BOOL=ON \
 -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
 -D Trilinos_ENABLE_Thyra:BOOL=ON \
 -D Trilinos_ENABLE_Rythmos:BOOL=ON \
 -D Trilinos_ENABLE_MOOCHO:BOOL=ON \
 -D Trilinos_ENABLE_Stokhos:BOOL=ON \
 -D Trilinos_ENABLE_Piro:BOOL=ON \
 -D Trilinos_ENABLE_Teko:BOOL=ON \
\
 -D Trilinos_ENABLE_STKIO:BOOL=ON \
 -D Trilinos_ENABLE_STKMesh:BOOL=ON \
 -D TPL_ENABLE_Boost:BOOL=ON \
 -D Boost_INCLUDE_DIRS:FILEPATH="$BOOST_DIR/include" \
 -D Boost_LIBRARY_DIRS:FILEPATH="$BOOST_DIR/lib" \
 -D TPL_ENABLE_BoostLib:BOOL=ON \
 -D BoostLib_INCLUDE_DIRS:FILEPATH="$BOOST_DIR/include" \
 -D BoostLib_LIBRARY_DIRS:FILEPATH="$BOOST_DIR/lib" \
\
 -D Trilinos_ENABLE_SEACAS:BOOL=ON \
 -D TPL_ENABLE_X11:BOOL=OFF \
 -D TPL_ENABLE_Matio:BOOL=OFF \
 -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
 -D Trilinos_ENABLE_SEACASExodus:BOOL=ON \
 -D TPL_ENABLE_Netcdf:BOOL=ON \
 -D Netcdf_INCLUDE_DIRS:PATH="$NETCDF_DIR/include" \
 -D Netcdf_LIBRARY_DIRS:PATH="$NETCDF_DIR/lib" \
 -D TPL_ENABLE_HDF5:BOOL=ON \
 -D HDF5_INCLUDE_DIRS:PATH="$HDF5_HOME/include" \
 -D HDF5_LIBRARY_DIRS:PATH="$HDF5_HOME/lib" \
\
 -D Trilinos_ENABLE_Tpetra:BOOL=ON \
 -D Trilinos_ENABLE_Kokkos:BOOL=ON \
 -D HAVE_INTREPID_KOKKOSCORE:BOOL=ON \
 -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
 -D Trilinos_ENABLE_Amesos2:BOOL=ON \
 -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
 -D Trilinos_ENABLE_MueLu:BOOL=ON \
 -D Amesos2_ENABLE_KLU2:BOOL=ON \
 -D Amesos2_ENABLE_MUMPS:BOOL=ON \
 -D Amesos2_ENABLE_SuperLU:BOOL=ON \
 -D Amesos2_ENABLE_SuperLUDist:BOOL=ON \
\
 -D TPL_ENABLE_MKL:BOOL=ON \
 -D MKL_INCLUDE_DIRS:PATH="$MKL_DIR/include" \
 -D MKL_LIBRARY_DIRS:PATH="$MKL_DIR/lib/intel64" \
 -D TPL_ENABLE_LAPACK:BOOL=ON \
 -D LAPACK_INCLUDE_DIRS:PATH="$MKL_DIR/include" \
 -D LAPACK_LIBRARY_DIRS:PATH="$MKL_DIR/lib/intel64" \
 -D LAPACK_LIBRARY_NAMES:STRING="mkl_lapack95_lp64" \
 -D TPL_ENABLE_BLAS:BOOL=ON \
 -D BLAS_INCLUDE_DIRS:PATH="$MKL_DIR/include" \
 -D BLAS_LIBRARY_DIRS:PATH="$MKL_DIR/lib/intel64" \
 -D BLAS_LIBRARY_NAMES:STRING="mkl_blas95_lp64" \
 -D TPL_ENABLE_HYPRE:STRING=ON \
 -D HYPRE_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D HYPRE_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
 -D TPL_ENABLE_SuperLU:STRING=ON \
 -D SuperLU_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D SuperLU_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
 -D SuperLU_LIBRARY_NAMES:STRING="superlu_4.3" \
 -D TPL_ENABLE_SuperLUDist:BOOL=ON \
 -D SuperLUDist_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D SuperLUDist_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
 -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist_3.3" \
 -D TPL_ENABLE_MUMPS:BOOL=ON \
 -D MUMPS_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D MUMPS_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
\
 -D Trilinos_ENABLE_SCOREC:BOOL=ON \
 -D PCU_COMPRESS:BOOL=ON \
 -D SCOREC_DISABLE_STRONG_WARNINGS:BOOL=ON \
 -D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
 -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
 -D TPL_ENABLE_ParMETIS:STRING=ON \
 -D ParMETIS_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D ParMETIS_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
 -D TPL_ENABLE_METIS:STRING=ON \
 -D METIS_INCLUDE_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/include" \
 -D METIS_LIBRARY_DIRS:PATH="$PETSC_DIR/$PETSC_ARCH/lib" \
 -D Zoltan_ENABLE_ULLONG_IDS:BOOL=OFF \
 -D Teuchos_ENABLE_LONG_LONG_INT:BOOL=OFF \
\
 -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
 -D Tpetra_INST_FLOAT=OFF \
 -D Tpetra_INST_INT_INT=ON \
 -D Tpetra_INST_DOUBLE=ON \
 -D Tpetra_INST_COMPLEX_FLOAT=OFF \
 -D Tpetra_INST_COMPLEX_DOUBLE=OFF \
 -D Tpetra_INST_INT_LONG=OFF \
 -D Tpetra_INST_INT_UNSIGNED=OFF \
 -D Tpetra_INST_INT_LONG_LONG=OFF \
\
../
