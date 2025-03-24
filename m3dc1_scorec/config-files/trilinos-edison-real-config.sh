# Modify these paths for your system.
setenv MPIHOME /opt/cray/mpt/7.2.1/gni/mpich2-intel/140
setenv BLAS_DIR /usr/common/usg/gsl/1.16/intel
setenv NETCDF_DIR /opt/cray/netcdf-hdf5parallel/4.3.3.1/INTEL/140
#/opt/cray/netcdf/4.3.3.1/INTEL/140
setenv HDF5_DIR /opt/cray/hdf5-parallel/1.8.14/INTEL/140
setenv BOOST_DIR /usr/common/usg/boost/1.54/intel
setenv TPSL_DIR /global/project/projectdirs/mp288/edison/petsc-3.5.4-real/cray-mpich-7.2 
setenv MKL_DIR /opt/intel/composer_xe_2013.5.192/mkl
setenv PNETCDF_DIR /opt/cray/parallel-netcdf/1.6.0/INTEL/14.0

#/global/project/projectdirs/mp288/edison/cmake-3.3.2/bin/cmake \
cmake \
\
 -D Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON \
 -D CMAKE_INSTALL_PREFIX:PATH="/global/project/projectdirs/mp288/edison/trilinos-12.2.1/real-install" \
 -D CMAKE_BUILD_TYPE:STRING=DEBUG \
 -D TPL_ENABLE_MPI:BOOL=ON \
 -D MPI_BASE_DIR:PATH=$MPIHOME \
 -D CMAKE_C_COMPILER:STRING="/opt/cray/craype/2.2.1/bin/cc" \
 -D CMAKE_CXX_COMPILER:STRING="/opt/cray/craype/2.2.1/bin/CC" \
 -D CMAKE_Fortran_COMPILER:STRING="/opt/cray/craype/2.2.1/bin/ftn" \
 -D CMAKE_C_FLAGS:STRING="-O3 -g " \
 -D CMAKE_CXX_FLAGS:STRING="-O3 -g " \
 -D Trilinos_ENABLE_CXX11:BOOL=ON \
 -D Trilinos_CXX11_FLAGS:STRING="-std=c++11" \
 -D CMAKE_ENABLE_Fortran:BOOL=ON \
 -D CMAKE_Fortran_FLAGS:STRING="-O3 -g -fpic" \
 -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
 -D BUILD_SHARED_LIBS:BOOL=OFF \
 -D TPL_FIND_SHARED_LIBS:BOOL=ON \
 -D Trilinos_LINK_SEARCH_START_STATIC:BOOL=ON \
# -D Trilinos_EXTRA_LINK_FLAGS:STRING="-ldl" \
# -D CMAKE_FIND_LIBRARY_SUFFIXES:STRING=".a" \
# -D CMAKE_EXE_LINKER_FLAGS:STRING="-static -static-libgcc -static-libstdc++" \
\
 -D Trilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=OFF \
\
 -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -D Trilinos_ENABLE_TESTS:BOOL=OFF \
 -D Trilinos_ENABLE_Gtest:BOOL=OFF \
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
 -D TPL_Netcdf_LIBRARIES:FILEPATH="$NETCDF_DIR/lib/libnetcdf.a" \
 -D TPL_ENABLE_Pnetcdf:BOOL=ON \
 -D Pnetcdf_INCLUDE_DIRS:PATH="$PNETCDF_DIR/include" \
 -D Pnetcdf_LIBRARY_DIRS:PATH="$PNETCDF_DIR/lib" \
 -D TPL_ENABLE_HDF5:BOOL=ON \
 -D HDF5_INCLUDE_DIRS:PATH="$HDF5_DIR/include" \
 -D HDF5_LIBRARY_DIRS:PATH="$HDF5_DIR/lib" \
 -D TPL_HDF5_LIBRARIES:FILEPATH="$HDF5_DIR/lib/libhdf5.a" \
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
 -D MKL_LIBRARY_NAMES:STRING="mkl_intel_lp64" \
 -D TPL_ENABLE_LAPACK:BOOL=ON \
 -D LAPACK_INCLUDE_DIRS:PATH="$MKL_DIR/include" \
 -D LAPACK_LIBRARY_DIRS:PATH="$MKL_DIR/lib/intel64" \
 -D LAPACK_LIBRARY_NAMES:STRING="mkl_lapack95_lp64" \
 -D TPL_ENABLE_BLAS:BOOL=ON \
 -D BLAS_INCLUDE_DIRS:PATH="$MKL_DIR/include" \
 -D BLAS_LIBRARY_DIRS:PATH="$MKL_DIR/lib/intel64" \
 -D BLAS_LIBRARY_NAMES:STRING="mkl_blas95_lp64" \
 -D TPL_ENABLE_HYPRE:STRING=OFF \
 -D HYPRE_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D HYPRE_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
 -D TPL_ENABLE_SuperLU:STRING=OFF \
 -D SuperLU_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D SuperLU_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
 -D SuperLU_LIBRARY_NAMES:STRING="superlu_4.3" \
 -D TPL_ENABLE_SuperLUDist:BOOL=OFF \
 -D SuperLUDist_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D SuperLUDist_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
 -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist_3.3" \
 -D TPL_ENABLE_MUMPS:BOOL=ON \
 -D MUMPS_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D MUMPS_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
\
 -D Trilinos_ENABLE_SCOREC:BOOL=ON \
 -D PCU_COMPRESS:BOOL=ON \
 -D SCOREC_DISABLE_STRONG_WARNINGS:BOOL=ON \
 -D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
 -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
 -D TPL_ENABLE_ParMETIS:STRING=ON \
 -D ParMETIS_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D ParMETIS_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
 -D TPL_ParMETIS_LIBRARIES:FILEPATH="$TPSL_DIR/lib/libparmetis.a" \
 -D TPL_ENABLE_METIS:STRING=ON \
 -D METIS_INCLUDE_DIRS:PATH="$TPSL_DIR/include" \
 -D METIS_LIBRARY_DIRS:PATH="$TPSL_DIR/lib" \
 -D TPL_METIS_LIBRARIES:FILEPATH="$TPSL_DIR/lib/libmetis.a" \
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
