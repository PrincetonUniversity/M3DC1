# M3DC1 Build Instructions

## Build Options

The essential build options are `M3DC1_ENABLE_3D`, `M3DC1_ENABLE_COMPLEX`, `M3DC1_ENABLE_PARTICLE`, `M3DC1_ENABLE_GPU`, `M3DC1_ENABLE_ST`, and `M3DC1_ENABLE_ADAS`. They are all **OFF** by default. To enable them, set the desired option to `ON`.

## Essential Dependencies

- GCC (C, C++, Fortran compilers)
- MPI (e.g., MPICH)
- HDF5 (with Fortran and HL support)
- CMake (>= 3.20)
- Ninja (optional, recommended)
- GSL
- FFTW
- NetCDF-C (with MPI support)
- PUMI (>= 2.2.8, with Zoltan)
- PETSc (>= 3.19, with MUMPS; add `+complex` for complex number support)

---

## 1. Building with Spack

One way to configure the dependencies is through [Spack](https://spack.io/). All essential dependencies can be installed via a Spack environment. An example `spack.yaml`:

```yaml
spack:
  specs:
  - gcc@14.2.0
  - mpich %gcc@14.2.0
  - hdf5+fortran+hl %gcc@14.2.0
  - ninja %gcc@14.2.0
  - cmake %gcc@14.2.0
  - gsl %gcc@14.2.0
  - fftw %gcc@14.2.0
  - netcdf-c+mpi %gcc@14.2.0
  - pumi@2.2.8+zoltan ^zoltan+parmetis %gcc@14.2.0
  - petsc ~hypre+mumps %gcc@14.2.0
  # Use petsc +complex if complex number support is needed
```

> **Note:** If you need complex number support, install PETSc with the `+complex` variant.

Once the Spack environment is activated, build M3DC1:

```bash
cmake -S M3DC1 -B build-m3dc1 \
  -DCMAKE_INSTALL_PREFIX=build-m3dc1/install \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DM3DC1_ENABLE_3D=OFF \
  -DM3DC1_ENABLE_COMPLEX=OFF \
  -DM3DC1_ENABLE_PARTICLE=OFF \
  -DM3DC1_ENABLE_GPU=OFF \
  -DM3DC1_ENABLE_ST=OFF \
  -DM3DC1_ENABLE_ADAS=OFF \
  -DCMAKE_BUILD_TYPE=Debug

cmake --build build-m3dc1 -j 8
```

---

## 2. Building on Scorec

On the Scorec system, most dependencies are available through environment modules. Load the required modules:

```bash
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno
module load mpich/4.1.1-xpoyz4t
module load cmake/3.26.3-2duxfcd
module load cuda/12.1.1-zxa4msk
module load parmetis/4.0.3-yyczvvl
module load openblas
module load ninja
module load gsl
module load fftw
```

### 2.1 Build PETSc, HDF5, Zoltan, and Core

The remaining dependencies (PETSc, HDF5, Zoltan, and Core) must be built manually:

```bash
export PARMETIS_HOME=/opt/scorec/spack/rhel9/v0201_4/install/linux-rhel9-x86_64/gcc-12.3.0/parmetis-4.0.3-yyczvvlsvs5skddbpr5vr4z63fkdt3ks
export METIS_HOME=/opt/scorec/spack/rhel9/v0201_4/install/linux-rhel9-x86_64/gcc-12.3.0/metis-5.1.0-65szzoyrtgauis34eop3w5zu6uarer
export PETSC_DIR=<path-to-petsc-source>
export PETSC_ARCH=cuda
export PREFIX=/users/yus9/lore.scorec.rpi.edu/issues_test/m3dc1_test

# Build PETSc
cd $PETSC_DIR
./configure \
  PETSC_ARCH=$PETSC_ARCH \
  --with-parmetis-dir=$PARMETIS_HOME \
  --with-metis-dir=$METIS_HOME \
  --with-cuda=1 \
  --with-shared-libraries=0 \
  --download-mumps \
  --download-scalapack \
  --with-openblas-dir="${OPENBLAS_RHEL9_ROOT}"
make all check
cd -

# Build HDF5
cmake -S hdf5 -B build-hdf5 \
  -DCMAKE_INSTALL_PREFIX=$PWD/build-hdf5/install \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_BUILD_HL_LIB=ON \
  -DHDF5_ENABLE_PARALLEL=ON \
  -DHDF5_ENABLE_ZLIB_SUPPORT=ON \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc
cmake --build build-hdf5 -j 8 --target install

# Build Zoltan (from Trilinos)
cmake -S Trilinos -B build-zoltan \
  -GNinja \
  -DCMAKE_INSTALL_PREFIX=build-zoltan/install \
  -DTPL_ENABLE_MPI:BOOL=ON \
  -DCMAKE_C_FLAGS="-O3 -fPIC" \
  -DCMAKE_CXX_FLAGS="-O3 -fPIC" \
  -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -DHDF5_ROOT=$PWD/build-hdf5/install \
  -DTrilinos_ENABLE_Zoltan:BOOL=ON \
  -DZoltan_ENABLE_EXAMPLES:BOOL=OFF \
  -DZoltan_ENABLE_TESTS:BOOL=OFF \
  -DZoltan_ENABLE_ParMETIS:BOOL=ON \
  -DParMETIS_INCLUDE_DIRS=$PARMETIS_HOME/include \
  -DParMETIS_LIBRARY_DIRS=$PARMETIS_HOME/lib
cmake --build build-zoltan -j 4 --target install

# Build Core
rm -rf build-core
cmake -S core -B build-core \
  -DCMAKE_C_COMPILER="$MPIHOME/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$MPIHOME/bin/mpicxx" \
  -DCMAKE_C_FLAGS="-g -O2" \
  -DCMAKE_CXX_FLAGS="-g -O2" \
  -DENABLE_THREADS=ON \
  -DENABLE_ZOLTAN=ON \
  -DCMAKE_INSTALL_PREFIX="build-core/install" \
  -DZOLTAN_INCLUDE_DIR="build-zoltan/install/include" \
  -DZOLTAN_LIBRARY="build-zoltan/install/lib64/libzoltan.a" \
  -DPARMETIS_PREFIX=$PARMETIS_HOME/lib \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-core -j 4 --target install
```

### 2.2 Build M3DC1

Since dependencies are built from scratch, CMake paths must be set manually:

```bash
cmake -S M3DC1 -B build-m3dc1 \
  -DCMAKE_INSTALL_PREFIX=build-m3dc1/install \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DPETSC_LINK_STATIC=ON \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DPETSC_LIBRARIES=$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so \
  -DCMAKE_PREFIX_PATH="$PWD/build-core/install;$PWD/build-hdf5/install" \
  -DHDF5_ROOT=$PWD/build-hdf5/install \
  -DHDF5_INCLUDE_DIRS=$PWD/build-hdf5/install/include \
  -DHDF5_LIBRARIES=$PWD/build-hdf5/install/lib/libhdf5.so \
  -DM3DC1_ENABLE_3D=OFF \
  -DM3DC1_ENABLE_COMPLEX=OFF \
  -DM3DC1_ENABLE_PARTICLE=OFF \
  -DM3DC1_ENABLE_GPU=ON \
  -DM3DC1_ENABLE_ST=OFF \
  -DM3DC1_ENABLE_ADAS=ON \
  -DCMAKE_BUILD_TYPE=Debug

cmake --build build-m3dc1 -j 8
```

---

## 3. Building on Perlmutter

On the Perlmutter system, most dependencies are available through environment modules. Load the required modules first, then follow the steps below.

### 3.1 Build ParMETIS (includes METIS and GKlib)

```bash
cd GKlib
make config prefix=$PWD/GKlib/build
make install
cd ..

cd METIS
make config prefix=$PWD/METIS/build gklib_path=$PWD/GKlib/build
make install
cd ..

cd ParMETIS
make config cc=mpicc prefix=$PWD/ParMETIS/build gklib_path=$PWD/GKlib/build metis_path=$PWD/METIS/build
make install
cd ..
```

### 3.2 Build PETSc, HDF5, Zoltan, and Core

```bash
# Build PETSc
cd petsc
./configure \
  PETSC_ARCH=cuda \
  --with-cuda=1 \
  --with-shared-libraries=0 \
  --download-mumps \
  --download-scalapack \
  --download-openblas=1
# use --with-scalar-type=complex for complex number support
make all
cd ..

# Build HDF5
cmake -S hdf5 -B build-hdf5 \
  -DCMAKE_INSTALL_PREFIX=$PWD/build-hdf5/install \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_BUILD_HL_LIB=ON \
  -DHDF5_ENABLE_PARALLEL=ON \
  -DHDF5_ENABLE_ZLIB_SUPPORT=ON \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc
cmake --build build-hdf5 -j 8 --target install

# Build Zoltan (from Trilinos)
cmake -S Trilinos -B build-zoltan \
  -GNinja \
  -DCMAKE_INSTALL_PREFIX=build-zoltan/install \
  -DTPL_ENABLE_MPI:BOOL=ON \
  -DCMAKE_C_FLAGS="-O3 -fPIC" \
  -DCMAKE_CXX_FLAGS="-O3 -fPIC" \
  -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -DHDF5_ROOT=$PWD/build-hdf5/install \
  -DTrilinos_ENABLE_Zoltan:BOOL=ON \
  -DZoltan_ENABLE_EXAMPLES:BOOL=OFF \
  -DZoltan_ENABLE_TESTS:BOOL=OFF \
  -DZoltan_ENABLE_ParMETIS:BOOL=ON \
  -DParMETIS_INCLUDE_DIRS="$PWD/ParMETIS/build/include;$PWD/METIS/build/include" \
  -DParMETIS_LIBRARY_DIRS="$PWD/ParMETIS/build/lib;$PWD/METIS/build/lib"
cmake --build build-zoltan -j 4 --target install

# Build Core
cmake -S core -B build-core \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_FLAGS="-g -O2" \
  -DCMAKE_CXX_FLAGS="-g -O2" \
  -DENABLE_THREADS=ON \
  -DENABLE_ZOLTAN=ON \
  -DCMAKE_INSTALL_PREFIX="build-core/install" \
  -DZOLTAN_INCLUDE_DIR="build-zoltan/install/include" \
  -DZOLTAN_LIBRARY="build-zoltan/install/lib64/libzoltan.a" \
  -DPARMETIS_PREFIX=$PWD/ParMETIS/build \
  -DMETIS_INCLUDE_DIR=$PWD/METIS/build/include \
  -DMETIS_LIBRARY=$PWD/METIS/build/lib/libmetis.so \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-core -j 4 --target install
```

### 3.3 Build M3DC1

```bash
cmake -S M3DC1 -B build-m3dc1 \
  -DCMAKE_INSTALL_PREFIX=build-m3dc1/install \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DPETSC_LINK_STATIC=ON \
  -DPETSC_DIR=$PWD/petsc \
  -DPETSC_ARCH=cuda \
  -DPETSC_LIBRARIES=$PWD/petsc/cuda/lib/libpetsc.a \
  -DCMAKE_PREFIX_PATH="$PWD/build-core/install;$PWD/build-hdf5/install;/opt/cray/pe/fftw/3.3.10.11/x86_milan;$PWD/GKlib/build/lib:$PWD/METIS/build/lib:$PWD/ParMETIS/build/lib" \
  -DHDF5_ROOT=$PWD/build-hdf5/install \
  -DHDF5_INCLUDE_DIRS=$PWD/build-hdf5/install/include \
  -DHDF5_LIBRARIES=$PWD/build-hdf5/install/lib/libhdf5.so \
  -DM3DC1_ENABLE_3D=OFF \
  -DM3DC1_ENABLE_COMPLEX=OFF \
  -DM3DC1_ENABLE_PARTICLE=OFF \
  -DM3DC1_ENABLE_GPU=ON \
  -DM3DC1_ENABLE_ST=OFF \
  -DM3DC1_ENABLE_ADAS=ON \
  -DCMAKE_BUILD_TYPE=Debug

cmake --build build-m3dc1 -j 8
```


4. ## 4. Building on Stellar
The following modules are required to build M3DC1 on Stellar:

```
module purge
module load \
  intel/2021.1.2 \
  intel-mpi/intel/2021.3.1 \
  fftw/intel-2021.1/intel-mpi/3.3.9 \
  hdf5/intel-2021.1/intel-mpi/1.10.6 \
  netcdf/intel-2021.1/hdf5-1.10.6/intel-mpi/4.7.4 \
  gsl/2.6
```

Then build M3DC1 with the following commands:

```
PETSC_DIR=/home/ur8212/sourceCodes/M3DC1/dependencies/petsc
PETSC_ARCH=real-intel2021.1.2-intelmpi2021.3.1
HDF5_DIR=/usr/local/hdf5/intel-2021.1/intel-mpi/1.10.6
CORE_DIR=/home/ur8212/sourceCodes/M3DC1/dependencies/build-core
NETCDF_DIR=/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/intel-mpi/4.7.4 
BUILD_DIR=/projects/M3DC1/PETSC/petsc-3.13.5/real-intel2021.1.2-intelmpi2021.3.1
PARMETIS_DIR=$BUILD_DIR
ZOLTAN_DIR=$BUILD_DIR
FFTW_DIR=/usr/local/fftw/intel-2021.1/intel-mpi/3.3.9
cmake -S . -B build-m3dc1 \
  -DCMAKE_INSTALL_PREFIX=build-m3dc1/install \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_FLAGS="-g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS="-g -O0 -DPETSCMASTER -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-r8 -fpp" \
  -DPETSC_LINK_STATIC=ON \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DPETSC_LIBRARIES=$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so \
  -DCMAKE_PREFIX_PATH="$CORE_DIR/install;$HDF5_DIR;$PARMETIS_DIR;$ZOLTAN_DIR" \
  -DHDF5_ROOT=$HDF5_DIR \
  -DHDF5_INCLUDE_DIRS=$HDF5_DIR/include \
  -DHDF5_LIBRARIES=$HDF5_DIR/lib/libhdf5.so \
  -DFFTW3_INCLUDE_DIR=$FFTW_DIR/include \
  -DFFTW3_LIBRARY=$FFTW_DIR/lib/libfftw3.so \
  -DNetCDF_ROOT=$NETCDF_DIR \
  -DM3DC1_ENABLE_3D=OFF \
  -DM3DC1_ENABLE_COMPLEX=OFF \
  -DM3DC1_ENABLE_PARTICLE=OFF \
  -DM3DC1_ENABLE_GPU=OFF \
  -DM3DC1_ENABLE_ST=OFF \
  -DM3DC1_ENABLE_ADAS=ON \
  -DCMAKE_BUILD_TYPE=Debug

cmake --build build-m3dc1 -j 8
```