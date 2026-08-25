# Configure the m3dc1_cudss library (CMakeLists.txt in this directory) on
# Perlmutter.  Usage, kairos.sh-style, from a build directory:
#
#   module load craype-x86-milan PrgEnv-gnu/8.5.0
#   module swap cray-mpich cray-mpich/8.1.30
#   module load cudatoolkit/12.9
#   module load nccl/2.24.3          # sets NCCL_DIR (used by CMakeLists.txt)
#
#   mkdir build && cd build
#   sh ../perlmutter.sh
#   make
#
# Same toolchain as build.sbatch: nvcc from cudatoolkit/12.9 for cudss_block.cu,
# Cray cc/CC wrappers for the C/C++ side, cuDSS vendored in deps/cudss.

CMAKETYPE=Release
PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606
PETSC_ARCH=perlmuttergpu-gnu

cmake .. \
  -DCMAKE_C_COMPILER="cc" \
  -DCMAKE_CXX_COMPILER="CC" \
  -DCMAKE_CUDA_COMPILER="nvcc" \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_CUDA_FLAGS=" -I$CRAY_MPICH_DIR/include" \
  -DCMAKE_C_FLAGS=" -fPIC -O2 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -fPIC -O2 -I$PETSC_DIR/include" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DCUDSS_DIR="$(cd .. && pwd)/deps/cudss" \
  -DCMAKE_INSTALL_PREFIX="$(cd .. && pwd)" \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE \
