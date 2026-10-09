# Configure the m3dc1_cudss library (CMakeLists.txt in this directory) on
# stellar-ai (Princeton).  Usage, perlmutter.sh-style, from a build directory:
#
#   module load gcc/14
#   module load openmpi/cuda-13.2/gcc/5.0.10   # CUDA-aware MPI (required for cuDSS MGMN)
#   module load cudatoolkit/13.2
#   # NCCL: no cluster module assumed; point NCCL_DIR at an install providing
#   # include/nccl.h and lib/libnccl.so (used by CMakeLists.txt)
#   export NCCL_DIR=/path/to/nccl
#
#   mkdir build && cd build
#   sh ../stellarai.sh
#   make
#
# Differences from perlmutter.sh:
#   - mpicc/mpicxx wrappers (PETSc petsc.20220201 was built with mpicc/mpicxx/mpif90),
#     not Cray cc/CC.
#   - Single CUDA toolkit (13.2) for both compile and runtime.  On Perlmutter two
#     toolkits were mixed (nvcc 12.9 to compile because 12.4's nvcc rejected gcc-14,
#     CUDA 12.4 libs at runtime to match PETSc).  Here nvcc 13.2 accepts gcc-14 and
#     PETSc arch real-stellarai-gcc14-cuda100-gO2 was linked against the same 13.2
#     cudart/cusparse/cublas, so no mixing -- and none allowed: keep 13.4 (or any
#     other toolkit) out of PATH/LD_LIBRARY_PATH or VecNorm dies with
#     cuBLAS EXECUTION_FAILED, as seen on Perlmutter.
#   - CUDA arch 100 (Blackwell, from the "cuda100" in PETSC_ARCH) instead of 80 (A100).
#   - MPI include for nvcc derived from the mpicc wrapper location instead of
#     $CRAY_MPICH_DIR.
#   - cuDSS MGMN comm layer: build deps/cudss/src/cudss_commlayer_openmpi.cu against
#     this OpenMPI (see cudss_build_commlayer.sh), not cudss_commlayer_craympi.c.

CMAKETYPE=Release
PETSC_DIR=/projects/M3DC1/PETSC/petsc.20220201
PETSC_ARCH=cplx-stellarai-gcc14-cuda100-gO2

MPI_HOME=${MPI_HOME:-$(dirname "$(dirname "$(which mpicc)")")}

cmake .. $COMPLEX_OPT \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_CUDA_COMPILER="nvcc" \
  -DCMAKE_CUDA_ARCHITECTURES=100 \
  -DCMAKE_CUDA_FLAGS=" -I$MPI_HOME/include" \
  -DCMAKE_C_FLAGS=" -fPIC -O2 -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -fPIC -O2 -I$PETSC_DIR/include" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DCUDSS_DIR="$(cd .. && pwd)/deps/cudss" \
  -DCMAKE_INSTALL_PREFIX="$(cd .. && pwd)" \
  -DCMAKE_BUILD_TYPE=$CMAKETYPE \
  -DENABLE_COMPLEX=ON
