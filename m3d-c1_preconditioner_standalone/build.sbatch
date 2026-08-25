#!/bin/bash
#SBATCH -A mp288_g
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH --gpus-per-node=1
#SBATCH -c 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:20:00
#SBATCH -J build_petsc_cudss
#SBATCH -o build-%j.out
#
# Build petsc_cudss_solve: PETSc LGMRES (aijcusparse/cuda) with a PCShell
# block-Jacobi preconditioner whose per-plane sub-solver is cuDSS (single-GPU
# for one rank per plane, MGMN+NCCL for several ranks per plane).  Split build:
#
#   nvcc (cudatoolkit/12.9, gcc-14 OK) -c cudss_block.cu  -> cudss_block.o
#   cc   (Cray wrapper + PETSc pkg-config) -c petsc_cudss_solve.c -> .o
#   cc   link both + PETSc static libs + cudss + nccl + cudart/cusparse/cublas
#
# Output: bin/ under the directory this job is submitted from (README §2).
#
# GPU node because nvcc + cuDSS want a real GPU on host for JIT/codegen checks.

set -euo pipefail

# Package root = the directory this script was submitted from (works for any
# user / any location).  Falls back to the script's own dir for bare `bash`.
PROJ_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJ_DIR"
SRC_DIR="$PROJ_DIR/src"
CU_SRC="$SRC_DIR/cudss_block.cu"
C_SRC="$SRC_DIR/petsc_cudss_solve.c"
C_SRC_MG="$SRC_DIR/petsc_cudss_solve_mg.c"
SHIM_SRC="$SRC_DIR/cudss_commlayer_craympi.c"
BIN_DIR="$PROJ_DIR/bin"
OUT="$BIN_DIR/petsc_cudss_solve"
OUT_MG="$BIN_DIR/petsc_cudss_solve_mg"
SHIM_OUT="$BIN_DIR/libcudss_commlayer_craympi.so"
OBJ_DIR="$BIN_DIR/obj"
mkdir -p "$BIN_DIR" "$OBJ_DIR"

PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606
PETSC_ARCH=perlmuttergpu-gnu

# cuDSS 0.8.0.10.  Vendored in deps/cudss so the package is self-contained;
# override with CUDSS_DIR=... (e.g. your own pip/wheel install) if you prefer.
CUDSS_DIR="${CUDSS_DIR:-$PROJ_DIR/deps/cudss}"
CUDSS_INC="$CUDSS_DIR/include"
CUDSS_LIB="$CUDSS_DIR/lib"
[[ -f "$CUDSS_INC/cudss.h"        ]] || { echo "ERROR: cudss.h not found at $CUDSS_INC"; exit 2; }
[[ -f "$CUDSS_LIB/libcudss.so.0"  ]] || { echo "ERROR: libcudss.so.0 not found at $CUDSS_LIB"; exit 2; }

# nvcc/cc -lcudss looks for libcudss.so (unversioned) — wheel doesn't ship it.
LINK_DIR="$BIN_DIR/cudss_link"
mkdir -p "$LINK_DIR"
ln -sf "$CUDSS_LIB/libcudss.so.0" "$LINK_DIR/libcudss.so"

# Modules.  cudatoolkit/12.9 nvcc supports gcc-14 (PrgEnv-gnu/8.5); 12.4 doesn't.
# nccl module needs cudatoolkit loaded first (prereq).
module load craype-x86-milan PrgEnv-gnu/8.5.0
module load cray-hdf5-parallel/1.14.3.7 cray-fftw/3.3.10.11
module swap cray-mpich cray-mpich/8.1.30
module load cudatoolkit/12.9
module load nccl/2.24.3

NCCL_INC="${NCCL_DIR}/include"
NCCL_LIB="${NCCL_DIR}/lib"
[[ -f "$NCCL_INC/nccl.h"        ]] || { echo "ERROR: nccl.h not found at $NCCL_INC"; exit 2; }
[[ -f "$NCCL_LIB/libnccl.so.2"  ]] || { echo "ERROR: libnccl.so.2 not found at $NCCL_LIB"; exit 2; }

# Link cudart/cusparse/cublas against the SAME cuda 12.4 PETSc was built with
# (NVHPC 24.5), NOT the 12.9 nvcc toolkit used to compile the .cu.  Mixing a
# 12.9 cublas/cudart into PETSc's process (built vs 12.4) causes cuBLAS
# EXECUTION_FAILED at the first VecNorm.  nvcc-12.9 objects are minor-version
# compatible with the 12.4 runtime (same cuda major 12).
NVHPC_MATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4
NVHPC_CUDA=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4

# Cray's GPU transport layer is reached through an UNVERSIONED symlink,
# /opt/cray/pe/lib64/libmpi_gtl_cuda.so.0, which follows the newest CPE rather
# than the pinned cray-mpich: since CPE 26.03 it points at mpich/9.1.0, built
# against libcudart.so.13, absent under CUDA 12.4.  Put the GTL matching
# cray-mpich/8.1.30 first, as the run scripts do, so the ldd check below sees
# what the runs will see.
export LD_LIBRARY_PATH="/opt/cray/pe/mpich/8.1.30/gtl/lib:${LD_LIBRARY_PATH:-}"

echo "=== build env ==="
module list 2>&1
echo "nvcc: $(which nvcc)"; nvcc --version | tail -1
echo "cc:   $(which cc)";   cc --version | head -1
echo "CUDSS_DIR=$CUDSS_DIR"
echo "NCCL_DIR=$NCCL_DIR"
echo ""

# ---- 1) compile the cuDSS/NCCL object with nvcc ----
# nvcc is not the Cray wrapper, so it doesn't know the MPICH include path;
# pass it explicitly ($CRAY_MPICH_DIR is set by the cray-mpich module).
echo "=== nvcc -c cudss_block.cu ==="
nvcc -O3 -std=c++14 -arch=sm_80 \
     -I"$CUDSS_INC" \
     -I"$NCCL_INC" \
     -I"$CRAY_MPICH_DIR/include" \
     -c "$CU_SRC" \
     -o "$OBJ_DIR/cudss_block.o"
echo "  -> $OBJ_DIR/cudss_block.o"
echo ""

# ---- 2) compile the PETSc driver with cc (PETSc pkg-config cflags) ----
export PKG_CONFIG_PATH="$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
PETSC_CFLAGS=$(pkg-config --cflags PETSc)
PETSC_LIBS=$(pkg-config --libs --static PETSc)
echo "=== cc -c petsc_cudss_solve.c ==="
echo "  PETSc cflags: $PETSC_CFLAGS"
cc -O2 $PETSC_CFLAGS -I"$SRC_DIR" -c "$C_SRC" -o "$OBJ_DIR/petsc_cudss_solve.o"
echo "  -> $OBJ_DIR/petsc_cudss_solve.o"
echo ""

# ---- 3) link ----
echo "=== cc link -> $OUT ==="
cc -O2 -o "$OUT" \
   "$OBJ_DIR/petsc_cudss_solve.o" "$OBJ_DIR/cudss_block.o" \
   -Wl,-rpath,"$PETSC_DIR/$PETSC_ARCH/lib" \
   -Wl,-rpath,"$CUDSS_LIB" \
   -Wl,-rpath,"$NCCL_LIB" \
   -Wl,-rpath,"$NVHPC_MATH/lib64" \
   -Wl,-rpath,"$NVHPC_CUDA/lib64" \
   -L"$LINK_DIR" -L"$NCCL_LIB" -L"$NVHPC_MATH/lib64" -L"$NVHPC_CUDA/lib64" \
   $PETSC_LIBS \
   -lcudss -lnccl -lcusparse -lcublas -lcudart -lstdc++ -lm

echo ""

# ---- 3b) compile + link the MG variant (single-process multi-GPU per plane) ----
# Same cudss_block.o; only the driver differs (adds -mg_devices -> CUDSS_MG_DEVICES).
echo "=== cc -c petsc_cudss_solve_mg.c ==="
cc -O2 $PETSC_CFLAGS -I"$SRC_DIR" -c "$C_SRC_MG" -o "$OBJ_DIR/petsc_cudss_solve_mg.o"
echo "  -> $OBJ_DIR/petsc_cudss_solve_mg.o"
echo "=== cc link -> $OUT_MG ==="
cc -O2 -o "$OUT_MG" \
   "$OBJ_DIR/petsc_cudss_solve_mg.o" "$OBJ_DIR/cudss_block.o" \
   -Wl,-rpath,"$PETSC_DIR/$PETSC_ARCH/lib" \
   -Wl,-rpath,"$CUDSS_LIB" \
   -Wl,-rpath,"$NCCL_LIB" \
   -Wl,-rpath,"$NVHPC_MATH/lib64" \
   -Wl,-rpath,"$NVHPC_CUDA/lib64" \
   -L"$LINK_DIR" -L"$NCCL_LIB" -L"$NVHPC_MATH/lib64" -L"$NVHPC_CUDA/lib64" \
   $PETSC_LIBS \
   -lcudss -lnccl -lcusparse -lcublas -lcudart -lstdc++ -lm
echo "  -> $OUT_MG"
echo ""

# ---- 4) build the custom cray-mpich cuDSS comm-layer shim ----
# The Cray `cc` wrapper provides GPU-aware cray-mpich + MPI includes/libs; only
# add the cuDSS dist-iface header and the CUDA include (cuda_runtime.h,
# library_types.h).  No explicit -lcudart: cudaStreamSynchronize/cudaGetErrorString
# resolve against the cudart already loaded by the host process at dlopen time
# (keeps the shim agnostic to the 12.4-vs-12.9 cudart the host links).
echo "=== cc -shared cudss_commlayer_craympi.c -> $SHIM_OUT ==="
cc -O2 -fPIC -shared \
   -I"$CUDSS_INC" \
   -I"$CUDA_HOME/include" \
   "$SHIM_SRC" \
   -o "$SHIM_OUT"
echo "  -> $SHIM_OUT"
echo ""

echo "=== built ==="
ls -lh "$OUT" "$OUT_MG" "$SHIM_OUT"
file "$OUT"
echo ""
echo "=== ldd (verify cudss/nccl resolve) ==="
ldd "$OUT" | grep -Ei 'cudss|nccl|cudart|cusparse|petsc' || true
echo ""
echo "=== shim: exported symbol + ldd (verify cray-mpich libmpi resolves) ==="
nm -D "$SHIM_OUT" | grep -i cudssDistributedInterface || true
ldd "$SHIM_OUT" | grep -Ei 'mpi|cudart|cuda' || true
