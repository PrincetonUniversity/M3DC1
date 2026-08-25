#!/bin/bash
#SBATCH -A mp288_g
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -c 32
#SBATCH -t 00:30:00
#SBATCH -o slurm-%j.out
#
# One GPU per plane ("Config B"): PETSc LGMRES (aijcusparse/cuda) with a cuDSS
# block-Jacobi PCShell, 4 ranks / 1 node, 1 plane per rank -> each plane block
# is a SINGLE-GPU cuDSS LU solve (no comm layer, no NCCL).  See README §3.1.
#
# Must be submitted from the package root: the paths below (bin/, data/) are
# relative to $SLURM_SUBMIT_DIR.  There is no -pc_type option because the
# PCShell installs the cuDSS block-Jacobi in code.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

BINARY="./bin/petsc_cudss_solve"
[[ -x "$BINARY" ]] || { echo "ERROR: $BINARY not found in $(pwd)"; exit 2; }

module load craype-x86-milan PrgEnv-gnu/8.5.0
module load cray-hdf5-parallel/1.14.3.7 cray-fftw/3.3.10.11
module load cudatoolkit/12.4
module swap cray-mpich cray-mpich/8.1.30
module load nccl/2.24.3

# PETSc 3.23 GPU shared lib was linked against libnvToolsExt.so.1, removed
# from CUDA >=12.6; prepend NVHPC 24.5 / cuda 12.4 lib dir (still has it).
export LD_LIBRARY_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64:${LD_LIBRARY_PATH:-}"
# cuDSS runtime lib (0.8.0.10) + its NCCL comm layer.
CUDSS_LIB="${CUDSS_DIR:-$PWD/deps/cudss}/lib"
export LD_LIBRARY_PATH="$CUDSS_LIB:${NCCL_DIR}/lib:${LD_LIBRARY_PATH}"
# Cray's GPU transport layer is reached through an UNVERSIONED symlink,
# /opt/cray/pe/lib64/libmpi_gtl_cuda.so.0, which follows the newest CPE rather
# than the pinned cray-mpich: since CPE 26.03 it points at mpich/9.1.0, built
# against libcudart.so.13, which does not exist under CUDA 12.4 (the binary then
# dies with "libcudart.so.13: cannot open shared object file").  Put the GTL that
# matches cray-mpich/8.1.30 (libcudart.so.12) first.
export LD_LIBRARY_PATH="/opt/cray/pe/mpich/8.1.30/gtl/lib:${LD_LIBRARY_PATH}"

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND="cores"

LOG_VIEW=":profile.txt,:profile.fg:ascii_flamegraph,:profile.xml:ascii_xml"

echo "=== petsc_cudss Config B: 4 ranks / 1 node, single-GPU cuDSS per plane, -log_view $LOG_VIEW ==="
echo "Binary: $BINARY"
echo "Dir:    $(pwd)"
echo "Nodes:  $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS  GPUs/node: 4"
echo ""

srun -n 4 -c 32 --cpu-bind=cores --gpu-bind=none --export=ALL \
    "$BINARY" \
    -A data/A.bin -b data/b.bin -x0 data/x0.bin -xref data/xout.bin \
    -nplanes 4 -nreps 2 \
    -mat_type aijcusparse -vec_type cuda \
    -ksp_type lgmres -ksp_lgmres_augment 4 -ksp_gmres_restart 220 \
    -ksp_max_it 10000 -ksp_rtol 1e-9 -ksp_atol 1e-20 \
    -ksp_converged_reason \
    -log_view "$LOG_VIEW"
