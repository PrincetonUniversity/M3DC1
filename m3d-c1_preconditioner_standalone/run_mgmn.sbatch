#!/bin/bash
#SBATCH -A mp288_g
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -c 32
#SBATCH -t 00:30:00
#SBATCH -o slurm-mgmn-%j.out
#
# SEVERAL RANKS PER PLANE (MGMN): 8 ranks / 2 nodes, 4 planes -> 2 ranks per
# plane, one GPU per rank, so each plane block is factored by cuDSS distributed
# across two GPUs through a comm layer.  See README §3.3.
#
# Copied from run_config_b.sbatch; the only changes are the node/rank count
# (nodes 1 -> 2, srun -n 4 -> 8), the output file names, and the comm-layer
# notes below.  -nplanes stays 4 because the matrix has 4 planes.
#
# COMM LAYER: this uses the Cray-MPICH shim built from
# src/cudss_commlayer_craympi.c, because the NCCL comm layer that ships inside
# cuDSS is linked against OpenMPI (libmpi.so.40) and therefore cannot load on
# Perlmutter -- cudssSetCommLayer then fails with CUDSS_STATUS_INVALID_VALUE (3).
# The path must be absolute: cuDSS dlopen()s it and rejects a relative path.
#
# NOTE ON CONVERGENCE: the distributed path is more restricted than the
# single-GPU one.  cuDSS does not support matching under MGMN -- with matching
# requested, ANALYSIS returns CUDSS_STATUS_NOT_SUPPORTED (4) -- so this script
# turns it off.  Matching is what puts large entries on the diagonal for this
# matrix, and without it the factor is inaccurate and the outer LGMRES needs
# orders of magnitude more iterations than one rank per plane does.  Hence this
# script exercises the distributed path; it is not the fast configuration.
# -ksp_max_it 10000 bounds the worst case.
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

LOG_VIEW=":profile_mgmn.txt,:profile_mgmn.fg:ascii_flamegraph,:profile_mgmn.xml:ascii_xml"

export CUDSS_MGMN_MATCH=0
export CUDSS_COMM_BACKEND=mpi
COMM_LIB=(-cudss_comm_lib "$PWD/bin/libcudss_commlayer_craympi.so")
[[ -f "${COMM_LIB[1]}" ]] || { echo "ERROR: ${COMM_LIB[1]} not found; run build.sbatch"; exit 2; }

echo "=== petsc_cudss MGMN: 8 ranks / 2 nodes, 2 ranks per plane, -log_view $LOG_VIEW ==="
echo "Binary: $BINARY"
echo "Dir:    $(pwd)"
echo "Nodes:  $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS  GPUs/node: 4"
echo "Comm:   ${CUDSS_COMM_BACKEND:-nccl} ${COMM_LIB[1]}"
echo ""

srun -n 8 -c 32 --cpu-bind=cores --gpu-bind=none --export=ALL \
    "$BINARY" \
    -A data/A.bin -b data/b.bin -x0 data/x0.bin -xref data/xout.bin \
    -nplanes 4 -nreps 2 "${COMM_LIB[@]}" \
    -mat_type aijcusparse -vec_type cuda \
    -ksp_type lgmres -ksp_lgmres_augment 4 -ksp_gmres_restart 220 \
    -ksp_max_it 10000 -ksp_rtol 1e-9 -ksp_atol 1e-20 \
    -ksp_converged_reason \
    -log_view "$LOG_VIEW"
