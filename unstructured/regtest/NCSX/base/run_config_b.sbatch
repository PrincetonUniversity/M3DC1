#!/bin/bash
#SBATCH -A mp288_g
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=4
#SBATCH -c 16
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

BINARY="/global/cfs/cdirs/mp288/jinchen/M3DC1/gitrepo/M3DC1/unstructured/_perlmutter_gpu_gcc-3d-st-opt-125/m3dc1_3d_st"
[[ -x "$BINARY" ]] || { echo "ERROR: $BINARY not found in $(pwd)"; exit 2; }

#module load craype-x86-milan PrgEnv-gnu/8.5.0
#module load cray-hdf5-parallel/1.14.3.7 cray-fftw/3.3.10.11
#module load cudatoolkit/12.4
#module swap cray-mpich cray-mpich/8.1.30
#module load nccl/2.24.3

# PETSc 3.23 GPU shared lib was linked against libnvToolsExt.so.1, removed
# from CUDA >=12.6; prepend NVHPC 24.5 / cuda 12.4 lib dir (still has it).
#export LD_LIBRARY_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/12.4/lib64:${LD_LIBRARY_PATH:-}"
# cuDSS runtime lib (0.8.0.10) + its NCCL comm layer.
CUDSS_LIB="${CUDSS_DIR:-/global/cfs/cdirs/mp288/jinchen/LIB/deps/cudss}/lib"
export LD_LIBRARY_PATH="$CUDSS_LIB:${NCCL_DIR}/lib:${LD_LIBRARY_PATH}"
# Cray's GPU transport layer is reached through an UNVERSIONED symlink,
# /opt/cray/pe/lib64/libmpi_gtl_cuda.so.0, which follows the newest CPE rather
# than the pinned cray-mpich: since CPE 26.03 it points at mpich/9.1.0, built
# against libcudart.so.13, which does not exist under CUDA 12.4 (the binary then
# dies with "libcudart.so.13: cannot open shared object file").  Put the GTL that
# matches cray-mpich/8.1.30 (libcudart.so.12) first.
#export LD_LIBRARY_PATH="/opt/cray/pe/mpich/8.1.30/gtl/lib:${LD_LIBRARY_PATH}"

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND="cores"

LOG_VIEW=":profile.txt,:profile.fg:ascii_flamegraph,:profile.xml:ascii_xml"

echo "=== petsc_cudss Config B: 4 ranks / 1 node, single-GPU cuDSS per plane, -log_view $LOG_VIEW ==="
echo "Binary: $BINARY"
echo "Dir:    $(pwd)"
echo "Nodes:  $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS  GPUs/node: 4"
echo ""

#  - cuDSS defaults to non-deterministic scheduling (faster). 
#    It only becomes bitwise-reproducible when CUDSS_DETERMINISTIC=1 is set — 
#    the env var is read at cudss_block_create time and passed to CUDSS_CONFIG_DETERMINISTIC_MODE (cudss_block.cu:235-238).
#  - If the var is unset (or =0), you get the default non-deterministic mode, 
#    where the recomputed true residual ||Ax-b|| differs in its last digits 
#    between otherwise-identical solves (per the comment at lines 231-234).
#  - It costs performance, so it's opt-in by design.
#  - Confirm it took effect by looking for this line on stderr: [cudss_block] deterministic_mode: 1.
#  One caveat: this sets the cuDSS config flag, but full run-to-run bitwise reproducibility for the whole solve 
#  also depends on the comm backend. NCCL reductions (the MGMN default) are not guaranteed deterministic across runs. 
#  For a single-GPU plane (Config B, plane_size==1, no NCCL) CUDSS_DETERMINISTIC=1 alone should give you reproducible factor/solve.
#export CUDSS_DETERMINISTIC=1
#-mat_cusparse_spmv_alg CSRMV_ALG2
#  a 15% slowdown


srun -n $SLURM_NTASKS -c 16 --cpu-bind=cores --gpu-bind=none --export=ALL \
	"$BINARY" \
	-cudsssolve 1 \
	-mat_type aijcusparse -vec_type cuda \
	-ksp_monitor_true_residual \
	-log_view "$LOG_VIEW" -log_view_gpu_time \
	-ksp_view \
	-options_left

mv C1ke "C1ke.${SLURM_JOB_ID}"

