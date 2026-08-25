#!/bin/bash
#
# srun launcher for cuDSS single-process multi-GPU (MG) mode: give each rank on
# this node a DISJOINT set of $MG_DEVICES GPUs, then exec the real binary.
#
#   srun ... ./mg_gpu_wrap.sh ./bin/petsc_cudss_solve_mg -mg_devices 2 ...
#
# With -mg_devices N a rank factors its plane block across the N GPUs visible to
# it, so two ranks on the same node must not see the same GPU.  Used by
# run_mg.sbatch (README §3.2).  MG_DEVICES defaults to 1, which is a no-op.

set -euo pipefail

N="${MG_DEVICES:-1}"
LOCAL="${SLURM_LOCALID:-0}"
AVAIL="${SLURM_GPUS_ON_NODE:-4}"

first=$(( N * LOCAL ))
last=$(( first + N - 1 ))

if (( last >= AVAIL )); then
    echo "ERROR: rank ${SLURM_PROCID:-?} (localid $LOCAL) needs GPUs $first-$last," \
         "but this node has $AVAIL. Reduce --ntasks-per-node or MG_DEVICES so that" \
         "MG_DEVICES * ntasks-per-node <= $AVAIL." >&2
    exit 2
fi

export CUDA_VISIBLE_DEVICES="$(seq -s, "$first" "$last")"
exec "$@"
