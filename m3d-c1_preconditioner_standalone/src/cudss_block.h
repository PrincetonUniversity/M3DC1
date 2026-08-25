// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * cudss_block.h  --  extern "C" interface to a cuDSS LU solve for one plane
 *                    diagonal block, used as the per-plane sub-solver inside
 *                    a PETSc PCShell block-Jacobi preconditioner.
 *
 * Two modes are supported transparently, selected by how many ranks own the
 * plane (i.e. the size of plane_comm):
 *   - plane_comm size == 1  -> single-GPU cuDSS (no comm layer, no Row1d).
 *                              This is Config B (4 ranks / 1 node, 1 plane/rank).
 *   - plane_comm size  > 1  -> cuDSS MGMN with the NCCL comm layer, the plane
 *                              block 1D row-distributed across the plane's
 *                              ranks/GPUs.  This is Config A.
 *
 * The .c driver builds, per local rank, the plane DIAGONAL block in CSR with:
 *   - local rows  = this rank's rows of the plane,
 *   - column indices rebased to PLANE-LOCAL space (0..n_global_plane-1),
 *   - off-plane couplings DROPPED (the block-Jacobi approximation).
 * Column indices are therefore "global within the plane", exactly what the
 * cuDSS Row1d distribution expects.
 */
#ifndef CUDSS_BLOCK_H
#define CUDSS_BLOCK_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CudssBlock CudssBlock;

/*
 * Create the cuDSS block solver and run ANALYSIS + FACTORIZATION.
 *
 *   plane_comm        MPI sub-communicator over the ranks owning this plane.
 *   n_global_plane    number of rows in the whole plane (global plane size).
 *   local_n           number of local rows on this rank.
 *   local_first_row   0-based, plane-relative index of this rank's first row.
 *   h_rowptr          int32[local_n+1], LOCAL CSR offsets (start at 0).
 *   h_colidx          int32[local_nnz], PLANE-LOCAL (global-within-plane) cols.
 *   h_vals            double[local_nnz], matching values.
 *   local_nnz         number of local nonzeros.
 *   comm_lib_path     path to libcudss_commlayer_nccl.so.0 (used only when the
 *                     plane spans >1 rank; may be NULL otherwise).
 *
 * The CUDA device must already be the one PETSc uses for this rank's vectors
 * (the caller does NOT change device); cuDSS runs on the current device.
 */
CudssBlock *cudss_block_create(MPI_Comm plane_comm,
                               int n_global_plane,
                               int local_n,
                               int local_first_row,
                               const int *h_rowptr,
                               const int *h_colidx,
                               const double *h_vals,
                               long long local_nnz,
                               const char *comm_lib_path);

/*
 * Apply M^{-1} for this plane block: solve (block) * sol = rhs.
 *   d_rhs   DEVICE pointer to the local rhs segment (length local_n).
 *   d_sol   DEVICE pointer to the local solution segment (length local_n).
 * Both are on the current device.  Synchronizes the device on return so the
 * result is visible to the caller's (PETSc) stream.
 */
void cudss_block_solve(CudssBlock *blk, const double *d_rhs, double *d_sol);

/*
 * Re-run the NUMERIC factorization only, reusing the symbolic ANALYSIS already
 * computed in cudss_block_create.  Mirrors an M3D-C1 timestep: the sparsity
 * pattern is fixed, only the values change (here we re-submit the same device
 * values, which has identical factorization cost).  Returns elapsed ms.
 */
double cudss_block_refactor(CudssBlock *blk);

/* One-time setup timings and accumulated solve timing (milliseconds). */
void cudss_block_get_timers(CudssBlock *blk,
                            double *analysis_ms,
                            double *factor_ms,
                            double *solve_ms,
                            int *n_solves);

/* Reset only the per-solve accumulators (keep analysis/factor). */
void cudss_block_reset_solve_timer(CudssBlock *blk);

void cudss_block_destroy(CudssBlock *blk);

#ifdef __cplusplus
}
#endif

#endif /* CUDSS_BLOCK_H */
