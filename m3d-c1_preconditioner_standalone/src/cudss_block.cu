// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * cudss_block.cu  --  cuDSS LU solve for one plane diagonal block.
 *
 * extern "C" wrappers compiled with nvcc and linked into the PETSc driver.
 * Single-GPU (Config B) when the plane is owned by one rank; cuDSS MGMN with
 * the NCCL comm layer (Config A) when the plane spans >1 rank/GPU.
 *
 * cuDSS config mirrors the proven hybrid_lgmres_left_cudss.cu proxy:
 *   pivot = CUDSS_PIVOT_NONE, matching ENABLED, reorder = CUDSS_ALG_DEFAULT,
 *   pivot_epsilon = 1e-12, ir_steps = 0.
 * (MGMN forbids CUDSS_ALG_1/ALG_2 reordering -- DEFAULT is what we use anyway.)
 */
 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>

#include <cuda_runtime.h>
#include <cudss.h>
#include <nccl.h>

#include "cudss_block.h"

#define CUDA_CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){fprintf(stderr,"[cudss_block] CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)
#define CUDSS_CHECK(call) do { cudssStatus_t s=(call); if(s!=CUDSS_STATUS_SUCCESS){fprintf(stderr,"[cudss_block] cuDSS %s:%d: status=%d\n",__FILE__,__LINE__,(int)s);exit(1);}}while(0)
#define NCCL_CHECK(call) do { ncclResult_t r=(call); if(r!=ncclSuccess){fprintf(stderr,"[cudss_block] NCCL %s:%d: %s\n",__FILE__,__LINE__,ncclGetErrorString(r));exit(1);}}while(0)

struct CudssBlock {
    int   n_global;       /* plane size (global)            */
    int   local_n;        /* local rows on this rank        */
    int   first_row;      /* plane-relative first row (0-based) */
    long long local_nnz;
    int   plane_size;     /* # ranks owning this plane      */
    int   plane_rank;     /* rank within plane_comm         */
    int   mgmn;           /* 1 if MGMN (plane_size>1)       */

    /* device CSR */
    int    *d_rowptr;
    int    *d_colidx;
    double *d_vals;
    /* device dense rhs/sol (local segment) */
    double *d_rhs;
    double *d_sol;

    cudssHandle_t handle;
    cudssConfig_t config;
    cudssData_t   data;
    cudssMatrix_t mat, rhs_mat, sol_mat;

    ncclComm_t *nccl_comm; /* heap-allocated for MGMN+NCCL; NULL otherwise */

    /* cray-mpich comm-layer backend (CUDSS_COMM_BACKEND=mpi): cuDSS reads the
     * communicator from a STABLE address, so keep the plane MPI_Comm here (a
     * struct member that outlives the solve), NOT a stack temporary. */
    int      use_mpi_comm;       /* 1 if MGMN over the cray-mpich shim */
    MPI_Comm plane_comm_stored;  /* plane communicator passed to cuDSS */

    /* timers */
    double analysis_ms, factor_ms, solve_ms;
    int    n_solves;
    cudaEvent_t ev0, ev1;
};

extern "C"
CudssBlock *cudss_block_create(MPI_Comm plane_comm,
                               int n_global_plane,
                               int local_n,
                               int local_first_row,
                               const int *h_rowptr,
                               const int *h_colidx,
                               const double *h_vals,
                               long long local_nnz,
                               const char *comm_lib_path)
{
    CudssBlock *blk = (CudssBlock *)calloc(1, sizeof(CudssBlock));
    blk->n_global   = n_global_plane;
    blk->local_n    = local_n;
    blk->first_row  = local_first_row;
    blk->local_nnz  = local_nnz;

    MPI_Comm_size(plane_comm, &blk->plane_size);
    MPI_Comm_rank(plane_comm, &blk->plane_rank);
    blk->mgmn = (blk->plane_size > 1) ? 1 : 0;

    /* Comm-layer backend for MGMN: NCCL (default) or the custom cray-mpich shim
     * (CUDSS_COMM_BACKEND=mpi). The shim reads the plane communicator from a
     * stable struct member; store it now (harmless when not used). */
    blk->use_mpi_comm = 0;
    {
        const char *be = getenv("CUDSS_COMM_BACKEND");
        if (be && strcmp(be, "mpi") == 0) blk->use_mpi_comm = 1;
    }
    blk->plane_comm_stored = plane_comm;
    if (blk->mgmn) {
        fprintf(stderr, "[cudss_block] MGMN comm backend: %s (plane_size=%d)\n",
                blk->use_mpi_comm ? "mpi (cray-mpich shim)" : "nccl",
                blk->plane_size);
    }

    /* Use the device PETSc already selected for this rank's vectors. */
    int dev = 0;
    CUDA_CHECK(cudaGetDevice(&dev));

    CUDA_CHECK(cudaEventCreate(&blk->ev0));
    CUDA_CHECK(cudaEventCreate(&blk->ev1));

    /* ---- upload local CSR + allocate rhs/sol ---- */
    CUDA_CHECK(cudaMalloc(&blk->d_rowptr, (size_t)(local_n + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&blk->d_colidx, (size_t)local_nnz * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&blk->d_vals,   (size_t)local_nnz * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&blk->d_rhs,    (size_t)local_n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&blk->d_sol,    (size_t)local_n * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(blk->d_rowptr, h_rowptr, (size_t)(local_n + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(blk->d_colidx, h_colidx, (size_t)local_nnz * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(blk->d_vals,   h_vals,   (size_t)local_nnz * sizeof(double), cudaMemcpyHostToDevice));

    /* ---- NCCL comm for MGMN (skipped for the cray-mpich shim) ---- */
    if (blk->mgmn && !blk->use_mpi_comm) {
        blk->nccl_comm = (ncclComm_t *)malloc(sizeof(ncclComm_t));
        ncclUniqueId id;
        if (blk->plane_rank == 0) NCCL_CHECK(ncclGetUniqueId(&id));
        MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, plane_comm);
        NCCL_CHECK(ncclCommInitRank(blk->nccl_comm, blk->plane_size, id, blk->plane_rank));
    }

    /* ---- cuDSS handle / config / data ----
     * Single-process multi-GPU (MG) mode is activated ONLY by creating the handle
     * with cudssCreateMg(device_count, device_indices).  Setting
     * CUDSS_CONFIG_DEVICE_COUNT/INDICES on a plain cudssCreate() handle is
     * silently ignored -> the factorization runs entirely on one GPU (confirmed
     * by nsys: with the config-only approach GPUs 1..N-1 stayed at 0% util /
     * 4 MiB).  CUDSS_MG_DEVICES=N selects N GPUs {0..N-1} of the (per-rank) set
     * visible via CUDA_VISIBLE_DEVICES.  Not compatible with MGMN (comm layer). */
    int mg_ndev = 0;
    if (!blk->mgmn) {
        if (const char *mg = getenv("CUDSS_MG_DEVICES")) mg_ndev = atoi(mg);
    }
    if (mg_ndev > 1) {
        if (mg_ndev > 16) mg_ndev = 16;
        /* cudssCreateMg requires the process's CURRENT CUDA device to be the
         * primary (device_indices[0]) -- that is also where PETSc placed the
         * input CSR / RHS.  Under MPI, PETSc may bind each rank to a different
         * device (rank%ngpu), so we must put the current device first and fill
         * the remaining slots with the node's other GPUs.  (Passing a fixed
         * {0..N-1} made cudssCreateMg return INVALID_VALUE on every rank whose
         * current device != 0.) */
        int cur_dev = 0, total_dev = 0;
        CUDA_CHECK(cudaGetDevice(&cur_dev));
        CUDA_CHECK(cudaGetDeviceCount(&total_dev));
        static int mg_dev_idx[16];
        mg_dev_idx[0] = cur_dev;
        int k = 1;
        for (int d = 0; d < total_dev && k < mg_ndev; d++)
            if (d != cur_dev) mg_dev_idx[k++] = d;
        mg_ndev = k; /* clamp to devices actually available on this node */
        fprintf(stderr, "[cudss_block] MG cudssCreateMg: cur_dev=%d device_count=%d indices=[", cur_dev, mg_ndev);
        for (int i = 0; i < mg_ndev; i++) fprintf(stderr, "%d%s", mg_dev_idx[i], i+1<mg_ndev?",":"");
        fprintf(stderr, "]\n");
        CUDSS_CHECK(cudssCreateMg(&blk->handle, mg_ndev, mg_dev_idx));
    } else {
        CUDSS_CHECK(cudssCreate(&blk->handle));
    }
    if (blk->mgmn) {
        CUDSS_CHECK(cudssSetCommLayer(blk->handle, comm_lib_path));
    }
    CUDSS_CHECK(cudssConfigCreate(&blk->config));
    CUDSS_CHECK(cudssDataCreate(blk->handle, &blk->data));

    cudssPivotType_t pivot_type = CUDSS_PIVOT_NONE;
    double pivot_epsilon = 1e-12;
    /* Root-cause experiment (2026-07-07): the default PIVOT_NONE relies on cuDSS
     * matching to permute large entries onto the diagonal; with matching OFF
     * (mandatory under MGMN) it hits tiny pivots -> inaccurate LU -> ~2566 iters.
     * SuperLU_DIST stays at 25 iters even with matching+equil off because it does
     * static pivoting. Let CUDSS_PIVOT select a real pivoting strategy so the
     * factor is accurate WITHOUT matching (the MGMN-legal path):
     *   0|none=NONE (default), 1|auto=AUTO, 2|local=LOCAL_BLOCK (general + ND). */
    if (const char *pv = getenv("CUDSS_PIVOT")) {
        if      (!strcmp(pv,"1") || !strcmp(pv,"auto"))  pivot_type = CUDSS_PIVOT_AUTO;
        else if (!strcmp(pv,"2") || !strcmp(pv,"local")) pivot_type = CUDSS_PIVOT_LOCAL_BLOCK;
        else if (!strcmp(pv,"gcol"))                     pivot_type = CUDSS_PIVOT_GLOBAL_COL;
        else if (!strcmp(pv,"grow"))                     pivot_type = CUDSS_PIVOT_GLOBAL_ROW;
        else                                             pivot_type = CUDSS_PIVOT_NONE;
        fprintf(stderr, "[cudss_block] pivot override: CUDSS_PIVOT=%s -> type=%d\n",
                pv, (int)pivot_type);
    }
    /* Reordering override (SuperLU-style column-ordering probe): COLAMD/BTF_COLAMD
     * are required for GLOBAL_COL/GLOBAL_ROW column pivoting.
     *   CUDSS_REORDER: default, btf_colamd, colamd, amd, nd, none */
    int reorder_env = 0;
    cudssReorderingAlg_t reorder_ovr = CUDSS_REORDERING_ALG_DEFAULT;
    if (const char *ro = getenv("CUDSS_REORDER")) {
        reorder_env = 1;
        if      (!strcmp(ro,"btf_colamd")) reorder_ovr = CUDSS_REORDERING_ALG_BTF_COLAMD;
        else if (!strcmp(ro,"colamd"))     reorder_ovr = CUDSS_REORDERING_ALG_COLAMD;
        else if (!strcmp(ro,"amd"))        reorder_ovr = CUDSS_REORDERING_ALG_AMD;
        else if (!strcmp(ro,"nd"))         reorder_ovr = CUDSS_REORDERING_ALG_NESTED_DISSECTION;
        else if (!strcmp(ro,"none"))       reorder_ovr = CUDSS_REORDERING_ALG_NONE;
        else                               reorder_ovr = CUDSS_REORDERING_ALG_DEFAULT;
        fprintf(stderr, "[cudss_block] reorder override: CUDSS_REORDER=%s -> alg=%d\n",
                ro, (int)reorder_ovr);
    }
    /* Iterative refinement (the classic partner to static pivoting; SuperLU does
     * it by default, cuDSS defaults to 0). Each solve then does k correction
     * steps x <- x + M^{-1}(b - A x), which can recover accuracy of an
     * ill-conditioned matching-off factor. MGMN-legal.  CUDSS_IR = n steps. */
    if (const char *ir = getenv("CUDSS_IR")) {
        int nsteps = atoi(ir);
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_IR_N_STEPS, &nsteps, sizeof(nsteps)));
        fprintf(stderr, "[cudss_block] IR_N_STEPS override: %d\n", nsteps);
    }
    /* Pivot threshold: how aggressively active pivoting (LOCAL_BLOCK/GLOBAL_*)
     * swaps to avoid small pivots (higher -> more stable, more fill). Only
     * meaningful when a pivoting strategy is active. CUDSS_PIVOT_THRESHOLD=<dbl>. */
    if (const char *pt = getenv("CUDSS_PIVOT_THRESHOLD")) {
        double thr = atof(pt);
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_PIVOT_THRESHOLD, &thr, sizeof(thr)));
        fprintf(stderr, "[cudss_block] pivot_threshold override: %g\n", thr);
    }
    /* Factorization algorithm variant. CUDSS_FACT_ALG: default/multiblock/general. */
    if (const char *fa = getenv("CUDSS_FACT_ALG")) {
        cudssFactorizationAlg_t a = CUDSS_FACTORIZATION_ALG_DEFAULT;
        if      (!strcmp(fa,"multiblock")) a = CUDSS_FACTORIZATION_ALG_MULTIBLOCK;
        else if (!strcmp(fa,"general"))    a = CUDSS_FACTORIZATION_ALG_GENERAL;
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_FACTORIZATION_ALG, &a, sizeof(a)));
        fprintf(stderr, "[cudss_block] factorization_alg override: %s -> %d\n", fa, (int)a);
    }
    /* Bitwise-reproducible factorization/solve. cuDSS defaults to the faster
     * non-deterministic scheduling, which makes the recomputed true residual
     * ||Ax-b|| differ in its last digits between otherwise identical solves.
     * Costs performance, so it stays opt-in. CUDSS_DETERMINISTIC=1. */
    if (const char *dt = getenv("CUDSS_DETERMINISTIC")) {
        int det = atoi(dt);
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_DETERMINISTIC_MODE, &det, sizeof(det)));
        fprintf(stderr, "[cudss_block] deterministic_mode: %d\n", det);
    }
    CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_PIVOT_TYPE,    &pivot_type,    sizeof(pivot_type)));
    /* Tiny-pivot replacement value + algorithm. cuDSS default is an ABSOLUTE
     * epsilon (1e-12); SuperLU_DIST replaces tiny pivots RELATIVE to the matrix
     * norm, which is why it stays accurate without matching. Expose both:
     *   CUDSS_PIVOT_EPS      = replacement value (double, e.g. 1e-8)
     *   CUDSS_PIVOT_EPS_ALG  = 0 DEFAULT, 1 SCALED (relative), 2 STATIC */
    if (const char *pe = getenv("CUDSS_PIVOT_EPS")) {
        pivot_epsilon = atof(pe);
        fprintf(stderr, "[cudss_block] pivot_epsilon override: %g\n", pivot_epsilon);
    }
    if (const char *pea = getenv("CUDSS_PIVOT_EPS_ALG")) {
        cudssPivotEpsilonAlg_t ea = CUDSS_PIVOT_EPSILON_ALG_DEFAULT;
        if      (!strcmp(pea,"1") || !strcmp(pea,"scaled")) ea = CUDSS_PIVOT_EPSILON_ALG_SCALED;
        else if (!strcmp(pea,"2") || !strcmp(pea,"static")) ea = CUDSS_PIVOT_EPSILON_ALG_STATIC;
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_PIVOT_EPSILON_ALG, &ea, sizeof(ea)));
        fprintf(stderr, "[cudss_block] pivot_epsilon_alg override: %s -> %d\n", pea, (int)ea);
    }
    CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_PIVOT_EPSILON, &pivot_epsilon, sizeof(pivot_epsilon)));
    if (blk->mgmn) {
        /* MGMN (distributed): the distributed solver needs the parallel
         * nested-dissection reordering (was CUDSS_ALG_3 in 0.7.1, now the
         * explicit CUDSS_REORDERING_ALG_NESTED_DISSECTION); ALG_DEFAULT/AMD
         * produced an inconsistent elimination tree across ranks in 0.7.1. */
        cudssReorderingAlg_t reorder = reorder_env ? reorder_ovr
                                                   : CUDSS_REORDERING_ALG_NESTED_DISSECTION;
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_REORDERING_ALG, &reorder, sizeof(reorder)));
        /* HEADLINE EXPERIMENT (cuDSS 0.8.0): in 0.7.1 MGMN forbade matching
         * (it forced ~2561 iters for this matrix). 0.8.0 replaced
         * CUDSS_CONFIG_USE_MATCHING with CUDSS_CONFIG_MATCHING_ALG. Try
         * enabling AUTO matching under MGMN -- gated by CUDSS_MGMN_MATCH
         * (default 1) so we can A/B matching-on vs -off. If ANALYSIS/FACTOR
         * now succeed and iters drop toward ~25, that is the breakthrough. */
        int mgmn_match = 1;
        if (getenv("CUDSS_MGMN_MATCH")) mgmn_match = atoi(getenv("CUDSS_MGMN_MATCH"));
        cudssMatchingAlg_t match = mgmn_match ? CUDSS_MATCHING_ALG_AUTO
                                              : CUDSS_MATCHING_ALG_NONE;
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_MATCHING_ALG, &match, sizeof(match)));
        fprintf(stderr, "[cudss_block] MGMN matching: %s (CUDSS_MGMN_MATCH=%d)\n",
                mgmn_match ? "AUTO" : "NONE", mgmn_match);
        /* Optional host-memory (hybrid) mode for dense layouts where several
         * ranks share a GPU (4 ranks/GPU at 16r/1node and 128r/8node): set
         * CUDSS_BLOCK_HYBRID=1 to spill factors to host and avoid GPU OOM.
         * 0.7.1 note: CUDSS_CONFIG_HYBRID_MODE (now HYBRID_MEMORY_MODE) alone
         * reaches the solve but still OOMs on transient solve buffers at 4
         * ranks/GPU; left as an experimentation knob. */
        if (getenv("CUDSS_BLOCK_HYBRID") && atoi(getenv("CUDSS_BLOCK_HYBRID"))) {
            int hyb = 1;
            CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_HYBRID_MEMORY_MODE, &hyb, sizeof(hyb)));
        }
    } else {
        /* Single-GPU: match the proven hybrid config (matching enabled,
         * default reordering).  CUDSS_NO_MATCH=1 disables matching for the
         * isolation experiment.  0.8.0: USE_MATCHING -> MATCHING_ALG
         * (AUTO when on, NONE when off). */
        int use_matching = 1;
        if (getenv("CUDSS_NO_MATCH") && atoi(getenv("CUDSS_NO_MATCH"))) use_matching = 0;
        cudssMatchingAlg_t match = use_matching ? CUDSS_MATCHING_ALG_AUTO
                                                : CUDSS_MATCHING_ALG_NONE;
        CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_MATCHING_ALG, &match, sizeof(match)));
        if (reorder_env) {
            CUDSS_CHECK(cudssConfigSet(blk->config, CUDSS_CONFIG_REORDERING_ALG, &reorder_ovr, sizeof(reorder_ovr)));
        }
        /* MG single-process multi-GPU is activated at handle creation via
         * cudssCreateMg (see above), NOT via CUDSS_CONFIG_DEVICE_* here. */
    }

    if (blk->mgmn) {
        if (blk->use_mpi_comm) {
            /* Hand cuDSS the address of our stable MPI_Comm member. cuDSS copies
             * `size` bytes into its own slot and later passes &slot as `comm`;
             * the shim reads it back as *(MPI_Comm*)comm. Under MPICH MPI_Comm
             * is a 4-byte int written into the 8-byte slot -- safe (the plan's
             * CommSplit caveat).  0.8.0 split CUDSS_DATA_COMM into DEVICE+HOST;
             * the cray-mpich shim drives both host and device buffers, so set
             * the same plane communicator for both. */
            CUDSS_CHECK(cudssDataSet(blk->handle, blk->data, CUDSS_DATA_COMM_DEVICE,
                                     &blk->plane_comm_stored, sizeof(MPI_Comm *)));
            CUDSS_CHECK(cudssDataSet(blk->handle, blk->data, CUDSS_DATA_COMM_HOST,
                                     &blk->plane_comm_stored, sizeof(MPI_Comm *)));
        } else {
            /* NCCL is device-only; set the device communicator. */
            CUDSS_CHECK(cudssDataSet(blk->handle, blk->data, CUDSS_DATA_COMM_DEVICE,
                                     blk->nccl_comm, sizeof(ncclComm_t *)));
        }
    }

    /* ---- matrix descriptors ----
     * MGMN (Row1d): all ranks pass GLOBAL nrows/ncols/nnz; the device arrays
     * hold only this rank's local rows (local 0-based rowptr, plane-local
     * column indices) and Row1d declares the local global-row range.
     * Single-GPU: everything is local == global. */
    long long global_nnz = local_nnz;
    if (blk->mgmn) {
        MPI_Allreduce(&local_nnz, &global_nnz, 1, MPI_LONG_LONG, MPI_SUM, plane_comm);
    }
    int64_t mat_n   = blk->mgmn ? (int64_t)n_global_plane : (int64_t)local_n;
    int64_t mat_nnz = (int64_t)global_nnz;
    CUDSS_CHECK(cudssMatrixCreateCsr(&blk->mat, mat_n, mat_n, mat_nnz,
        blk->d_rowptr, NULL, blk->d_colidx, blk->d_vals,
        CUDSS_R_32I, CUDSS_R_32I, CUDSS_R_64F, CUDSS_MTYPE_GENERAL, CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO));

    int64_t dn_n = blk->mgmn ? (int64_t)n_global_plane : (int64_t)local_n;
    CUDSS_CHECK(cudssMatrixCreateDn(&blk->rhs_mat, dn_n, 1, (int64_t)local_n,
        blk->d_rhs, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));
    CUDSS_CHECK(cudssMatrixCreateDn(&blk->sol_mat, dn_n, 1, (int64_t)local_n,
        blk->d_sol, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));

    if (blk->mgmn) {
        int64_t first = (int64_t)local_first_row;
        int64_t last  = (int64_t)(local_first_row + local_n - 1); /* INCLUDED */
        CUDSS_CHECK(cudssMatrixSetDistributionRow1d(blk->mat,     first, last));
        CUDSS_CHECK(cudssMatrixSetDistributionRow1d(blk->rhs_mat, first, last));
        CUDSS_CHECK(cudssMatrixSetDistributionRow1d(blk->sol_mat, first, last));
    }

    /* ---- ANALYSIS ---- */
    {
        auto t0 = std::chrono::steady_clock::now();
        CUDSS_CHECK(cudssExecute(blk->handle, CUDSS_PHASE_ANALYSIS, blk->config, blk->data,
                                 blk->mat, blk->sol_mat, blk->rhs_mat));
        CUDA_CHECK(cudaDeviceSynchronize());
        blk->analysis_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
    /* ---- FACTORIZATION ---- */
    {
        auto t0 = std::chrono::steady_clock::now();
        CUDSS_CHECK(cudssExecute(blk->handle, CUDSS_PHASE_FACTORIZATION, blk->config, blk->data,
                                 blk->mat, blk->sol_mat, blk->rhs_mat));
        CUDA_CHECK(cudaDeviceSynchronize());
        blk->factor_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }

    return blk;
}

extern "C"
double cudss_block_refactor(CudssBlock *blk)
{
    /* Re-submit the (unchanged) device values -- an M3D-C1 timestep would upload
     * new values here; the factorization cost is identical for a fixed pattern.
     * The symbolic ANALYSIS from cudss_block_create is reused (we call only the
     * FACTORIZATION phase). */
    CUDSS_CHECK(cudssMatrixSetValues(blk->mat, blk->d_vals));
    auto t0 = std::chrono::steady_clock::now();
    CUDSS_CHECK(cudssExecute(blk->handle, CUDSS_PHASE_FACTORIZATION, blk->config, blk->data,
                             blk->mat, blk->sol_mat, blk->rhs_mat));
    CUDA_CHECK(cudaDeviceSynchronize());
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    blk->factor_ms = ms;
    return ms;
}

extern "C"
void cudss_block_solve(CudssBlock *blk, const double *d_rhs, double *d_sol)
{
    /* Stage the PETSc-provided local segment into our owned buffers so the
     * cuDSS dense descriptors keep stable device pointers across applies. */
    CUDA_CHECK(cudaMemcpy(blk->d_rhs, d_rhs, (size_t)blk->local_n * sizeof(double), cudaMemcpyDeviceToDevice));

    CUDA_CHECK(cudaEventRecord(blk->ev0));
    CUDSS_CHECK(cudssExecute(blk->handle, CUDSS_PHASE_SOLVE, blk->config, blk->data,
                             blk->mat, blk->sol_mat, blk->rhs_mat));
    CUDA_CHECK(cudaEventRecord(blk->ev1));
    CUDA_CHECK(cudaEventSynchronize(blk->ev1));
    float ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, blk->ev0, blk->ev1));
    blk->solve_ms += ms;
    blk->n_solves++;

    CUDA_CHECK(cudaMemcpy(d_sol, blk->d_sol, (size_t)blk->local_n * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
}

extern "C"
void cudss_block_get_timers(CudssBlock *blk, double *analysis_ms, double *factor_ms,
                            double *solve_ms, int *n_solves)
{
    if (analysis_ms) *analysis_ms = blk->analysis_ms;
    if (factor_ms)   *factor_ms   = blk->factor_ms;
    if (solve_ms)    *solve_ms    = blk->solve_ms;
    if (n_solves)    *n_solves    = blk->n_solves;
}

extern "C"
void cudss_block_reset_solve_timer(CudssBlock *blk)
{
    blk->solve_ms = 0.0;
    blk->n_solves = 0;
}

extern "C"
void cudss_block_destroy(CudssBlock *blk)
{
    if (!blk) return;
    cudssMatrixDestroy(blk->mat);
    cudssMatrixDestroy(blk->rhs_mat);
    cudssMatrixDestroy(blk->sol_mat);
    cudssDataDestroy(blk->handle, blk->data);
    cudssConfigDestroy(blk->config);
    cudssDestroy(blk->handle);
    if (blk->nccl_comm) { ncclCommDestroy(*blk->nccl_comm); free(blk->nccl_comm); }
    cudaFree(blk->d_rowptr);
    cudaFree(blk->d_colidx);
    cudaFree(blk->d_vals);
    cudaFree(blk->d_rhs);
    cudaFree(blk->d_sol);
    cudaEventDestroy(blk->ev0);
    cudaEventDestroy(blk->ev1);
    free(blk);
}
