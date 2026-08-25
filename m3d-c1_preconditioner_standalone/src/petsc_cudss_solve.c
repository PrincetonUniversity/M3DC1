// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * PETSc standalone solver for the M3D-C1 extracted matrix (matrix_id=5),
 * using cuDSS as the per-plane block-Jacobi sub-solver via a PCShell.
 *
 * Derived from petsc_standalone_solve.c, the SuperLU_DIST baseline driver,
 * which is not part of this package.  The Krylov outer solve (LGMRES,
 * aijcusparse/cuda) is identical; the preconditioner is replaced by a custom
 * PCShell whose Apply runs a cuDSS LU solve on this rank's piece of its plane
 * diagonal block.  When a plane is owned by a single rank (Config B) the cuDSS
 * solve is single-GPU; when a plane spans >1 rank (Config A) cuDSS runs in MGMN
 * mode over the plane's sub-communicator, using either the NCCL comm layer
 * (default) or the cray-mpich shim in cudss_commlayer_craympi.c
 * (CUDSS_COMM_BACKEND=mpi).  See cudss_block.{h,cu}.
 *
 * Requires plane-aligned distribution: nplanes>0, nprocs%nplanes==0.
 *
 * Mat/Vec types are controlled by -mat_type / -vec_type (use aijcusparse/cuda).
 *
 * Build: split nvcc (.cu) + Cray cc (.c) + link; see build.sbatch.
 *
 * Usage (Config B, the headline configuration -- 4 planes, 1 rank/GPU each):
 *   srun -n 4 ./petsc_cudss_solve \
 *     -A A.bin -b b.bin [-x0 x0.bin] [-xref xout.bin] -nplanes 4 \
 *     -mat_type aijcusparse -vec_type cuda \
 *     -ksp_type lgmres ... -log_view [...]
 *
 * Options owned by this driver (all others are PETSc's):
 *   -A <f> -b <f>       matrix / RHS in PETSc binary format (both required)
 *   -x0 <f> -xref <f>   optional initial guess / reference solution
 *   -nplanes <N>        plane-aligned row distribution (required, >0).  Must be the
 *                       toroidal plane count OF THE MATRIX (4 for the shipped A.bin),
 *                       not a tuning knob: the Jacobi blocks are the plane diagonal
 *                       blocks, so a divisible-but-wrong value preconditions with
 *                       something else, converges, and reports meaningless numbers.
 *   -nreps <N>          repeat the solve N times (default 1).  With N>1 the
 *                       last rep is logged in its own "Benchmark" stage and the
 *                       earlier ones in "Warmup".  Each rep is self-contained
 *                       (extract + ANALYSIS + FACTORIZATION + solve), matching
 *                       the library entry point petsc_cudss_solve() that
 *                       m3dc1_scorec calls per timestep.
 *   -cudss_comm_lib <f> cuDSS comm-layer .so; MGMN (Config A) only.  Defaults
 *                       to ${CUDSS_DIR:-deps/cudss}/lib/libcudss_commlayer_nccl.so.0,
 *                       i.e. the vendored NCCL layer.  For the cray-mpich shim
 *                       instead: CUDSS_COMM_BACKEND=mpi with
 *                       -cudss_comm_lib bin/libcudss_commlayer_craympi.so.
 *
 * -mat_cusparse_spmv_alg is a PETSc option, but this driver additionally
 * forwards it to the MPIAIJCUSPARSE sub-blocks so that the reproducible
 * CSRMV_ALG2 is reachable on more than one rank (see ReadMatrixParallel).
 * cuDSS itself is tuned through CUDSS_* environment variables
 * (CUDSS_DETERMINISTIC, CUDSS_REORDER, CUDSS_PIVOT, ...) -- see cudss_block.cu.
 */
#include <petsc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include "cudss_block.h"
#include "petsc_cudss_solve.h"

static void swap4(void *p) { char *c=(char*)p, t=c[0]; c[0]=c[3]; c[3]=t; t=c[1]; c[1]=c[2]; c[2]=t; }
static void swap8(void *p) { char *c=(char*)p; for(int i=0;i<4;i++){char t=c[i];c[i]=c[7-i];c[7-i]=t;} }
static int read_i32(FILE *f) { int v; fread(&v,4,1,f); swap4(&v); return v; }

static void compute_row_range(PetscInt n_global, PetscMPIInt rank, PetscMPIInt nprocs,
                              PetscInt nplanes, PetscInt *out_start, PetscInt *out_n)
{
    if (nplanes > 0 && nprocs >= nplanes && (nprocs % nplanes) == 0) {
        PetscInt rows_per_plane = n_global / nplanes;
        PetscInt ranks_per_plane = nprocs / nplanes;
        PetscInt plane = rank / ranks_per_plane;
        PetscInt local_rank = rank % ranks_per_plane;
        PetscInt plane_start = plane * rows_per_plane;
        PetscInt base = rows_per_plane / ranks_per_plane;
        PetscInt rem = rows_per_plane % ranks_per_plane;
        *out_start = plane_start + base * local_rank + ((local_rank < rem) ? local_rank : rem);
        *out_n = base + ((local_rank < rem) ? 1 : 0);
    } else {
        PetscInt base = n_global / nprocs;
        PetscInt rem = n_global % nprocs;
        *out_start = base * rank + ((rank < rem) ? rank : rem);
        *out_n = base + ((rank < rem) ? 1 : 0);
    }
}

static PetscErrorCode ReadMatrixParallel(const char *path, Mat *pA, MPI_Comm comm, PetscInt nplanes)
{
    PetscFunctionBeginUser;
    PetscMPIInt rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    PetscInt n_global, nnz_global;
    PetscInt *row_lengths = NULL;

    FILE *f = fopen(path, "rb");
    if (!f) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open matrix file");

    int classid = read_i32(f);
    if (classid != 1211216) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Bad classid for matrix");
    n_global = read_i32(f);
    read_i32(f);
    nnz_global = read_i32(f);

    PetscPrintf(comm, "  %d x %d, nnz = %d\n", (int)n_global, (int)n_global, (int)nnz_global);

    PetscCall(PetscMalloc1(n_global, &row_lengths));
    for (PetscInt i = 0; i < n_global; i++) {
        int v; fread(&v, 4, 1, f); swap4(&v);
        row_lengths[i] = v;
    }

    PetscInt local_start, local_n;
    compute_row_range(n_global, rank, nprocs, nplanes, &local_start, &local_n);
    PetscInt local_end = local_start + local_n;

    long long nnz_before = 0;
    for (PetscInt i = 0; i < local_start; i++) nnz_before += row_lengths[i];
    long long local_nnz = 0;
    for (PetscInt i = local_start; i < local_end; i++) local_nnz += row_lengths[i];

    if (rank == 0)
        PetscPrintf(comm, "  Rank 0: rows [0, %d), local_nnz = %lld\n", (int)local_n, local_nnz);

    long long col_base = 16LL + (long long)n_global * 4;
    long long val_base = col_base + (long long)nnz_global * 4;

    PetscInt *col_idx;
    PetscCall(PetscMalloc1(local_nnz, &col_idx));
    fseek(f, col_base + nnz_before * 4, SEEK_SET);
    {
        long long remaining = local_nnz, offset = 0;
        while (remaining > 0) {
            long long chunk = (remaining > 100000000LL) ? 100000000LL : remaining;
            if ((long long)fread(col_idx + offset, 4, (size_t)chunk, f) != chunk)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "col read fail");
            for (long long i = 0; i < chunk; i++) swap4(&col_idx[offset + i]);
            offset += chunk; remaining -= chunk;
        }
    }

    PetscScalar *val;
    PetscCall(PetscMalloc1(local_nnz, &val));
    fseek(f, val_base + nnz_before * 8, SEEK_SET);
    {
        long long remaining = local_nnz, offset = 0;
        while (remaining > 0) {
            long long chunk = (remaining > 100000000LL) ? 100000000LL : remaining;
            if ((long long)fread(val + offset, 8, (size_t)chunk, f) != chunk)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "val read fail");
            for (long long i = 0; i < chunk; i++) swap8(&val[offset + i]);
            offset += chunk; remaining -= chunk;
        }
    }
    fclose(f);

    PetscInt *d_nnz, *o_nnz;
    PetscCall(PetscMalloc2(local_n, &d_nnz, local_n, &o_nnz));
    {
        long long pos = 0;
        for (PetscInt i = 0; i < local_n; i++) {
            PetscInt global_row = local_start + i;
            PetscInt rl = row_lengths[global_row];
            PetscInt d = 0, o = 0;
            for (PetscInt j = 0; j < rl; j++) {
                if (col_idx[pos + j] >= local_start && col_idx[pos + j] < local_end) d++;
                else o++;
            }
            d_nnz[i] = d;
            o_nnz[i] = o;
            pos += rl;
        }
    }

    PetscCall(MatCreate(comm, pA));
    PetscCall(MatSetSizes(*pA, local_n, local_n, n_global, n_global));
    PetscCall(MatSetFromOptions(*pA));
    PetscCall(MatSeqAIJSetPreallocation(*pA, 0, d_nnz));
    PetscCall(MatMPIAIJSetPreallocation(*pA, 0, d_nnz, 0, o_nnz));
    PetscCall(PetscFree2(d_nnz, o_nnz));

    {
        long long pos = 0;
        for (PetscInt i = 0; i < local_n; i++) {
            PetscInt global_row = local_start + i;
            PetscInt rl = row_lengths[global_row];
            PetscCall(MatSetValues(*pA, 1, &global_row, rl, &col_idx[pos], &val[pos], INSERT_VALUES));
            pos += rl;
        }
    }

    PetscCall(MatAssemblyBegin(*pA, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*pA, MAT_FINAL_ASSEMBLY));

    /* PETSc reads -mat_cusparse_spmv_alg only in MatSetFromOptions_SeqAIJCUSPARSE,
     * so on more than one rank the MPIAIJCUSPARSE sub-blocks that actually run the
     * SpMV never see it and CSRMV_ALG2 (the reproducible one) is unreachable.
     * Re-issue it to the sub-blocks under a private prefix: calling a bare
     * MatSetFromOptions on them would re-read -mat_type and destroy the values we
     * just assembled. Has to happen before the first MatMult sizes the cuSPARSE
     * buffer for the current algorithm. */
    {
        char        alg[64];
        PetscBool   flg;
        PetscMPIInt size;

        PetscCall(PetscOptionsGetString(NULL, NULL, "-mat_cusparse_spmv_alg", alg, sizeof(alg), &flg));
        PetscCallMPI(MPI_Comm_size(comm, &size));
        if (flg && size > 1) {
            Mat sub[2];
            PetscCall(MatMPIAIJGetSeqAIJ(*pA, &sub[0], &sub[1], NULL));
            PetscCall(PetscOptionsSetValue(NULL, "-m3dc1_sub_mat_cusparse_spmv_alg", alg));
            for (int s = 0; s < 2; s++) {
                if (!sub[s]) continue;
                PetscCall(MatSetOptionsPrefix(sub[s], "m3dc1_sub_"));
                PetscCall(MatSetFromOptions(sub[s]));
                PetscCall(MatSetOptionsPrefix(sub[s], NULL));
            }
        }
    }

    PetscCall(PetscFree(col_idx));
    PetscCall(PetscFree(val));
    PetscCall(PetscFree(row_lengths));

    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode ReadVectorParallel(const char *path, Vec *pv, MPI_Comm comm, PetscInt nplanes)
{
    PetscFunctionBeginUser;
    PetscMPIInt rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    FILE *f = fopen(path, "rb");
    if (!f) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open vector file");

    read_i32(f);
    PetscInt vn = read_i32(f);

    PetscInt local_start, local_n;
    compute_row_range(vn, rank, nprocs, nplanes, &local_start, &local_n);

    PetscScalar *data;
    PetscCall(PetscMalloc1(local_n, &data));
    fseek(f, 8 + local_start * 8, SEEK_SET);
    if ((PetscInt)fread(data, 8, local_n, f) != local_n)
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Vector read fail");
    for (PetscInt i = 0; i < local_n; i++) swap8(&data[i]);
    fclose(f);

    PetscCall(VecCreate(comm, pv));
    PetscCall(VecSetSizes(*pv, local_n, vn));
    PetscCall(VecSetFromOptions(*pv));
    PetscInt *idx;
    PetscCall(PetscMalloc1(local_n, &idx));
    for (PetscInt i = 0; i < local_n; i++) idx[i] = local_start + i;
    PetscCall(VecSetValues(*pv, local_n, idx, data, INSERT_VALUES));
    PetscCall(PetscFree(idx));
    PetscCall(PetscFree(data));
    PetscCall(VecAssemblyBegin(*pv));
    PetscCall(VecAssemblyEnd(*pv));

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------------------
 * cuDSS block-Jacobi PCShell
 * ------------------------------------------------------------------------- */
typedef struct {
    CudssBlock *blk;     /* persistent cuDSS solver for this rank's plane block */
    MPI_Comm    plane_comm;
} CudssPCCtx;

/* Extract this rank's local rows of the plane DIAGONAL block from the
 * assembled (distributed) matrix A: keep entries whose GLOBAL column lies in
 * [plane_start, plane_end); rebase columns to plane-local (col-plane_start).
 * Off-plane couplings are dropped (the block-Jacobi approximation).  Returns
 * a host int32 CSR (offsets start at 0, local rows only). */
static PetscErrorCode ExtractPlaneBlock(Mat A,
                                        PetscInt plane_start, PetscInt plane_end,
                                        PetscInt local_start, PetscInt local_n,
                                        int **out_rowptr, int **out_colidx,
                                        double **out_vals, long long *out_nnz)
{
    PetscFunctionBeginUser;
    int *rowptr;
    PetscCall(PetscMalloc1(local_n + 1, &rowptr));
    rowptr[0] = 0;
    for (PetscInt i = 0; i < local_n; i++) {
        PetscInt gr = local_start + i, ncols;
        const PetscInt *cols;
        const PetscScalar *vals;
        PetscCall(MatGetRow(A, gr, &ncols, &cols, &vals));
        int cnt = 0;
        for (PetscInt j = 0; j < ncols; j++)
            if (cols[j] >= plane_start && cols[j] < plane_end) cnt++;
        rowptr[i + 1] = rowptr[i] + cnt;
        PetscCall(MatRestoreRow(A, gr, &ncols, &cols, &vals));
    }
    long long nnz = rowptr[local_n];
    int *colidx;
    double *vout;
    PetscCall(PetscMalloc1(nnz, &colidx));
    PetscCall(PetscMalloc1(nnz, &vout));
    long long pos = 0;
    for (PetscInt i = 0; i < local_n; i++) {
        PetscInt gr = local_start + i, ncols;
        const PetscInt *cols;
        const PetscScalar *vals;
        PetscCall(MatGetRow(A, gr, &ncols, &cols, &vals));
        for (PetscInt j = 0; j < ncols; j++) {
            if (cols[j] >= plane_start && cols[j] < plane_end) {
                colidx[pos] = (int)(cols[j] - plane_start);
                vout[pos]   = (double)PetscRealPart(vals[j]);
                pos++;
            }
        }
        PetscCall(MatRestoreRow(A, gr, &ncols, &cols, &vals));
    }
    *out_rowptr = rowptr;
    *out_colidx = colidx;
    *out_vals   = vout;
    *out_nnz    = nnz;
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CudssPCSetUp(PC pc)
{
    PetscFunctionBeginUser;   /* block built in setKspType_cudss(); nothing to do */
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CudssPCApply(PC pc, Vec xin, Vec xout)
{
    PetscFunctionBeginUser;
    CudssPCCtx *ctx;
    PetscCall(PCShellGetContext(pc, &ctx));
    const PetscScalar *d_in;
    PetscScalar *d_out;
    PetscCall(VecCUDAGetArrayRead(xin, &d_in));
    PetscCall(VecCUDAGetArrayWrite(xout, &d_out));
    cudss_block_solve(ctx->blk, (const double *)d_in, (double *)d_out);
    PetscCall(VecCUDARestoreArrayRead(xin, &d_in));
    PetscCall(VecCUDARestoreArrayWrite(xout, &d_out));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Owns the context (allocated in setKspType_cudss): destroys the cuDSS block
 * and the plane communicator when the KSP/PC is destroyed. */
static PetscErrorCode CudssPCDestroy(PC pc)
{
    PetscFunctionBeginUser;
    CudssPCCtx *ctx;
    PetscCall(PCShellGetContext(pc, &ctx));
    if (ctx) {
        cudss_block_destroy(ctx->blk);
        MPI_Comm_free(&ctx->plane_comm);
        PetscCall(PetscFree(ctx));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------------------
 * setKspType_cudss -- one-time KSP + cuDSS preconditioner setup
 * (see petsc_cudss_solve.h)
 *
 * The cuDSS analogue of matrix_solve::setKspType() in m3dc1_matrix.cc:
 * extracts the plane diagonal block from A, runs cuDSS
 * ANALYSIS+FACTORIZATION, creates the KSP (operators + tolerances +
 * KSPSetFromOptions) and installs the cuDSS block-Jacobi PCShell on it.
 * The returned KSP owns the cuDSS block and the plane communicator
 * (freed by CudssPCDestroy inside KSPDestroy).
 *
 * Called from m3dc1_scorec (matrix_solve::solve_cudss /
 * solve_cudss_with_guess) when _kspSet == 0; the caller stores the KSP in
 * _ksp, reuses it across solves, and destroys it in ~matrix_solve.  The
 * standalone main() below calls it per rep, so each rep re-runs the full
 * setup as before.
 * ------------------------------------------------------------------------- */
PetscErrorCode setKspType_cudss(Mat A, PetscInt nplanes, KSP *out_ksp)
{
    PetscFunctionBeginUser;
    MPI_Comm comm;
    PetscCall(PetscObjectGetComm((PetscObject)A, &comm));
    PetscMPIInt rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (nplanes <= 0 || nprocs < nplanes || (nprocs % nplanes) != 0) {
        SETERRQ(comm, PETSC_ERR_ARG_WRONG,
                "setKspType_cudss requires plane-aligned layout: nplanes>0 and nprocs%%nplanes==0 (got nprocs=%d, nplanes=%d)",
                (int)nprocs, (int)nplanes);
    }
    PetscInt n_global;
    PetscCall(MatGetSize(A, &n_global, NULL));
    /* nplanes must be the toroidal plane count OF THIS MATRIX, not a free knob:
     * the block-Jacobi blocks are the plane diagonal blocks.  A divisible but
     * wrong value cannot be detected in general -- it just preconditions with
     * something else and still converges -- so catch at least the indivisible
     * case rather than silently truncating rows. */
    PetscCheck(n_global % nplanes == 0, comm, PETSC_ERR_ARG_INCOMP,
               "nplanes %d does not divide the matrix dimension %d; nplanes must equal the "
               "number of toroidal planes in the matrix",
               (int)nplanes, (int)n_global);
    PetscInt rows_per_plane  = n_global / nplanes;
    PetscInt ranks_per_plane = nprocs / nplanes;
    PetscInt plane           = rank / ranks_per_plane;
    PetscInt plane_start     = plane * rows_per_plane;
    PetscInt plane_end       = plane_start + rows_per_plane;

    /* The matrix owns its row layout (PETSc- or m3dc1-decided); take it from
     * the Mat rather than recomputing, and require plane alignment. */
    PetscInt local_start, local_end, local_n;
    PetscCall(MatGetOwnershipRange(A, &local_start, &local_end));
    local_n = local_end - local_start;
    PetscCheck(local_start >= plane_start && local_end <= plane_end,
               PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
               "rank %d owns rows [%d, %d) which are not contained in its plane %d rows [%d, %d); "
               "the matrix row distribution must be plane-aligned",
               (int)rank, (int)local_start, (int)local_end, (int)plane,
               (int)plane_start, (int)plane_end);

    /* Path to the cuDSS comm-layer .so (used only for MGMN, i.e. when a plane
     * spans >1 rank; one rank per plane never touches it).  The default names
     * the NCCL layer shipped inside cuDSS, following the same
     * ${CUDSS_DIR:-deps/cudss} convention as build.sbatch and the run scripts.
     * On Perlmutter that layer CANNOT load: it is linked against OpenMPI
     * (libmpi.so.40), which does not exist here, and cudssSetCommLayer then
     * fails with CUDSS_STATUS_INVALID_VALUE.  Use CUDSS_COMM_BACKEND=mpi with
     * -cudss_comm_lib <abs>/bin/libcudss_commlayer_craympi.so, the cray-mpich
     * shim built from cudss_commlayer_craympi.c (see run_mgmn.sbatch). */
    char comm_lib[1024];
    PetscBool has_comm_lib;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-cudss_comm_lib", comm_lib, sizeof(comm_lib), &has_comm_lib));
    if (!has_comm_lib) {
        const char *cudss_dir = getenv("CUDSS_DIR");
        snprintf(comm_lib, sizeof(comm_lib), "%s/lib/libcudss_commlayer_nccl.so.0",
                 (cudss_dir && cudss_dir[0]) ? cudss_dir : "deps/cudss");
    }

    /* cuDSS dlopen()s the comm layer and reports every problem with it as an
     * opaque CUDSS_STATUS_INVALID_VALUE from cudssSetCommLayer -- a relative
     * path included -- so resolve and check it here where we can say why. */
    if (ranks_per_plane > 1) {
        char resolved[PATH_MAX];
        PetscCheck(realpath(comm_lib, resolved), PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "-cudss_comm_lib '%s': %s.  MGMN needs a comm layer, cuDSS wants an "
                   "absolute path, and the NCCL layer inside cuDSS is linked against "
                   "OpenMPI (libmpi.so.40) so it cannot load on Perlmutter -- pass "
                   "<pkg>/bin/libcudss_commlayer_craympi.so (see run_mgmn.sbatch).",
                   comm_lib, strerror(errno));
        PetscCall(PetscStrncpy(comm_lib, resolved, sizeof(comm_lib)));
    }

    MPI_Comm plane_comm;
    MPI_Comm_split(comm, (int)plane, (int)(rank % ranks_per_plane), &plane_comm);

    PetscPrintf(comm,
        "cuDSS block-Jacobi PC: %d ranks/plane -> %s\n",
        (int)ranks_per_plane, (ranks_per_plane > 1) ? "MGMN (distributed cuDSS)" : "single-GPU/plane");
    if (ranks_per_plane > 1)
        PetscPrintf(comm, "  comm layer: %s\n", comm_lib);

    int *blk_rowptr, *blk_colidx;
    double *blk_vals;
    long long blk_nnz;
    PetscLogDouble t_ext0, t_ext1;
    PetscCall(PetscTime(&t_ext0));
    PetscCall(ExtractPlaneBlock(A, plane_start, plane_end, local_start, local_n,
                                &blk_rowptr, &blk_colidx, &blk_vals, &blk_nnz));
    PetscCall(PetscTime(&t_ext1));
    PetscPrintf(comm, "  plane block extracted (rank0 local_n=%d, local_nnz=%lld) in %.2f s\n",
                (int)local_n, blk_nnz, t_ext1 - t_ext0);

    CudssPCCtx *pcctx;
    PetscCall(PetscNew(&pcctx));
    PetscLogDouble t_setup0, t_setup1;
    PetscCall(PetscTime(&t_setup0));
    pcctx->blk = cudss_block_create(plane_comm, (int)rows_per_plane, (int)local_n,
                                    (int)(local_start - plane_start),
                                    blk_rowptr, blk_colidx, blk_vals, blk_nnz, comm_lib);
    pcctx->plane_comm = plane_comm;
    PetscCall(PetscTime(&t_setup1));
    PetscCall(PetscFree(blk_rowptr));
    PetscCall(PetscFree(blk_colidx));
    PetscCall(PetscFree(blk_vals));

    {
        double an, fa, so; int ns;
        cudss_block_get_timers(pcctx->blk, &an, &fa, &so, &ns);
        PetscPrintf(comm,
            "  cuDSS setup: ANALYSIS %.1f ms, FACTORIZATION %.1f ms (wall %.2f s)\n",
            an, fa, t_setup1 - t_setup0);
    }

    /* Create the KSP: same pattern as matrix_solve::setKspType(). */
    KSP ksp;
    PetscCall(KSPCreate(comm, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetTolerances(ksp, .000001, .000000001, PETSC_DEFAULT, 1000));
    PetscCall(KSPSetFromOptions(ksp));

    /* Install our cuDSS block-Jacobi as a PCShell (overrides -pc_type).
     * The shell owns pcctx (freed by CudssPCDestroy when the KSP dies). */
    {
        PC pc;
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCSHELL));
        PetscCall(PCShellSetContext(pc, pcctx));
        PetscCall(PCShellSetSetUp(pc, CudssPCSetUp));
        PetscCall(PCShellSetApply(pc, CudssPCApply));
        PetscCall(PCShellSetDestroy(pc, CudssPCDestroy));
        PetscCall(PCShellSetName(pc, "cudss_bjacobi"));
    }

    *out_ksp = ksp;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------------------
 * petsc_cudss_solve -- solve entry point (see petsc_cudss_solve.h)
 *
 * Solves A x = b with the caller-provided KSP, which must have been set up
 * by setKspType_cudss() (cuDSS PCShell installed, factorization done).  The
 * KSP and its cuDSS state are owned by the caller (persistent _ksp in
 * m3dc1_scorec's matrix_solve::solve_cudss; per-rep in main() below) and are
 * NOT destroyed here.
 * ------------------------------------------------------------------------- */
PetscErrorCode petsc_cudss_solve(KSP ksp, Mat A, Vec b, Vec x,
                                 PetscInt nplanes,
                                 PetscBool use_initial_guess,
                                 PetscInt *out_its)
{
    PetscFunctionBeginUser;
    (void)nplanes;   /* plane setup now lives in setKspType_cudss() */
    MPI_Comm comm;
    PetscCall(PetscObjectGetComm((PetscObject)A, &comm));

    if (!use_initial_guess)
        PetscCall(VecZeroEntries(x));
    PetscCall(KSPSetInitialGuessNonzero(ksp, use_initial_guess));

    PetscLogDouble t0, t1;
    PetscCall(PetscTime(&t0));
    PetscCall(KSPSolve(ksp, b, x));
    PetscCall(PetscTime(&t1));

    PetscInt its;
    KSPConvergedReason reason;
    PetscReal rnorm;
    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));
    if (out_its) *out_its = its;

    PetscPrintf(comm, "  cuDSS solve: reason %d, iterations %d, KSP ||r|| %.6e, %.3f s\n",
                reason, (int)its, rnorm, t1 - t0);

    {   /* cuDSS preconditioner-apply timing (max across ranks; cumulative
         * over the KSP's lifetime, since the block persists in the PCShell) */
        PC pc;
        CudssPCCtx *pcctx;
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCShellGetContext(pc, &pcctx));
        double an, fa, so; int ns;
        cudss_block_get_timers(pcctx->blk, &an, &fa, &so, &ns);
        double an_max, fa_max, so_max;
        int ns_max;
        MPI_Reduce(&an, &an_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&fa, &fa_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&so, &so_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&ns, &ns_max, 1, MPI_INT,    MPI_MAX, 0, comm);
        PetscPrintf(comm, "  cuDSS timing (max/rank): ANALYSIS %.1f ms, FACTOR %.1f ms, "
                    "SOLVE %.1f ms (%d applies, %.2f ms/apply)\n",
                    an_max, fa_max, so_max, ns_max, ns_max ? so_max/ns_max : 0.0);
    }

    /* KSP (and the cuDSS block it owns) is caller-owned; nothing torn down. */
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------------------
 * File-based benchmark driver.  Compiled out (-DCUDSS_NO_MAIN) when this file
 * is built into the m3dc1_scorec library.
 * ------------------------------------------------------------------------- */
#ifndef CUDSS_NO_MAIN
int main(int argc, char **argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL,
        "PETSc + cuDSS block-Jacobi solver for the M3D-C1 extracted matrix\n"
        "Usage: srun -n N petsc_cudss_solve -A mat.bin -b rhs.bin -nplanes P [options]\n"));

    /* No PetscLogSetThreshold call here on purpose: it would SEGV on PETSc
     * 3.20+ for plain -log_view runs, where no nested handler is installed.
     * Thresholding is only reachable via -log_threshold, which PETSc parses
     * inside PetscInitialize and honors solely for nested log handlers. */

    char mat_path[512], rhs_path[512], x0_path[512], xref_path[512];
    PetscBool has_mat, has_rhs, has_x0, has_xref;
    PetscInt nplanes = 0;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-A", mat_path, sizeof(mat_path), &has_mat));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-b", rhs_path, sizeof(rhs_path), &has_rhs));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-x0", x0_path, sizeof(x0_path), &has_x0));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-xref", xref_path, sizeof(xref_path), &has_xref));
    PetscInt nreps = 1;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-nplanes", &nplanes, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-nreps", &nreps, NULL));

    /* Fail with a nonzero exit status, so a wrapper script cannot mistake a
     * missing argument for a successful run. */
    PetscCheck(has_mat && has_rhs, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Required: -A <mat.bin> -b <rhs.bin> -nplanes <planes in the matrix>");

    PetscMPIInt rank, nprocs;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    PetscPrintf(PETSC_COMM_WORLD, "=== PETSc Standalone Solver ===\n");
    PetscPrintf(PETSC_COMM_WORLD, "MPI ranks: %d\n", nprocs);
    if (nplanes > 0)
        PetscPrintf(PETSC_COMM_WORLD, "Plane-aligned distribution: %d planes, %d ranks/plane\n",
                    (int)nplanes, nprocs / (int)nplanes);
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    PetscLogDouble t_start, t_load;
    PetscCall(PetscTime(&t_start));

    PetscPrintf(PETSC_COMM_WORLD, "Reading matrix: %s\n", mat_path);
    Mat A;
    PetscCall(ReadMatrixParallel(mat_path, &A, PETSC_COMM_WORLD, nplanes));

    PetscPrintf(PETSC_COMM_WORLD, "Loading RHS: %s\n", rhs_path);
    Vec b;
    PetscCall(ReadVectorParallel(rhs_path, &b, PETSC_COMM_WORLD, nplanes));

    Vec x;
    PetscCall(VecDuplicate(b, &x));

    PetscCall(PetscTime(&t_load));
    PetscPrintf(PETSC_COMM_WORLD, "  Data loaded in %.1f s\n", t_load - t_start);

    PetscReal bnorm;
    PetscCall(VecNorm(b, NORM_2, &bnorm));
    PetscPrintf(PETSC_COMM_WORLD, "\n||b|| = %.6e\n\n", bnorm);

    if (nplanes > 0) {
        PetscInt n_global;
        PetscCall(MatGetSize(A, &n_global, NULL));
        PetscPrintf(PETSC_COMM_WORLD, "Plane size: %d rows/plane\n\n", (int)(n_global / nplanes));
    }

    PetscLogStage stage_warmup, stage_bench;
    if (nreps > 1) {
        PetscCall(PetscLogStageRegister("Warmup", &stage_warmup));
        PetscCall(PetscLogStageRegister("Benchmark", &stage_bench));
    }

    PetscPrintf(PETSC_COMM_WORLD, "Repetitions: %d%s\n\n", (int)nreps,
                nreps > 1 ? " (first is warmup)" : "");

    for (PetscInt rep = 0; rep < nreps; rep++) {
        if (nreps > 1)
            PetscCall(PetscLogStagePush((rep == nreps - 1) ? stage_bench : stage_warmup));

        PetscPrintf(PETSC_COMM_WORLD, "--- Rep %d/%d %s---\n", (int)rep+1, (int)nreps,
                    (nreps > 1 && rep < nreps-1) ? "(warmup) " : "");

        if (has_x0) {
            Vec x0;
            PetscCall(ReadVectorParallel(x0_path, &x0, PETSC_COMM_WORLD, nplanes));
            PetscCall(VecCopy(x0, x));
            PetscCall(VecDestroy(&x0));
        }

        /* Per-rep KSP + cuDSS setup; KSPDestroy tears the cuDSS block down,
         * so each rep re-runs extract + ANALYSIS + FACTORIZATION as before. */
        PetscInt its = 0;
        PetscLogDouble t0, t1;
        PetscCall(PetscTime(&t0));
        KSP ksp;
        PetscCall(setKspType_cudss(A, nplanes, &ksp));
        PetscCall(petsc_cudss_solve(ksp, A, b, x, nplanes, has_x0, &its));
        PetscCall(PetscTime(&t1));
        PetscCall(KSPDestroy(&ksp));

        Vec r;
        PetscCall(VecDuplicate(b, &r));
        PetscCall(MatMult(A, x, r));
        PetscCall(VecAXPY(r, -1.0, b));
        PetscReal true_rnorm;
        PetscCall(VecNorm(r, NORM_2, &true_rnorm));
        PetscCall(VecDestroy(&r));

        PetscPrintf(PETSC_COMM_WORLD, "\n========================================\n");
        PetscPrintf(PETSC_COMM_WORLD, "    RESULTS  (rep %d/%d)\n", (int)rep+1, (int)nreps);
        PetscPrintf(PETSC_COMM_WORLD, "========================================\n");
        PetscPrintf(PETSC_COMM_WORLD, "Iterations:   %d\n", (int)its);
        PetscPrintf(PETSC_COMM_WORLD, "Load time:    %.3f s\n", t_load - t_start);
        PetscPrintf(PETSC_COMM_WORLD, "Setup+solve:  %.3f s\n", t1 - t0);
        PetscPrintf(PETSC_COMM_WORLD, "||b||         = %.6e\n", bnorm);
        PetscPrintf(PETSC_COMM_WORLD, "True ||r||    = %.6e\n", true_rnorm);
        PetscPrintf(PETSC_COMM_WORLD, "||r||/||b||   = %.6e\n", true_rnorm / bnorm);
        PetscPrintf(PETSC_COMM_WORLD, "========================================\n");

        if (has_xref && rep == nreps - 1) {
            Vec xref;
            PetscPrintf(PETSC_COMM_WORLD, "\nReference solution comparison:\n");
            PetscCall(ReadVectorParallel(xref_path, &xref, PETSC_COMM_WORLD, nplanes));
            PetscReal xnorm, diff;
            PetscCall(VecNorm(xref, NORM_2, &xnorm));
            PetscCall(VecAXPY(xref, -1.0, x));
            PetscCall(VecNorm(xref, NORM_2, &diff));
            PetscPrintf(PETSC_COMM_WORLD, "||x_ref||         = %.6e\n", xnorm);
            PetscPrintf(PETSC_COMM_WORLD, "||x - x_ref||     = %.6e\n", diff);
            PetscPrintf(PETSC_COMM_WORLD, "Relative diff:     %.6e\n", diff / xnorm);
            PetscCall(VecDestroy(&xref));
        }

        if (nreps > 1)
            PetscCall(PetscLogStagePop());
    }

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    return PetscFinalize();
}
#endif /* CUDSS_NO_MAIN */
