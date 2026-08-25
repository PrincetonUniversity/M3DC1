// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * PETSc standalone solver for the M3D-C1 extracted matrix (matrix_id=5),
 * using cuDSS as the per-plane block-Jacobi sub-solver via a PCShell.
 *
 * MG VARIANT of petsc_cudss_solve.c: it adds -mg_devices N, which factors each
 * plane block across N GPUs inside a SINGLE process (cuDSS single-process
 * multi-GPU over NVLink, via cudssCreateMg).  That path uses no comm layer and
 * is therefore single-node only.  Everything else matches petsc_cudss_solve.c
 * except that this variant deliberately does NOT carry the two later additions
 * made there: -refactor_each_rep and the -mat_cusparse_spmv_alg forwarding to
 * the MPIAIJCUSPARSE sub-blocks.  Use petsc_cudss_solve for the cadence and
 * determinism measurements; use this one only to explore >1 GPU per plane.
 *
 * The Krylov outer solve (LGMRES, aijcusparse/cuda) is that of the
 * SuperLU_DIST baseline driver petsc_standalone_solve.c, which is not part of
 * this package; the preconditioner is replaced by a custom PCShell whose Apply
 * runs a cuDSS LU solve on this rank's piece of its plane diagonal block.
 * With one rank per plane the solve is single-GPU (or MG, above); when a plane
 * spans >1 rank cuDSS runs in MGMN mode over the plane's sub-communicator,
 * using the NCCL comm layer or the cray-mpich shim in
 * cudss_commlayer_craympi.c (CUDSS_COMM_BACKEND=mpi).  MG and MGMN are
 * mutually exclusive.  See cudss_block.{h,cu}.
 *
 * Requires plane-aligned distribution: nplanes>0, nprocs%nplanes==0.
 *
 * Mat/Vec types are controlled by -mat_type / -vec_type (use aijcusparse/cuda).
 *
 * Build: split nvcc (.cu) + Cray cc (.c) + link; see build.sbatch.
 *
 * Usage (4 planes on 1 node, 1 rank/plane, 1 GPU/plane -- raise -mg_devices
 * only with a launcher that gives each rank a disjoint GPU set, see below):
 *   srun -n 4 ./petsc_cudss_solve_mg \
 *     -A A.bin -b b.bin [-x0 x0.bin] [-xref xout.bin] -nplanes 4 \
 *     -mg_devices 1 \
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
 *                       earlier ones in "Warmup".
 *   -mg_devices <N>     GPUs per plane block in single-process MG mode
 *                       (default 1 = plain single-GPU).  Forwarded to
 *                       cudss_block as CUDSS_MG_DEVICES.
 *   -cudss_comm_lib <f> cuDSS comm-layer .so; MGMN only (>1 rank/plane).
 *                       Defaults to the vendored NCCL layer at
 *                       ${CUDSS_DIR:-deps/cudss}/lib/libcudss_commlayer_nccl.so.0.
 *                       For the cray-mpich shim instead: CUDSS_COMM_BACKEND=mpi
 *                       with -cudss_comm_lib bin/libcudss_commlayer_craympi.so.
 *
 * cuDSS is further tuned through CUDSS_* environment variables (CUDSS_REORDER,
 * CUDSS_PIVOT, CUDSS_DETERMINISTIC, ...) -- see cudss_block.cu.
 */
#include <petsc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include "cudss_block.h"

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
    PetscFunctionBeginUser;            /* block built once in main(); nothing to do */
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

static PetscErrorCode CudssPCDestroy(PC pc)
{
    PetscFunctionBeginUser;            /* block destroyed in main() */
    PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL,
        "PETSc standalone solver for M3D-C1 extracted matrix\n"
        "Usage: mpirun -np N petsc_standalone_solve -A mat.bin -b rhs.bin [options]\n"));

    /* NOTE: thresholding is controlled via the -log_threshold 0.0 CLI option
     * (parsed inside PetscInitialize, gated on `if (any_nested)`).  Calling
     * PetscLogSetThreshold here unconditionally would SEGV on PETSc 3.20+
     * for plain -log_view runs (no nested handler installed). */

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

    /* MG variant: each plane block (1 rank/plane) is factored by cuDSS across
     * -mg_devices GPUs in a single process (single-process multi-GPU, NVLink,
     * NO comm layer -> single-node only).  We forward this to cudss_block via
     * the CUDSS_MG_DEVICES env var (read in its single-GPU path).  The launcher
     * MUST give each rank a disjoint CUDA_VISIBLE_DEVICES set of mg_devices
     * GPUs, via an srun wrapper script keyed on SLURM_LOCALID, so that
     * co-located ranks do not collide.  No such wrapper ships in this package;
     * with the default -mg_devices 1 none is needed. */
    PetscInt mg_devices = 1;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-mg_devices", &mg_devices, NULL));
    if (mg_devices > 1 && !getenv("CUDSS_MG_DEVICES")) {
        char mgstr[16];
        snprintf(mgstr, sizeof(mgstr), "%d", (int)mg_devices);
        setenv("CUDSS_MG_DEVICES", mgstr, 1);
    }

    /* Path to the cuDSS comm-layer .so (used only for MGMN, i.e. when a plane
     * spans >1 rank; unused in MG or single-GPU mode).  Defaults to the NCCL
     * layer that ships with the vendored cuDSS, following the same
     * ${CUDSS_DIR:-deps/cudss} convention as build.sbatch and the run scripts;
     * the relative fallback assumes the package root is the cwd, as it is for
     * the -A data/A.bin paths.  Use CUDSS_COMM_BACKEND=mpi plus
     * -cudss_comm_lib bin/libcudss_commlayer_craympi.so for the cray-mpich
     * shim built from cudss_commlayer_craympi.c. */
    char comm_lib[1024];
    PetscBool has_comm_lib;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-cudss_comm_lib", comm_lib, sizeof(comm_lib), &has_comm_lib));
    if (!has_comm_lib) {
        const char *cudss_dir = getenv("CUDSS_DIR");
        snprintf(comm_lib, sizeof(comm_lib), "%s/lib/libcudss_commlayer_nccl.so.0",
                 (cudss_dir && cudss_dir[0]) ? cudss_dir : "deps/cudss");
    }

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
    if (has_x0) {
        PetscPrintf(PETSC_COMM_WORLD, "Loading x0: %s\n", x0_path);
        Vec x0;
        PetscCall(ReadVectorParallel(x0_path, &x0, PETSC_COMM_WORLD, nplanes));
        PetscCall(VecCopy(x0, x));
        PetscCall(VecDestroy(&x0));
    }

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

    /* ---- cuDSS block-Jacobi preconditioner setup (once) ---- */
    if (nplanes <= 0 || nprocs < nplanes || (nprocs % nplanes) != 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "petsc_cudss_solve requires plane-aligned layout: nplanes>0 and nprocs%%nplanes==0 (got nprocs=%d, nplanes=%d)",
                (int)nprocs, (int)nplanes);
    }
    PetscInt n_global;
    PetscCall(MatGetSize(A, &n_global, NULL));
    /* -nplanes must be the toroidal plane count OF THIS MATRIX, not a free knob:
     * the block-Jacobi blocks are the plane diagonal blocks.  A divisible but
     * wrong value cannot be detected in general -- it just preconditions with
     * something else and still converges -- so catch at least the indivisible
     * case rather than silently truncating rows. */
    PetscCheck(n_global % nplanes == 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
               "-nplanes %d does not divide the matrix dimension %d; -nplanes must equal the "
               "number of toroidal planes in the matrix (4 for the A.bin shipped with this package)",
               (int)nplanes, (int)n_global);
    PetscInt rows_per_plane  = n_global / nplanes;
    PetscInt ranks_per_plane = nprocs / nplanes;
    PetscInt plane           = rank / ranks_per_plane;
    PetscInt plane_start     = plane * rows_per_plane;
    PetscInt plane_end       = plane_start + rows_per_plane;

    PetscInt local_start, local_n;
    compute_row_range(n_global, rank, nprocs, nplanes, &local_start, &local_n);

    MPI_Comm plane_comm;
    MPI_Comm_split(PETSC_COMM_WORLD, (int)plane, (int)(rank % ranks_per_plane), &plane_comm);

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

    PetscPrintf(PETSC_COMM_WORLD,
        "cuDSS block-Jacobi PC: %d ranks/plane -> %s\n",
        (int)ranks_per_plane, (ranks_per_plane > 1) ? "MGMN (distributed cuDSS)"
                              : (mg_devices > 1) ? "MG single-process multi-GPU/plane"
                                                 : "single-GPU/plane");
    if (ranks_per_plane > 1)
        PetscPrintf(PETSC_COMM_WORLD, "  comm layer: %s\n", comm_lib);
    else if (mg_devices > 1)
        PetscPrintf(PETSC_COMM_WORLD, "  MG devices/plane: %d (CUDA_VISIBLE_DEVICES=%s)\n",
                    (int)mg_devices, getenv("CUDA_VISIBLE_DEVICES") ? getenv("CUDA_VISIBLE_DEVICES") : "(unset)");

    int *blk_rowptr, *blk_colidx;
    double *blk_vals;
    long long blk_nnz;
    PetscLogDouble t_ext0, t_ext1;
    PetscCall(PetscTime(&t_ext0));
    PetscCall(ExtractPlaneBlock(A, plane_start, plane_end, local_start, local_n,
                                &blk_rowptr, &blk_colidx, &blk_vals, &blk_nnz));
    PetscCall(PetscTime(&t_ext1));
    PetscPrintf(PETSC_COMM_WORLD, "  plane block extracted (rank0 local_n=%d, local_nnz=%lld) in %.2f s\n",
                (int)local_n, blk_nnz, t_ext1 - t_ext0);

    CudssPCCtx pcctx;
    PetscLogDouble t_setup0, t_setup1;
    PetscCall(PetscTime(&t_setup0));
    pcctx.blk = cudss_block_create(plane_comm, (int)rows_per_plane, (int)local_n,
                                   (int)(local_start - plane_start),
                                   blk_rowptr, blk_colidx, blk_vals, blk_nnz, comm_lib);
    PetscCall(PetscTime(&t_setup1));
    PetscCall(PetscFree(blk_rowptr));
    PetscCall(PetscFree(blk_colidx));
    PetscCall(PetscFree(blk_vals));

    {
        double an, fa, so; int ns;
        cudss_block_get_timers(pcctx.blk, &an, &fa, &so, &ns);
        PetscPrintf(PETSC_COMM_WORLD,
            "  cuDSS setup: ANALYSIS %.1f ms, FACTORIZATION %.1f ms (wall %.2f s)\n\n",
            an, fa, t_setup1 - t_setup0);
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
        } else {
            PetscCall(VecZeroEntries(x));
        }

        cudss_block_reset_solve_timer(pcctx.blk);

        KSP ksp;
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
        PetscCall(KSPSetOperators(ksp, A, A));
        PetscCall(KSPSetFromOptions(ksp));
        if (has_x0)
            PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));

        /* Install our cuDSS block-Jacobi as a PCShell (overrides -pc_type). */
        {
            PC pc;
            PetscCall(KSPGetPC(ksp, &pc));
            PetscCall(PCSetType(pc, PCSHELL));
            PetscCall(PCShellSetContext(pc, &pcctx));
            PetscCall(PCShellSetSetUp(pc, CudssPCSetUp));
            PetscCall(PCShellSetApply(pc, CudssPCApply));
            PetscCall(PCShellSetDestroy(pc, CudssPCDestroy));
            PetscCall(PCShellSetName(pc, "cudss_bjacobi"));
        }

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
        PetscPrintf(PETSC_COMM_WORLD, "Reason:       %d\n", reason);
        PetscPrintf(PETSC_COMM_WORLD, "Iterations:   %d\n", its);
        PetscPrintf(PETSC_COMM_WORLD, "Load time:    %.3f s\n", t_load - t_start);
        PetscPrintf(PETSC_COMM_WORLD, "Solve time:   %.3f s\n", t1 - t0);
        PetscPrintf(PETSC_COMM_WORLD, "||b||         = %.6e\n", bnorm);
        PetscPrintf(PETSC_COMM_WORLD, "KSP ||r||     = %.6e\n", rnorm);
        PetscPrintf(PETSC_COMM_WORLD, "True ||r||    = %.6e\n", true_rnorm);
        PetscPrintf(PETSC_COMM_WORLD, "||r||/||b||   = %.6e\n", true_rnorm / bnorm);
        PetscPrintf(PETSC_COMM_WORLD, "========================================\n");

        {   /* cuDSS preconditioner-apply timing (max across ranks) */
            double an, fa, so; int ns;
            cudss_block_get_timers(pcctx.blk, &an, &fa, &so, &ns);
            double an_max, fa_max, so_max;
            int ns_max;
            MPI_Reduce(&an, &an_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
            MPI_Reduce(&fa, &fa_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
            MPI_Reduce(&so, &so_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
            MPI_Reduce(&ns, &ns_max, 1, MPI_INT,    MPI_MAX, 0, PETSC_COMM_WORLD);
            PetscPrintf(PETSC_COMM_WORLD, "cuDSS timing (max/rank): ANALYSIS %.1f ms, FACTOR %.1f ms, "
                        "SOLVE %.1f ms (%d applies, %.2f ms/apply)\n",
                        an_max, fa_max, so_max, ns_max, ns_max ? so_max/ns_max : 0.0);
            PetscPrintf(PETSC_COMM_WORLD, "========================================\n");
        }

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

        PetscCall(KSPDestroy(&ksp));

        if (nreps > 1)
            PetscCall(PetscLogStagePop());
    }

    cudss_block_destroy(pcctx.blk);
    MPI_Comm_free(&plane_comm);
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    return PetscFinalize();
}
