# PETSc + cuDSS block-Jacobi solver - build and run on Perlmutter

A standalone driver that solves one linear system extracted from M3D-C1 with PETSc
LGMRES (`aijcusparse`/`cuda`) preconditioned by a **PCShell block-Jacobi** whose
per-block sub-solve is a **cuDSS LU** factorization. The matrix is block-structured by
toroidal plane; each plane's diagonal block is factored by cuDSS on the GPU(s) that own
it, and the off-plane coupling is handled by the outer Krylov iteration.

`run_config_b.sbatch` is the reference for every invocation - Slurm header, modules,
`LD_LIBRARY_PATH`, `srun` line and PETSc options.


**Package files:**
```
├── build.sbatch                  # split build on a GPU node
├── run_config_b.sbatch           # one GPU per plane, the reference invocation
├── run_mg.sbatch                 # several GPUs per plane, one process
├── run_mgmn.sbatch               # several ranks per plane
├── mg_gpu_wrap.sh                # disjoint-GPU srun launcher used by run_mg.sbatch
├── src/
│   ├── petsc_cudss_solve.c       # main driver
│   ├── petsc_cudss_solve_mg.c    # -mg_devices variant
│   ├── cudss_block.cu / .h       # extern "C" per-block cuDSS LU solver
│   └── cudss_commlayer_craympi.c # cuDSS MGMN comm layer on Cray-MPICH
├── deps/cudss/                   # vendored cuDSS 0.8.0.10
├── example_output/               # logs + profiles from a known-good build and run
└── data/                         # symlinks, one per option:
    ├── A.bin                     #   -A     matrix, 15.8 GiB
    ├── b.bin                     #   -b     right-hand side
    ├── x0.bin                    #   -x0    initial guess (optional)
    └── xout.bin                  #   -xref  reference solution
```

## 0. Quick Start

- Copy this folder with `cp -a` so you have write access.
- Adapt the account setting in the `*.sbatch` scripts.
- Build:
```
sbatch build.sbatch
```
- Run (after the build job finishes):
```
sbatch run_config_b.sbatch
sbatch run_mg.sbatch
```
- Compare results with those in `./example_output`.

## 1. Config

The scripts build in place and resolve `src/`, `data/` and `bin/` relative to the
directory you submit from, so **you need a writable copy of this directory** and then
submit from its root. If you copy this directory **use `cp -a`**, not `-r` to keep the
symlinks to the large data files.

**Point the `*.sbatch` scripts at YOUR Slurm account (they ship with -A m4642_g).**
All scripts use the `debug` queue, which allows only two jobs queued or running at
once - change as needed.

### 1.1 Dependencies

- cuDSS 0.8.0.10
	- in: `deps/cudss/`
	- override with `CUDSS_DIR`
- PETSc 3.23, GPU build
	- in: `/global/cfs/cdirs/mp288/jinchen/PETSC/production/petsc.20230606`, arch `perlmuttergpu-gnu`
	- override with: `PETSC_DIR`/`PETSC_ARCH` in `build.sbatch`
		- must be built `--with-cuda` against CUDA 12.4 with `aijcusparse`, for `sm_80`.
- Matrix and vector data
	- in:  `data/`, four symlinks into `/pscratch/sd/m/mnoack/m3d-c1/refs/m3dc1_bench_data/`:
		- `A.bin` (15.8 GiB),
		- `b|x0|xout.bin` (15 MB each)
	- override: pass different filenames via command line, preconditioner block size is computed by
	            matrix size divided by number of planes

Modules are loaded by the scripts; nothing to set up by hand. They are **pinned**
(`PrgEnv-gnu/8.5.0`, `cray-mpich/8.1.30`, `cray-hdf5-parallel/1.14.3.7`,
`cray-fftw/3.3.10.11`, `nccl/2.24.3`), so every job prints:
`The following have been reloaded with a version change:`

Two CUDA toolkits are used on purpose: `nvcc` **12.9** to compile (12.4's `nvcc` rejects
gcc-14 from `PrgEnv-gnu/8.5.0`), but `cudart`/`cusparse`/`cublas` at runtime must be the
NVHPC 24.5 / CUDA **12.4** ones PETSc was linked against. Mixing 12.9 cuBLAS into the
process fails with `cuBLAS EXECUTION_FAILED` at the first `VecNorm`.

## 2. Build

Run (on a GPU node):
```bash
sbatch build.sbatch
```

`nvcc` compiles `src/cudss_block.cu`, Cray `cc` (PETSc `pkg-config` flags) compiles and
links the drivers, `cc -shared` builds the Cray-MPICH comm layer. Creates `bin/`:

| binary                              | description                                  |
|-------------------------------------|----------------------------------------------|
| `bin/petsc_cudss_solve`             | main driver                                  |
| `bin/petsc_cudss_solve_mg`          | variant adding `-mg_devices`                 |
| `bin/libcudss_commlayer_craympi.so` | cuDSS MGMN comm layer over Cray-MPICH        |

alongside `bin/obj/` (objects) and `bin/cudss_link/`, which holds the unversioned
`libcudss.so` symlink the cuDSS wheel does not ship but the linker wants.

Check `build-<jobid>.out` for errors and for the
`=== ldd (verify cudss/nccl resolve) ===`
section, where `libcudss.so.0` must resolve to `deps/cudss/lib` (or `$CUDSS_DIR`).


## 3. Run

| script                | mode                                 | geometry                          |
|-----------------------|--------------------------------------|-----------------------------------|
| `run_config_b.sbatch` | one GPU per plane                    | 4 ranks, 1 node, 4 GPUs           |
| `run_mg.sbatch`       | multiple GPUs per plane, one process | 4 ranks, 2 nodes, 2 GPUs per rank |
| `run_mgmn.sbatch`     | multiple ranks per plane             | 8 ranks, 2 nodes, 1 GPU per rank  |

All scripts run in less than 5 minutes on Perlmutter.

Rows of the matrix are distributed plane-aligned, so `ranks_per_plane = nprocs / nplanes`
decides the mode; `nprocs % nplanes == 0` is required.
`-nplanes 4` must match the input matrix.
Each mode has its own script, and each writes its own `slurm-*` and `profile*` files so
they do not overwrite one another. Subsequent runs of the same script will, though.

## 3.1. Config: 1 node, 4 ranks, 4 GPUs, 1 GPU per rank and plane
```bash
sbatch run_config_b.sbatch
```

`ranks_per_plane == 1`: each rank owns one whole plane block and factors it on its
own GPU. No comm layer is involved. This is what `run_config_b.sbatch` does with
`-n 4 -nplanes 4` on 4 GPUs.

## 3.2. Config: 2 nodes, 4 ranks, 8 GPUs, 2 GPUs per rank and plane
```bash
sbatch run_mg.sbatch
```

`bin/petsc_cudss_solve_mg -mg_devices N` keeps one rank per plane but factors that
plane's block across `N` GPUs, inside one process, over NVLink via `cudssCreateMg`. No
comm layer is used, so **a rank's GPUs must all be on the rank's node**.

Two consequences the script encodes for you. Each rank must see a *disjoint* set of `N`
GPUs, which is what `mg_gpu_wrap.sh` does - it derives `CUDA_VISIBLE_DEVICES` from
`SLURM_LOCALID` and `$MG_DEVICES`, and refuses to start if
`MG_DEVICES × ntasks-per-node` exceeds the GPUs on a node. And since a node has 4 GPUs,
four planes at `MG_DEVICES=2` need 8 GPUs, so the ranks spread 2 per node across 2 nodes;
`MG_DEVICES=4` would need `--nodes=4 --ntasks-per-node=1`. Change `MG_DEVICES` in
`run_mg.sbatch` and the Slurm geometry together.

## 3.3. Config: 2 nodes, 8 ranks, 8 GPUs, 1 GPU per rank, 2 ranks per plane
```bash
sbatch run_mgmn.sbatch
```

`ranks_per_plane > 1`: cuDSS factors each plane block distributed across that plane's
ranks, one GPU per rank, and needs a comm layer.
**This part is work in progress within cuDSS**; it is included here, as it might be relevant
for the M3D-C1 integration. This mode would only be required if a single plane is distributed 
on more GPUs than a single node provides. So please treat this script as a preview on how
to execute and profile this distribution pattern, but not yet as a production-ready configuration.
The root cause of the long runtimes/slow convergence is a missing feature:

- **No distributed matching, yet.** cuDSS does not support matching in the distributed path;
  requesting it makes `CUDSS_PHASE_ANALYSIS` return `CUDSS_STATUS_NOT_SUPPORTED`, so the script
  sets `CUDSS_MGMN_MATCH=0`. For this matrix, matching is what puts large entries on the
  diagonal, and without it the block factor is inaccurate: the solve still converges, but
  the outer LGMRES needs roughly two orders of magnitude more iterations, which
  dominates the run time. (`CUDSS_PIVOT` and `CUDSS_IR` are the MGMN-legal
  alternatives to matching; neither recovers the iteration count.)

Perlmutter requires the following for MGMN:

- **Comm layer.** The NCCL comm layer shipped inside cuDSS is linked against OpenMPI
  (`libmpi.so.40`), which does not exist on Perlmutter, so it cannot be loaded here and
  `cudssSetCommLayer` rejects it. The script therefore uses the Cray-MPICH shim built
  from `src/cudss_commlayer_craympi.c` (`CUDSS_COMM_BACKEND=mpi` plus
  `-cudss_comm_lib`), and passes an **absolute** path because cuDSS `dlopen()`s it.
  This is why `-cudss_comm_lib` has no useful default: the built-in one names the NCCL
  layer.

## 4. Output [AI generated]

In the scripts, the benchmark is repeated twice, once for warm-up. The timings of the second
run are representative.

```
=== PETSc Standalone Solver ===
MPI ranks: 4
Plane-aligned distribution: 4 planes, 1 ranks/plane         # planes must match the matrix
...
Reading matrix: data/A.bin
  1879920 x 1879920, nnz = 1409524416                       # 4 planes x 469980 rows
  Rank 0: rows [0, 469980), local_nnz = ...
Loading RHS: data/b.bin
Loading x0: data/x0.bin
  Data loaded in 14.0 s
||b|| = 5.930404e-03
Plane size: 469980 rows/plane
cuDSS block-Jacobi PC: 1 ranks/plane -> single-GPU/plane    # which mode is running
  plane block extracted ... in 1.4 s
  cuDSS setup: ANALYSIS ... ms, FACTORIZATION ... ms (wall ...)
```

`Plane-aligned distribution` tells you `-nplanes` took effect and with how many ranks per
plane; `cuDSS block-Jacobi PC: ... ->` tells you which kind of run you are actually in.
If either disagrees with what you intended, please stop here, as everything below will
still look healthy.

The output has one `RESULTS (rep i/N)` block per rep - `Reason` is PETSc's `KSPConvergedReason`
(`2` is `KSP_CONVERGED_RTOL`, negative means diverged), then iteration count, load and
solve times, and residual norms - each followed by a `cuDSS timing (max/rank)` line.
After the **last** rep, a single `Reference solution comparison` block reports
`||x_ref||`, `||x - x_ref||` and `Relative diff` against `-xref`. A healthy run converges
with the same iteration count in every rep, is slower in rep 1 than later ones (warmup),
and has a small `Relative diff`.

In the `cuDSS timing` line, **only** `SOLVE`, the apply count and the per-apply cost are
per-rep. `ANALYSIS` and `FACTOR` are the one-time setup costs incurred before rep 1 and
are reprinted unchanged in every rep, so do not read them as work rep 2 did - that would
contradict its much shorter `Solve time`. With `-refactor_each_rep`, `FACTOR` does
become per-rep. The apply count exceeds the iteration count by a couple (27 applies for 25
iterations) because LGMRES also applies the preconditioner outside the counted iterations,
for the initial residual and the augmentation vectors.

**`||r||/||b||` is not what `-ksp_rtol` tests.** With `-x0` supplied and a left
preconditioner, PETSc's convergence test compares the *preconditioned* residual against
`rtol` times its value at the initial guess; `KSP ||r||` is that preconditioned norm,
while `True ||r||` is the recomputed `||Ax - b||`. On a correctly converged run the
printed `||r||/||b||` is therefore orders of magnitude larger than `rtol`. Judge accuracy
by `Relative diff` against the supplied reference solution, and by `True ||r||`.

Console output goes to `slurm-<jobid>.out`, `slurm-mg-<jobid>.out` or`slurm-mgmn-<jobid>.out`

`-log_view` additionally writes a `.txt`, a `.fg` (flamegraph) and a `.xml` profile into
the submit directory: `profile.*`, `profile_mg.*` or `profile_mgmn.*`.
The names are **fixed per script**, so rerunning the same script overwrites
its own profiles. With `-nreps > 1` the last rep is timed separately in the `Benchmark`
stage of the `.txt` profile.

## 5. Details [AI generated]

### 5.1 Using a different inputs

`data/*.bin` are unmodified PETSc binary dumps of the M3D-C1 velocity (`hard_`) solve, so
another case needs no conversion - dump it from the M3D-C1 run itself:

```bash
srun ... m3dc1_3d -hard_ksp_view_mat      binary:A.bin \
                  -hard_ksp_view_rhs      binary:b.bin \
                  -hard_ksp_view_solution binary:xout.bin
```

Then point `-A`/`-b`/`-xref` at the new files and set `-nplanes` to that case's `nplanes`
from `C1input`. That is the only thing you must change: there is no block-size setting,
because the `-matload_block_size 36` in PETSc's `.info` sidecar is the per-node DOF
blocking (12 DOFs/node × 3 velocity components), which this driver ignores and which
always divides the plane size anyway. `-x0` is optional; the shipped `x0.bin` is a zero
vector.

The dump must satisfy what the driver's reader assumes, none of which it can check for
you: written by a PETSc with **32-bit indices and real double** scalars, square, total
`nnz` below 2^31, and rows ordered **plane-major** in equal contiguous planes - M3D-C1
already numbers rows that way, which is what makes the plane diagonal blocks extractable.

### 5.2 Driver options

Everything else on the command line is ordinary PETSc - see `run_config_b.sbatch` for
the `-ksp_*`, `-mat_type`, `-vec_type` and `-log_view` settings actually used.

There is deliberately **no `-pc_type`**: the driver installs the cuDSS block-Jacobi as a
`PCShell` in code, so the preconditioner is not selectable from the command line. Swapping
in a different preconditioner for comparison means editing `src/petsc_cudss_solve.c`.

- `-A <f>` - no default, required
    - matrix, PETSc binary format.
- `-b <f>` - no default, required
    - right-hand side, PETSc binary format.
- `-x0 <f>` - no default, optional
    - initial guess, PETSc binary format. Without it every rep starts from a zero vector.
- `-xref <f>` - no default, optional
    - reference solution, PETSc binary format. Used only after the last rep, for the
      `Reference solution comparison` block and its `Relative diff`; without it that block
      is not printed.
- `-nplanes <N>` - default `0`, required, and not a tuning knob
    - It must equal the number of toroidal planes in the matrix you pass - 4 for the
      shipped `data/A.bin` - because the block-Jacobi blocks *are* the plane diagonal blocks.
    - It must also be > 0 and divide the rank count, and every configuration in §3 keeps
      `-nplanes 4` while varying `-n`.
    - A value that divides the rank count but does not match the matrix cannot be detected:
      the run converges and prints a plausible iteration count, but it preconditioned with
      something that is not the plane blocks, so the numbers mean nothing. Confirm the two
      lines named in §4 before trusting any output.
- `-nreps <N>` - default `1`
    - repeat the solve. Rep 1 is warmup and the last rep gets its own `Benchmark` log stage,
      so use `-nreps 2` or more for any timing you intend to quote.
- `-refactor_each_rep` - default off
    - `petsc_cudss_solve` only, silently ignored by `petsc_cudss_solve_mg`: numeric
      refactorization before each rep, symbolic analysis reused.
- `-cudss_comm_lib <f>` - default cuDSS-internal NCCL layer, unusable on Perlmutter
    - comm layer, MGMN only; must be an absolute path (§3.3).
- `-mg_devices <N>` - default `1`
    - `petsc_cudss_solve_mg` only: GPUs per plane in MG mode (§3.2). `run_mg.sbatch` sets it
      to 2 via `$MG_DEVICES`.
    - That driver also ignores `-refactor_each_rep` and does not forward
      `-mat_cusparse_spmv_alg` to the sub-blocks, so use `petsc_cudss_solve` for
      refactor-cadence and determinism work.

### 5.3 Environment variables

Set by the batch scripts:

- `CUDSS_DIR` - cuDSS install to build and run against
    - Read by the scripts, not set by them, defaulting to `deps/cudss`.
    - If you override it, pass it to *both* jobs (`CUDSS_DIR=... sbatch build.sbatch` and
      `CUDSS_DIR=... sbatch run_config_b.sbatch`), or the run will load a different cuDSS
      than the build linked against.
- `LD_LIBRARY_PATH` - NVHPC 24.5 / CUDA 12.4 first (for `libnvToolsExt.so.1`), then cuDSS
  and NCCL
- `MPICH_GPU_SUPPORT_ENABLED=1` - GPU-aware Cray-MPICH
- `OMP_NUM_THREADS=1`, `SLURM_CPU_BIND=cores` - one thread per rank, pinned

Read by the solver. All are **opt-in**: unset, cuDSS runs the configuration the scripts
were written for, except where a script sets one itself (noted below). Any that you do set
echoes a `[cudss_block] ...` line to stderr.

- `CUDSS_DETERMINISTIC` = `1`
    - bitwise-reproducible factorization and solve, at a performance cost
- `CUDSS_COMM_BACKEND` = `nccl` (default), `mpi`
    - MGMN comm layer family; `run_mgmn.sbatch` sets `mpi` (§3.3)
- `CUDSS_SHIM_TRACE` = `1`
    - trace every Cray-MPICH comm-layer call to stderr; affects only
      `libcudss_commlayer_craympi.so` (§3.3)
- `CUDSS_MG_DEVICES` = `N`
    - GPUs per block in MG mode; `-mg_devices` sets this for you.
    - Not to be confused with `MG_DEVICES`, which `run_mg.sbatch` uses to drive both
      `-mg_devices` and `mg_gpu_wrap.sh`.
- `CUDSS_NO_MATCH` = `1`
    - disable MC64 matching on the single-GPU/MG path
- `CUDSS_MGMN_MATCH` = `1` (default), `0`
    - matching on the MGMN path; `1` fails, so `run_mgmn.sbatch` sets `0` (§3.3)
- `CUDSS_REORDER` = `default`, `nd`, `amd`, `colamd`, `btf_colamd`, `none`
    - fill-reducing ordering; MGMN defaults to `nd`
- `CUDSS_PIVOT` = `none` (default), `auto`, `local`, `gcol`, `grow`
    - pivoting strategy
- `CUDSS_PIVOT_THRESHOLD` = double
    - pivoting aggressiveness, only with an active strategy
- `CUDSS_PIVOT_EPS` = double
    - replacement value for tiny pivots
- `CUDSS_PIVOT_EPS_ALG` = `0` default, `1` scaled, `2` static
    - absolute or norm-relative replacement
- `CUDSS_IR` = `N`
    - `N` steps of iterative refinement per solve (default 0)
- `CUDSS_FACT_ALG` = `default`, `multiblock`, `general`
    - factorization variant
- `CUDSS_BLOCK_HYBRID` = `1`
    - MGMN only: spill factors to host memory to avoid GPU OOM

PETSc's `-mat_cusparse_spmv_alg` (`CSRMV_ALG2` is the reproducible one) is forwarded by
`petsc_cudss_solve` to the `MPIAIJCUSPARSE` sub-blocks, so it also takes effect on more
than one rank.

