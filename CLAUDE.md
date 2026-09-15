# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Hard rules

- **Never compile code (run `make`, or anything that builds `unstructured` or `m3dc1_scorec`) without the
  user's explicit permission first.** Builds on this repo are expensive (many `ARCH`/flag variants, HPC
  compute) — always ask before running a build command, even to "check that it compiles."
- **Never submit jobs** (e.g. `sbatch`, `salloc`, `srun`, or the `regtest/run` script, which submits batch
  jobs). These consume shared HPC allocation and must only be launched by the user.

## What this is

M3D-C1 is a parallel, unstructured-mesh extended-MHD code for simulating tokamak/stellarator plasmas
(finite elements, implicit time advance). The physics/solver code is almost entirely Fortran 90
(`unstructured/*.f90`), built on top of `m3dc1_scorec` — a C++ wrapper around the SCOREC `core` mesh/FE
library (PUMI, plus PETSc for linear algebra) that is vendored as a subdirectory and built as a separate
CMake target.

Repo layout:
- `unstructured/` — the M3D-C1 physics code itself (this is where almost all development happens).
  Built with its own hand-written `make`-based system, *not* the top-level CMake.
- `m3dc1_scorec/` — C++/Fortran interface layer to the SCOREC unstructured-mesh library. Built via CMake,
  normally only touched when changing the mesh/linear-algebra interface.
- `CMakeLists.txt` (top level) — thin wrapper that just adds both subdirectories; day-to-day builds do not
  go through this.
- `skeleton/` — minimal example/template.
- `doc/` — LaTeX source for the physics/numerics documentation (`M3D-C1.pdf`).

## Building `unstructured` (the code you'll normally be changing)

Builds are driven entirely by `make` from the `unstructured/` directory, keyed off an `ARCH` (a.k.a.
`M3DC1_ARCH`) that selects a `$(ARCH).mk` machine config file in that same directory (e.g.
`perlmutter_cpu.mk`, `cori.mk`, `stellar.mk`; NERSC Perlmutter is the primary machine used in this
checkout). Per-machine setup notes live in `unstructured/README/readme.<system>`.

```
cd unstructured
module use $M3DC1_CODE_DIR/unstructured/modules/perlmutter
module load m3dc1/devel-cpu          # or m3dc1/devel-cpu-gcc

# build one variant (each combination of flags produces a separate _<arch>-<flags> object dir)
make OPT=1 MAX_PTS=25 ARCH=perlmutter_cpu                       # 2D real
make OPT=1 COM=1 MAX_PTS=25 ARCH=perlmutter_cpu                 # 2D complex (linear stability)
make OPT=1 3D=1 MAX_PTS=60 ARCH=perlmutter_cpu                  # 3D real
make OPT=1 3D=1 MAX_PTS=125 ARCH=perlmutter_cpu ST=1            # 3D stellarator

# or build every standard variant at once
make ARCH=perlmutter_cpu all

# assemble all built executables + SCOREC mesh utilities into unstructured/_<ARCH>/bin
make bin ARCH=perlmutter_cpu
```

Key make flags: `OPT=1` (optimized vs. debug build), `3D=1` (3D vs. 2D), `COM=1` (complex arithmetic, used
for linear/eigenvalue runs), `ST=1` (stellarator geometry, implies 3D), `PAR=1` (kinetic particle module),
`OMP=1` / `ACC=1` (OpenMP / GPU offload), `MAX_PTS=N` (max quadrature sampling points; use 25 for 2D, 60 for
3D, 125 for 3D+ST). Each combination builds into its own `unstructured/_$(ARCH)<flags>-<opt>-<MAX_PTS>/`
object directory and produces a distinctly-named binary (`m3dc1_2d`, `m3dc1_2d_complex`, `m3dc1_3d`,
`m3dc1_3d_st`, ...). `make clean ARCH=<arch>` removes the object dir for one variant;
`make cleanall ARCH=<arch>` removes all variants for that arch.

There is no separate lint/format/unit-test tooling — correctness is validated by the regression suite below.

## Regression tests

```
cd unstructured
export M3DC1_MPIRUN=srun M3DC1_VERSION=local M3DC1_ARCH=perlmutter_cpu
make bin ARCH=$M3DC1_ARCH
export PATH=$PWD/_$M3DC1_ARCH/bin:$PATH
cd regtest
./run $M3DC1_ARCH          # submits every test as a batch job
./check $M3DC1_ARCH         # compares C1ke output against regtest/<test>/base/C1ke

./run $M3DC1_ARCH pellet    # run/check a single named test (see subdirectory names in regtest/)
./check $M3DC1_ARCH pellet
./clean $M3DC1_ARCH         # remove generated test run directories
```

Each test works by copying `regtest/<test>/base/` to a new run directory, executing the code there, and
diffing the resulting `C1ke` (energy/diagnostic) file against `base/C1ke`. Test scripts assume
`m3dc1_2d[_complex]`, `m3dc1_3d[_st]`, and the SCOREC mesh utilities are on `$PATH` (satisfied by the
`make bin` step above). If a change legitimately alters `C1ke` output, overwrite `base/C1ke` with the new
result and add an explanatory entry to `regtest/CHANGELOG`.

## Architecture notes (unstructured/)

- **Physics/discretization**: `metricterms_new.f90` assembles the weak-form finite-element operators;
  `ludef_t.f90` builds the linear systems solved each timestep; `time_step*.f90` drive the (split/unsplit)
  implicit time advance. `M3Dmodules.f90`-family global state modules are shared broadly across the code.
- **Equilibrium/initial conditions**: each `init_*.f90` file is a self-contained initial-condition/equilibrium
  generator (e.g. `init_eqdsk.f90` reads G-EQDSK, `init_solovev.f90` is analytic Solov'ev, `init_vmec.f90`
  reads VMEC stellarator equilibria); `gradshafranov.f90` solves the GS equation for 2D equilibria.
- **External fields/coils**: `coils.f90` / `coil_sets.f90` model external coil geometry and Biot-Savart
  fields; `rmp.f90` / `read_schaffer_field.f90` handle resonant magnetic perturbations; `fit_magnetics.f90`
  fits magnetic diagnostics.
- **Mesh/solver abstraction**: `field.o`, mesh, vector, and matrix modules are selected at compile time via
  the `mesh_mod`/`vector_mod`/`matrix_mod` preprocessor macros in `unstructured/makefile` (SCOREC backend
  by default; a PETSc-only backend exists via `USEPETSC=1`) — this indirection is why grepping for e.g.
  `scorec_matrix` won't show all callers of the matrix interface.
- **Kinetic/atomic physics extensions**: `particle.f90` (kinetic particle module, `PAR=1`), `kprad*.f90`
  (KPRAD atomic/radiation model), `adas_m3dc1.f90` (ADAS atomic data, `ADAS=1`), `bootstrap.f90` /
  `gyroviscosity.f90` (closure terms), `hot_tail.f90` / `runaway*.f90` (runaway-electron physics).
- **I/O**: `hdf5_output.f90` / `restart_hdf5.f90` write simulation output and restart files (HDF5 is the
  on-disk format); `iterdb.f90` / `readgeqdsk.f90` read external transport/equilibrium files.
- Post-processing/analysis lives outside `unstructured/`: `unstructured/idl/*.pro` (IDL analysis scripts,
  the traditional M3D-C1 analysis toolchain) and `unstructured/python/m3dc1/` (Python analysis package).

## Working conventions

- Source is fixed/free-form Fortran 90 compiled with `-fdefault-real-8` on most machine configs (i.e. bare
  `real` is 8 bytes) — match existing real-kind conventions in a file rather than introducing explicit kind
  parameters ad hoc.
- New physics/config options are typically threaded through as: a namelist variable in
  `unstructured/input.f90` (or the relevant `init_*`/physics module) → used in the relevant `.f90` physics
  routine → documented in `doc/`. Follow this pattern rather than hardcoding new behavior.
- Machine-specific build settings belong in the relevant `unstructured/<ARCH>.mk`, never inline in
  `unstructured/makefile`, which is meant to stay machine-independent.
