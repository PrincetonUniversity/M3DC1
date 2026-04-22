# AGENTS.md

## Project Scope
- This repository ports IDL `.pro` routines to Python under the `m3dc1` module.
- Keep Python file/function structure close to IDL call structure (1-to-1 where practical).

## User Preferences (Persistent)
- Use Homebrew Python 3.14 for commands:
  - `/opt/homebrew/bin/python3.14`
- Implement module-style code only:
  - No `main` function / no script-only entry requirement.
  - Must support `import m3dc1.<module>` usage.
- When IDL calls a subroutine/function, create a separate Python file for it and call it from other modules.
- Reuse existing Python implementations if already present (do not reimplement duplicates), e.g. `get_normalizations`, `convert_units`.
- For missing field names in HDF5 reads:
  - `read_field.py` should raise `KeyError`.
  - `plot_field.py` should raise `KeyError` when underlying field is missing.

## Plot/Return Behavior
- `plot_mesh.py` should return `(figure, axis)`.
- `plot_field.py` should return `(figure, axis)`.
- `plot_flux_average.py` should return `(figure, axis)`.
- Support `iso` option in plotting functions via `ax.set_aspect('equal')`.
- Default `contourf` levels should be `100`.
- In `contour_and_legend.py`, use:
  - `plasma` colormap if data is all positive or all negative.
  - `turbo` colormap if data crosses zero.

## String/Label Conventions
- Remove IDL text formatting escapes like `!X`, `!6` in Python label strings.
- Convert IDL-style math labels to LaTeX-like strings directly in-place (no extra wrapper function), using `$...$` where needed.
- Preserve intended subscript semantics (e.g. `I_{tot}`), not flattened alternatives.
- Apply similar label conversion for `read_scalar.py`, `read_field.py`, and unit parsing labels.

## Global Migration Rules
- Preserve IDL-style logging/print behavior across all migrated files.
- This rule applies to all future IDL-to-Python migration work in this repository unless explicitly overridden by the user.

## read_field / eval_field Notes
- `eval_field.py` is currently Numba-based (`njit` + `prange`) for element-loop acceleration.
- Worker control is environment-driven:
  - `M3DC1_EVAL_WORKERS` controls requested threads.
  - Clamp to valid Numba thread range internally.

## Performance Testing Conventions
- Prefer direct benchmark comparisons with fixed `points`, reporting:
  - runtime for each worker setting
  - speedup ratio
  - max absolute difference between outputs

## Geometry Mapping Notes (Persistent)
- `read_field.py` now implements `igeometry=1` mapping behavior:
  - For primitive reads, when `igeometry > 0` and `logical=False`, it reads `rst`/`zst`, builds map indices with `create_map`, and remaps data with `map_field`.
  - `logical=True` must bypass geometry remapping.
- IDL subroutine split rule is required here:
  - Keep mapping logic in separate modules: `create_map.py` and `map_field.py`, called from `read_field.py`.
- `create_map.py` implementation requirement:
  - Keep the Numba-decorated implementation.
  - It must still run when Numba is unavailable by using a no-op `njit` fallback (same function body executes in Python mode).
- `map_field.py` implementation requirement:
  - Use `scipy.interpolate.RegularGridInterpolator` (not custom `_interp2`).
  - Use vectorized interpolation calls (e.g., interpolate all `(mx,my)` points at once), not per-point Python-loop interpolation calls.

## Session Decisions (Persistent)
- `eval_field.py`:
  - Keep Numba thread chunksize fixed at `16` via `set_parallel_chunksize(16)`.
  - Keep `debug_spin` support active in `eval_field` and kernel path.
  - In the `debug_spin` loop, call `_eval_poly_numba(...)` (current debug workload behavior).
- `read_field.py` / `plot_field.py` interface updates:
  - Use `timeslices` argument name (not `slices`/`time`).
  - `read_field.py` no longer uses `x,y,t` positional arguments.
  - `plot_field.py` no longer uses `x,y` arguments and uses `timeslices`.
- Multi-timeslice behavior:
  - `read_field.py` supports multiple `timeslices` including `diff` behavior.
  - `plot_field.py` should pass multi-`timeslices` through to `read_field.py`.
- Numba runtime/testing conventions used in this repo:
  - For Numba tests in this environment, use:
    - `source /etc/profile`
    - `module load py-numba`
  - Common benchmark settings used in this session include larger `points` (e.g., `1000`/`10000`) and worker comparisons (`workers=1` vs `workers=4`).
- Plot mesh remap details:
  - In `plot_mesh.py`, when handling `rst_meta`/`zst_meta`, fill masked points using 8-neighbor averaging.
  - Fill target condition is `mask[i, j] == 1` (use `out_mask` from meta).
  - Neighbor inclusion condition is `mask[ii, jj] != 1` (with finite-value check).
- `plot_mesh.py` interface update:
  - Add `phi` argument and pass it to `read_field("rst", ...)` / `read_field("zst", ...)`.
- `plot_field.py` interface update:
  - Add explicit `logical` argument in `plot_field(...)`.
  - Pass `logical` through to `read_field(...)`.
  - Pass `logical` through to all internal `plot_mesh(...)` calls.
- Magnetic probe plotting naming:
  - Python entry is `plot_mag_probes(...)` (plural), corresponding to IDL `plot_mag_probes.pro` in the upper-level IDL directory.
  - It returns `(tdata, data)` and wraps `plot_signals("mag_probes", ...)`.
- Plot printing helpers:
  - `plot_scalar.py` supports `print=True` to print the plotted 1D y-data.
  - `plot_flux_average.py` supports `print=True` to print the plotted 1D y-data.
  - Printed numeric output should use fixed-width aligned columns with 8 values per row and general-format precision equivalent to `"{value:12.6g}"`.
- Plot auto-ylim conventions:
  - `plot_scalar.py`, `plot_flux_average.py`, `plot_signals.py`, and `plot_field_spectrum.py` use 10% padding on both sides of the plotted y-range when auto-setting `ylim`.
  - Add a horizontal line at `y=0` when the padded y-limits cross zero.
- `plot_scalar.py` axis conventions:
  - Public arguments use `xrange` and `yrange`.
  - If `xrange` is not set and `xlog` is false, the x-axis should start at `0`.
  - If `xrange` is set and `yrange` is not, derive y-limits from the visible x-window only.
  - Use a tolerance of `1e-9` when deciding whether data is effectively all positive or all negative.
- `flux_average.py` / `plot_flux_average.py` interface updates:
  - Public argument name is `timeslices`.
  - Accept `psi_norm`, `phi_norm`, and `rho` to select the x-axis coordinate.
  - In `plot_flux_average.py`, force the x-axis to `[0, 1]` on linear scale.
  - In `plot_flux_average.py`, use a tolerance of `0.001` for sign/near-zero decisions.
  - `plot_flux_average('q', timeslices)` should build flux coordinates using the requested `timeslices` value, not silently fall back to slice `0`.
  - For `eqsubtract=1` files, flux-coordinate construction for `plot_flux_average('q', timeslices)` should use total fields at the requested slice, which means `read_field` combines the perturbation at `timeslices` with equilibrium slice `-1`.
  - `flux_average.py` should preload `psi` for flux-coordinate construction with `equilibrium=False` and pass `slice=int(timeslices or 0)` into `flux_coordinates(...)`.
- `flux_average_field.py` correctness fix:
  - Do not transpose the `field_at_point(...)` output before flux-surface averaging.
- `field_spectrum.py` / `read_field_spectrum.py` / `plot_field_spectrum.py` conventions:
  - Default to reading fields with `complex=True` for spectrum calculations, with fallback to real data if the complex companion field is missing.
  - `read_field_spectrum.py` defaults to `pest=True` if none of `pest`, `boozer`, `hamada`, or `fast` is selected.
  - `field_spectrum.py` and `read_field_spectrum.py` support `m_val`; `read_field_spectrum.py` also accepts `m_vals` as an alias.
  - `field_spectrum.py` supports x-axis selection via `psi_norm`, `phi_norm`, and `rho`.
  - `plot_field_spectrum.py` should use Matplotlib default `axes.prop_cycle` colors.
  - If `m_val` is not provided to `plot_field_spectrum.py`, choose the 5 `m` values with the largest maximum absolute amplitudes and plot those.
  - When auto-selecting `m` values in `plot_field_spectrum.py`, ignore `NaN` amplitudes and rank using finite amplitudes across both positive and negative `m`.
  - For each `q_target` in `plot_field_spectrum.py`, draw vertical resonance lines for all profile crossings, not just one interpolated crossing.
- `read_scalar.py` / `plot_scalar.py` interface update:
  - Use public argument name `growth` instead of `growth_rate`.
  - Keep `growth_rate` only as a backward-compatible alias when needed.
- `read_hmn.py` / `plot_hmn.py` conventions:
  - Keep harmonics reading and plotting split across two files:
    - `read_hmn.py` for data loading / metadata assembly
    - `plot_hmn.py` for plotting only
  - `read_hmn.py` should return raw array data or metadata via `return_meta=True`.
  - `plot_hmn.py` should call `read_hmn(...)` rather than duplicating harmonics read logic.
  - Respect `cgs` and `mks` in both functions, and pass them through to `read_scalar("time", ...)`.
  - `plot_hmn.py` should use Matplotlib default `axes.prop_cycle` colors.
  - `plot_hmn.py` x-axis should start at `0` by default and accept explicit `xrange`.
  - When `me=True`, auto `ylim` should be computed from all components except `n=0`.
  - Auto `ylim` in `plot_hmn.py` should extend the upper bound by 10%.
  - When `growth=True`, add a horizontal line at `y=0`.
  - `plot_hmn.py` supports `print=<int>` to print the selected harmonic component `n` values with aligned formatting, 8 values per row.
- `contour_and_legend.py` level handling:
  - In `contour_and_legend_single(...)`, if explicit contour level values are provided, expand that existing level span by 1% and regenerate the same number of levels over the expanded span.
  - If contour levels are not provided, or only a level count is provided, build levels from the panel min/max and expand that span by 1%.
  - `contour_and_legend.py` supports `fill=True/False`:
    - `fill=True` keeps the filled contour behavior.
    - `fill=False` should draw contour lines only and skip the colorbar path.
- `plot_field.py` interface update:
  - Add `colorbar` argument; when `colorbar=False`, do not create a colorbar for the contour plot.
- LCFS helper behavior:
  - `get_lcfs.py` must filter kwargs before forwarding to `read_field(...)` and `read_lcfs(...)`, so plotting-only kwargs like `lines` do not break LCFS reads.
- `flux_average.py` / `flux_average_field.py` return-meta behavior:
  - `return_meta=True` should return structured Python objects, not bare tuples.
  - `flux_average.py` returns an object with `data`, `title`, `symbol`, `units`, and `fc`.
  - `flux_average_field.py` returns an object with `data`, `flux`, `nflux`, `area`, `volume`, and `r0`.
- `flux_coordinates.py` toroidal-flux normalization:
  - When building `phi_norm` / `rho`, use the last finite toroidal-flux value instead of blindly using `phi[-1]`.
  - Build `rho` from a clipped nonnegative normalized toroidal flux to avoid all-`NaN` output when the last toroidal-flux point is invalid.
  - When `slice >= 0`, `flux_coordinates.py` should read `psi`, `psi_r`, `psi_z`, and `I` with `equilibrium=False` so total-field reconstruction works correctly for `eqsubtract=1`.
  - When `slice < 0`, `flux_coordinates.py` should still use `equilibrium=True` for equilibrium-only reads.
- `read_field.py` equilibrium handling:
  - For `equilibrium=True` with `eqsubtract=1`, force reads to slice `-1` to avoid recursive total-field reconstruction loops.
  - `flux_coordinates.py` should pass its `slice` argument through to `lcfs(...)` so LCFS quantities are read from the matching time slice.
- `read_field.py` equilibrium handling:
  - For `equilibrium=True` with `eqsubtract=1`, force reads to slice `-1` to avoid recursive total-field reconstruction loops.
- `lcfs.py` / `read_lcfs.py` slice behavior:
  - `lcfs.py` should accept an explicit `slice` argument and pass it directly to `read_lcfs(...)`.
  - `read_lcfs.py` should respect negative slices in Python-style form, e.g. `slice=-1` means the last available slice.
  - `read_lcfs.py` should reproduce the IDL print behavior:
    - `slice time = ...`
    - `time step time: ...`
    - `time slice: ...`
- Legend helper behavior:
  - `plot_legend.py` should not hardcode legend font size; it should follow the active Matplotlib rc settings.
  - `plot_legend.py` currently uses `frameon=False`.
- `a2cc.f90` migration notes:
  - Python port lives in `m3dc1/a2cc.py` and related EQDSK-A parsing logic lives in a separate module.
  - Preserve Fortran-style comments and split helper modules when the original Fortran calls another file.
  - Map Fortran `write(*,...)` to Python `print(...)` on stdout.
  - Map Fortran `write(0,...)` to Python `print(..., file=sys.stderr)`.
