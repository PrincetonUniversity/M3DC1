from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .field_at_point import field_at_point
from .flux_coordinates import flux_coordinates
from .plot_poincare import run_trace
from .read_field import read_field


def _read_trace_q(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    try:
        data = np.loadtxt(path, ndmin=2)
    except ValueError as error:
        raise ValueError(f"{path} does not contain valid numeric trace q data.") from error

    if data.ndim != 2 or data.shape[0] == 0 or data.shape[1] < 4:
        raise ValueError(f"{path} has shape {data.shape}; four columns are required.")

    data = np.asarray(data[:, :4], dtype=float)
    finite = np.all(np.isfinite(data), axis=1)
    if not np.any(finite):
        raise ValueError(f"{path} does not contain any finite trace q data.")

    radius, qmean, qmin, qmax = data[finite].T
    return radius, qmean, qmin, qmax, np.flatnonzero(finite)


def _flux_xcoordinate(
    coordinate: str,
    *,
    filename: str | Path,
    seed_indices: np.ndarray,
    dR: float,
    dZ: float,
    dR0: float | None,
    dZ0: float | None,
    points: int,
) -> tuple[np.ndarray, str]:
    fc = flux_coordinates(filename=filename, points=points, slice=-1)
    offset_r = float(dR) if dR0 is None or float(dR0) == 0.0 else float(dR0)
    offset_z = 0.0 if dZ0 is None else float(dZ0)
    indices = np.asarray(seed_indices, dtype=float)
    seed_r = float(fc.r0) + offset_r + indices * float(dR)
    seed_z = float(fc.z0) + offset_z + indices * float(dZ)

    psi = read_field(
        "psi",
        filename=filename,
        timeslices=-1,
        equilibrium=True,
        points=points,
        return_meta=True,
    )
    psi_data = np.asarray(psi.data)
    if psi_data.ndim == 3:
        psi_data = psi_data[0]
    psi_seed = np.asarray(
        field_at_point(psi_data, psi.r, psi.z, seed_r, seed_z),
        dtype=float,
    )

    psi_grid = np.asarray(fc.psi, dtype=float).reshape(-1)
    psin_grid = np.asarray(fc.psi_norm, dtype=float).reshape(-1)
    finite = np.isfinite(psi_grid) & np.isfinite(psin_grid)
    if np.count_nonzero(finite) < 2:
        raise ValueError("Could not determine normalized poloidal flux for trace seeds.")
    slope, intercept = np.polyfit(psin_grid[finite], psi_grid[finite], 1)
    if not np.isfinite(slope) or abs(slope) <= np.finfo(float).tiny:
        raise ValueError("The poloidal-flux span is zero or invalid.")
    psi_norm_seed = (psi_seed - intercept) / slope

    if coordinate == "psi_norm":
        return psi_norm_seed, r"$\psi_N$"

    phi_grid = np.asarray(fc.phi_norm, dtype=float).reshape(-1)
    finite = np.isfinite(psin_grid) & np.isfinite(phi_grid)
    xp = np.r_[0.0, psin_grid[finite], 1.0]
    fp = np.r_[0.0, phi_grid[finite], 1.0]
    order = np.argsort(xp)
    xp = xp[order]
    fp = fp[order]
    xp, unique = np.unique(xp, return_index=True)
    fp = fp[unique]
    phi_norm_seed = np.interp(psi_norm_seed, xp, fp, left=np.nan, right=np.nan)
    if coordinate == "phi_norm":
        return phi_norm_seed, r"$\Phi_N$"
    return np.sqrt(np.clip(phi_norm_seed, 0.0, None)), r"$\rho=\sqrt{\Phi_N}$"


def plot_trace_q(
    *,
    qfile: str | Path = "q.out",
    directory: str | Path = ".",
    filename: str | Path = "C1.h5",
    slice: int = 0,
    run_trace_if_missing: bool = True,
    trace_executable: str | Path | None = None,
    trace_dR: float = 0.1,
    trace_dZ: float = 0.1,
    trace_dR0: float | None = None,
    trace_dZ0: float | None = None,
    trace_surfaces: int = 51,
    trace_transits: int = 100,
    trace_steps_per_transit: int = 100,
    trace_angle: float = 0.0,
    trace_phi0: float = 0.0,
    trace_field_scale: float = 1.0,
    trace_field_phase: float = 0.0,
    trace_reverse: bool = False,
    trace_extra_args: Sequence[str | int | float] | None = None,
    psi: bool = False,
    phi: bool = False,
    psi_norm: bool = False,
    phi_norm: bool = False,
    rho: bool = False,
    flux_points: int = 200,
    bounds: bool = False,
    fill_bounds: bool = False,
    color=None,
    linestyle: str = "-",
    linewidth: float | None = None,
    marker=None,
    label: str = "mean q",
    xrange=None,
    yrange=None,
    xscale: float = 1.0,
    yscale: float = 1.0,
    overplot: bool = False,
    title: str | None = None,
    outfile: str | Path | None = None,
):
    """Plot the safety-factor profile written by fusion-io ``trace``."""
    output_directory = Path(directory).expanduser().resolve()
    path = Path(qfile).expanduser()
    if not path.is_absolute():
        path = output_directory / path

    if not path.is_file() and run_trace_if_missing:
        if Path(qfile).name != "q.out":
            raise FileNotFoundError(
                "Automatic trace generation writes q.out; use qfile='q.out' or generate the custom file first."
            )
        run_trace(
            filename=filename,
            directory=output_directory,
            timeslice=slice,
            trace_executable=trace_executable,
            dR=trace_dR,
            dZ=trace_dZ,
            dR0=trace_dR0,
            dZ0=trace_dZ0,
            surfaces=trace_surfaces,
            transits=trace_transits,
            steps_per_transit=trace_steps_per_transit,
            angle=trace_angle,
            phi0=trace_phi0,
            field_scale=trace_field_scale,
            field_phase=trace_field_phase,
            reverse=trace_reverse,
            trace_extra_args=trace_extra_args,
        )

    if not path.is_file():
        raise FileNotFoundError(f"Trace q-profile file does not exist: {path}")

    radius, qmean, qmin, qmax, seed_indices = _read_trace_q(path)
    use_psi = bool(psi or psi_norm)
    use_phi = bool(phi or phi_norm)
    ncoordinates = int(use_psi) + int(use_phi) + int(rho)
    if ncoordinates > 1:
        raise TypeError("plot_trace_q() accepts only one of psi, phi, or rho.")
    if use_psi:
        xdata, xlabel = _flux_xcoordinate(
            "psi_norm",
            filename=filename,
            seed_indices=seed_indices,
            dR=trace_dR,
            dZ=trace_dZ,
            dR0=trace_dR0,
            dZ0=trace_dZ0,
            points=flux_points,
        )
    elif use_phi:
        xdata, xlabel = _flux_xcoordinate(
            "phi_norm",
            filename=filename,
            seed_indices=seed_indices,
            dR=trace_dR,
            dZ=trace_dZ,
            dR0=trace_dR0,
            dZ0=trace_dZ0,
            points=flux_points,
        )
    elif rho:
        xdata, xlabel = _flux_xcoordinate(
            "rho",
            filename=filename,
            seed_indices=seed_indices,
            dR=trace_dR,
            dZ=trace_dZ,
            dR0=trace_dR0,
            dZ0=trace_dZ0,
            points=flux_points,
        )
    else:
        xdata = radius
        xlabel = "Initial distance from magnetic axis (m)"

    finite = np.isfinite(xdata)
    if not np.any(finite):
        raise ValueError("No trace seed points have a finite requested flux coordinate.")
    xdata = np.asarray(xdata[finite], dtype=float)
    qmean = qmean[finite]
    qmin = qmin[finite]
    qmax = qmax[finite]
    order = np.argsort(xdata)
    xdata = xdata[order] * float(xscale)
    qmean = qmean[order]
    qmin = qmin[order]
    qmax = qmax[order]
    qmean = qmean * float(yscale)
    qmin = qmin * float(yscale)
    qmax = qmax * float(yscale)

    if outfile is not None:
        np.savetxt(outfile, np.column_stack([xdata, qmean, qmin, qmax]), fmt="%16.6e")

    if overplot:
        ax = plt.gca()
        fig = ax.figure
    else:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("q")

    (mean_line,) = ax.plot(
        xdata,
        qmean,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        marker=marker,
        label=label,
    )
    line_color = mean_line.get_color()
    if bounds:
        ax.plot(xdata, qmin, color=line_color, linestyle="--", linewidth=linewidth, label="q min")
        ax.plot(xdata, qmax, color=line_color, linestyle=":", linewidth=linewidth, label="q max")
    if fill_bounds:
        ax.fill_between(xdata, qmin, qmax, color=line_color, alpha=0.2, linewidth=0)

    if title is None and not overplot:
        title = "Trace safety factor"
    if title:
        ax.set_title(str(title))
    if xrange is not None:
        ax.set_xlim(xrange)
    elif xdata.size:
        ax.set_xlim(left=0.0)
    if yrange is not None:
        ax.set_ylim(yrange)
    if bounds or label:
        ax.legend(frameon=False)
    if not overplot:
        fig.tight_layout()
    return fig, ax
