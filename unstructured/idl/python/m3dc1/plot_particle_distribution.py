from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .read_lcfs import read_lcfs
from .read_parameter import read_parameter
from .read_particles import read_particles


_XI_YRANGE = (-1.0, 1.0)


def _orient_xi(xi: np.ndarray, field_filename: str | Path) -> np.ndarray:
    bzsign = float(read_parameter("bzsign", filename=field_filename))
    if np.isclose(bzsign, -1.0):
        return -xi
    return xi


def _normalize_pphi(
    pphi: np.ndarray,
    field_filename: str | Path,
    timeslices: int,
) -> np.ndarray:
    lcfs = read_lcfs(
        filename=field_filename,
        slice=int(timeslices),
        mks=True,
        return_meta=True,
    )
    edge_flux = float(lcfs.psilim)
    axis_flux = float(lcfs.flux0)
    flux_span = axis_flux - edge_flux
    if not np.isfinite(flux_span) or np.isclose(flux_span, 0.0):
        raise ValueError(
            "A finite nonzero psi0-psi_edge is required to normalize P_phi."
        )
    return (pphi - axis_flux) / flux_span


def _plot_range(values: np.ndarray, requested, name: str) -> tuple[float, float]:
    if requested is None:
        lower = float(np.min(values))
        upper = float(np.max(values))
    else:
        bounds = np.asarray(requested, dtype=float).reshape(-1)
        if bounds.size != 2:
            raise ValueError(f"{name} must contain exactly two values.")
        lower, upper = map(float, bounds)
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        raise ValueError(f"{name} must be a finite increasing range, got ({lower}, {upper}).")
    return lower, upper


def _particle_kde(
    values: np.ndarray,
    positions: np.ndarray,
    weights: np.ndarray,
    *,
    deltaf: bool,
    absolute_value: bool,
    bandwidth,
) -> np.ndarray:
    if not deltaf:
        kernel = stats.gaussian_kde(values, bw_method=bandwidth)
        return np.asarray(kernel(positions))

    absolute_weights = np.abs(weights)
    if not np.any(absolute_weights > 0.0):
        raise ValueError("All selected delta-f particle weights are zero.")

    kernel = stats.gaussian_kde(
        values,
        weights=absolute_weights,
        bw_method=bandwidth,
    )
    if not absolute_value:
        # SciPy disallows signed input weights. Use |w| for the covariance,
        # then evaluate the same Gaussian kernels with signed weights.
        kernel._weights = weights / np.sum(absolute_weights)
    return np.asarray(kernel(positions))


def _coordinate_labels(
    coordinate_mode: str,
    xlabel: str | None,
    ylabel: str | None,
) -> tuple[str, str]:
    if xlabel is None:
        xlabel = {
            "energy_xi": "Energy (keV)",
            "momentum": r"$p_{\parallel}$ (ion_mass m/s)",
            "com": r"$(P_{\phi}-\psi_0)/(\psi_0-\psi_{edge})$",
        }[coordinate_mode]
    if ylabel is None:
        ylabel = {
            "energy_xi": r"$\xi=v_{\parallel}/v$",
            "momentum": r"$p_{\perp}$ (ion_mass m/s)",
            "com": r"$\Lambda=\mu B_0/E$",
        }[coordinate_mode]
    return xlabel, ylabel


def plot_particle_distribution(
    timeslices: int = 0,
    *,
    filename: str | Path | None = None,
    field_filename: str | Path = "C1.h5",
    sps: int | None = None,
    minor_radius: float | None = None,
    minor_radius_width: float = 0.05,
    deltaf: bool = False,
    absolute_value: bool = True,
    momentum: bool = False,
    coordinates: str | None = None,
    sigma: int = 0,
    energy: float | None = None,
    energy_width: float = 1.0,
    max_particles: int | None = 10_000,
    field_points: int = 200,
    field_phi: float = 0.0,
    points: int = 100,
    levels=100,
    xrange=None,
    yrange=_XI_YRANGE,
    xscale: float = 1.0,
    yscale: float = 1.0,
    bandwidth=None,
    cmap=None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    colorbar: bool = False,
    colorbar_label: str | None = None,
    overplot: bool = False,
    outfile: str | Path | None = None,
):
    """Plot a full-f or delta-f 2D KDE from M3D-C1 marker information.

    ``timeslices`` selects ``ions_NNNN.h5`` beside ``field_filename``. Energy
    and ``xi = v_parallel / v`` are calculated using the magnetic-field
    magnitude at each marker location. ``field_phi`` selects the single field
    plane and is in degrees when ``itor=1``. Use ``sps=1`` for thermal ions,
    ``sps=2`` for fast ions, or ``None`` for all particles.
    If the field file parameter ``bzsign`` is -1, plotted xi is sign-reversed.
    ``minor_radius`` selects a normalized ``sqrt(psi_norm)`` shell.
    Set ``deltaf=True`` to apply particle weights. The default ``deltaf=False``
    plots the unweighted full marker distribution.
    For delta-f plots, ``absolute_value=True`` uses the absolute particle
    weights; set it to ``False`` to retain their signs.
    Set ``momentum=True`` to plot ``ppar`` against ``pperp`` instead of energy
    against xi. Momentum values are in ion_mass m/s before applying axis scales.
    ``coordinates`` may be ``"energy_xi"``, ``"momentum"``, or ``"com"``;
    the latter plots canonical momentum ``pphi`` against ``muB0overE``.
    ``sigma=1`` selects positive parallel velocity, ``sigma=-1`` selects
    negative parallel velocity, and the default ``sigma=0`` combines both.
    ``energy`` selects markers within ``energy_width`` keV of that energy.
    """
    sigma_value = int(sigma)
    if sigma_value not in (-1, 0, 1) or sigma_value != sigma:
        raise ValueError(f"sigma must be -1, 0, or 1, got {sigma!r}.")

    coordinate_mode = "momentum" if momentum else "energy_xi"
    if coordinates is not None:
        requested_mode = str(coordinates).strip().lower()
        if requested_mode not in {"energy_xi", "momentum", "com"}:
            raise ValueError(
                "coordinates must be 'energy_xi', 'momentum', 'com', or None, "
                f"got {coordinates!r}."
            )
        if momentum and requested_mode != "momentum":
            raise ValueError("momentum=True cannot be combined with another coordinates mode.")
        coordinate_mode = requested_mode

    physical_columns = {
        "energy_xi": ["energy", "xi"],
        "momentum": ["ppar", "pperp"],
        "com": ["pphi", "muB0overE"],
    }[coordinate_mode]
    read_columns = [*physical_columns]
    vparallel_column = None
    if sigma_value != 0:
        vparallel_column = len(read_columns)
        read_columns.append("v_parallel")
    result = read_particles(
        timeslices,
        filename=filename,
        columns=[*read_columns, "weight"],
        field_filename=field_filename,
        sps=sps,
        minor_radius=minor_radius,
        minor_radius_width=minor_radius_width,
        energy=energy,
        energy_width=energy_width,
        max_particles=max_particles,
        field_points=field_points,
        field_phi=field_phi,
        return_meta=True,
    )
    particle_data = np.asarray(result.data, dtype=float)
    xvalues = particle_data[:, 0]
    yvalues = particle_data[:, 1]
    vparallel = particle_data[:, vparallel_column] if vparallel_column is not None else None
    raw_weight = particle_data[:, -1]
    if coordinate_mode == "energy_xi":
        yvalues = _orient_xi(yvalues, field_filename)
    elif coordinate_mode == "com":
        xvalues = _normalize_pphi(xvalues, field_filename, timeslices)
    xdata = xvalues * float(xscale)
    ydata = yvalues * float(yscale)
    weights = raw_weight if deltaf else np.ones_like(raw_weight)
    finite = np.isfinite(xdata) & np.isfinite(ydata) & np.isfinite(weights)
    if vparallel is not None:
        finite &= np.isfinite(vparallel)
    range_xdata = xdata[finite]
    range_ydata = ydata[finite]
    if sigma_value > 0:
        finite &= vparallel > 0.0
    elif sigma_value < 0:
        finite &= vparallel < 0.0
    xdata = xdata[finite]
    ydata = ydata[finite]
    weights = weights[finite]
    if xdata.size < 3:
        raise ValueError("At least three finite particle markers are required for a 2D KDE.")
    if np.ptp(xdata) == 0.0 or np.ptp(ydata) == 0.0:
        raise ValueError("Particle x and y values must each span a nonzero range for a 2D KDE.")

    auto_physical_range = coordinate_mode != "energy_xi"
    requested_xrange = xrange
    requested_yrange = None if auto_physical_range and yrange is _XI_YRANGE else yrange
    if coordinate_mode == "energy_xi" and requested_xrange is None:
        xbounds = _plot_range(
            range_xdata,
            (0.0, float(np.max(range_xdata))),
            "xrange",
        )
    else:
        xbounds = _plot_range(range_xdata, requested_xrange, "xrange")
    if auto_physical_range and requested_yrange is None:
        ybounds = _plot_range(
            range_ydata,
            (0.0, float(np.max(range_ydata))),
            "yrange",
        )
    else:
        ybounds = _plot_range(range_ydata, requested_yrange, "yrange")
    grid_points = int(points)
    if grid_points < 2:
        raise ValueError("points must be at least 2.")

    xgrid, ygrid = np.mgrid[
        xbounds[0] : xbounds[1] : complex(grid_points),
        ybounds[0] : ybounds[1] : complex(grid_points),
    ]
    positions = np.vstack([xgrid.ravel(), ygrid.ravel()])

    xlabel, ylabel = _coordinate_labels(coordinate_mode, xlabel, ylabel)
    try:
        density = _particle_kde(
            np.vstack([xdata, ydata]),
            positions,
            weights=weights,
            deltaf=deltaf,
            absolute_value=absolute_value,
            bandwidth=bandwidth,
        )
    except (ValueError, np.linalg.LinAlgError) as exc:
        raise ValueError(f"Could not construct the particle KDE: {exc}") from exc
    density = density.reshape(xgrid.shape)

    if overplot:
        axis = plt.gca()
        figure = axis.figure
    else:
        figure, axis = plt.subplots(figsize=(6, 4.5))
    contour = axis.contourf(xgrid, ygrid, density, levels=levels, cmap=cmap)
    for collection in contour.collections:
        collection.set_edgecolor("face")

    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_xlim(xbounds)
    axis.set_ylim(ybounds)
    if coordinate_mode == "momentum":
        axis.axvline(0.0, color="black", linestyle="--", linewidth=0.8)
        axis.axhline(0.0, color="black", linestyle="--", linewidth=0.8)
    if title is not None:
        axis.set_title(str(title))
    if colorbar:
        figure.colorbar(contour, ax=axis, label=colorbar_label)
    if not overplot:
        figure.tight_layout()
    if outfile is not None:
        figure.savefig(str(outfile))
    return figure, axis
