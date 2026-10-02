from __future__ import annotations

from pathlib import Path

import numpy as np

from .particle_loss_boundary import particle_loss_boundary
from .plot_particle_distribution import plot_particle_distribution


def plot_particle_distribution_com(
    timeslices: int = 0,
    *,
    filename: str | Path | None = None,
    field_filename: str | Path = "C1.h5",
    sps: int | None = None,
    deltaf: bool = False,
    absolute_value: bool = True,
    sigma: int = 0,
    energy: float | None = None,
    energy_width: float = 1.0,
    loss_boundary: bool = True,
    boundary_energy: float | None = None,
    boundary_sps: int | None = None,
    boundary_points: int = 200,
    boundary_color: str = "black",
    boundary_linewidth: float = 1.0,
    boundary_legend: bool = False,
    max_particles: int | None = 10_000,
    field_points: int = 200,
    field_phi: float = 0.0,
    points: int = 100,
    levels=100,
    xrange=None,
    yrange=None,
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
    """Plot all-species KDE in normalized canonical momentum and pitch.

    By default, uniformly sample at most 10,000 markers. Set
    ``max_particles=None`` to use every marker row. Canonical momentum is
    normalized as ``(P_phi-psi0)/(psi0-psi_edge)``.
    ``sigma=1`` plots ``v_parallel>0``, ``sigma=-1`` plots
    ``v_parallel<0``, and the default ``sigma=0`` combines both populations.
    ``energy`` selects particles within ``energy_width`` keV of the requested
    energy; the default width is 1 keV. When an energy is selected, dashed
    co- and counter-passing confined/loss boundaries are overlaid by default.
    Each passing direction includes separate outboard- and inboard-midplane
    LCFS-touch curves because an orbit is lost if it reaches either side.
    ``boundary_energy`` can instead specify the boundary energy independently.
    Use ``sps=1`` for thermal ions, ``sps=2`` for fast ions, or ``sps=None``
    for all particles. Unless explicitly set, ``boundary_sps`` follows ``sps``;
    when plotting all particles, the boundary defaults to fast-ion parameters.
    """
    figure, axis = plot_particle_distribution(
        timeslices,
        filename=filename,
        field_filename=field_filename,
        sps=sps,
        deltaf=deltaf,
        absolute_value=absolute_value,
        coordinates="com",
        sigma=sigma,
        energy=energy,
        energy_width=energy_width,
        max_particles=max_particles,
        field_points=field_points,
        field_phi=field_phi,
        points=points,
        levels=levels,
        xrange=xrange,
        yrange=yrange,
        xscale=xscale,
        yscale=yscale,
        bandwidth=bandwidth,
        cmap=cmap,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        colorbar=colorbar,
        colorbar_label=colorbar_label,
        overplot=overplot,
        outfile=outfile,
    )
    loss_energy = energy if boundary_energy is None else boundary_energy
    if loss_boundary and loss_energy is not None:
        yscale_value = float(yscale)
        if not np.isfinite(yscale_value) or np.isclose(yscale_value, 0.0):
            raise ValueError("yscale must be finite and nonzero for loss boundaries.")
        xscale_value = float(xscale)
        if not np.isfinite(xscale_value) or np.isclose(xscale_value, 0.0):
            raise ValueError("xscale must be finite and nonzero for loss boundaries.")
        plot_xlim = axis.get_xlim()
        plot_ylim = axis.get_ylim()
        ylimits = np.asarray(plot_ylim, dtype=float) / yscale_value
        lambdas = np.linspace(float(np.min(ylimits)), float(np.max(ylimits)), boundary_points)
        xlimits = np.asarray(plot_xlim, dtype=float) / xscale_value
        pphi_values = np.linspace(float(np.min(xlimits)), float(np.max(xlimits)), boundary_points)
        if boundary_sps is None:
            boundary_species = 2 if sps is None else int(sps)
        else:
            boundary_species = int(boundary_sps)
        boundary = particle_loss_boundary(
            loss_energy,
            lambda_values=lambdas,
            pphi_values=pphi_values,
            timeslices=timeslices,
            field_filename=field_filename,
            sps=boundary_species,
            points=boundary_points,
            field_points=field_points,
            field_phi=field_phi,
        )
        line_options = {
            "color": boundary_color,
            "linewidth": float(boundary_linewidth),
        }

        def _passing_curve(values, lambda_max, endpoint):
            values = np.asarray(values, dtype=float)
            valid = np.isfinite(values) & (boundary.lambda_values < lambda_max)
            curve_x = np.concatenate([values[valid], [endpoint]]) * xscale_value
            curve_y = np.concatenate(
                [boundary.lambda_values[valid], [lambda_max]]
            ) * yscale_value
            return curve_x, curve_y

        co_outboard_x, co_outboard_y = _passing_curve(
            boundary.co_passing_outboard,
            boundary.outboard_loss_lambda_max,
            -1.0,
        )
        co_inboard_x, co_inboard_y = _passing_curve(
            boundary.co_passing_inboard,
            boundary.inboard_loss_lambda_max,
            -1.0,
        )
        counter_outboard_x, counter_outboard_y = _passing_curve(
            boundary.counter_passing_outboard,
            boundary.outboard_loss_lambda_max,
            -1.0,
        )
        counter_inboard_x, counter_inboard_y = _passing_curve(
            boundary.counter_passing_inboard,
            boundary.inboard_loss_lambda_max,
            -1.0,
        )
        co_axis_x, co_axis_y = _passing_curve(
            boundary.co_axis,
            boundary.axis_lambda_max,
            0.0,
        )
        counter_axis_x, counter_axis_y = _passing_curve(
            boundary.counter_axis,
            boundary.axis_lambda_max,
            0.0,
        )

        first_boundary_line = len(axis.lines)
        axis.plot(
            co_outboard_x,
            co_outboard_y,
            label="co-passing outboard touch",
            linestyle=(0, (6, 2)),
            **line_options,
        )
        axis.plot(
            co_inboard_x,
            co_inboard_y,
            label="co-passing inboard touch",
            linestyle=(0, (6, 2, 1, 2)),
            **line_options,
        )
        axis.plot(
            counter_outboard_x,
            counter_outboard_y,
            label="counter-passing outboard touch",
            linestyle=(0, (2, 2)),
            **line_options,
        )
        axis.plot(
            counter_inboard_x,
            counter_inboard_y,
            label="counter-passing inboard touch",
            linestyle=(0, (2, 2, 1, 2)),
            **line_options,
        )
        axis.plot(
            np.full(2, -1.0 * xscale_value),
            np.asarray(
                [boundary.inboard_loss_lambda_max, boundary.outboard_loss_lambda_max]
            )
            * yscale_value,
            label="trapped-particle loss boundary",
            linestyle=(0, (4, 2)),
            **line_options,
        )
        axis.plot(
            co_axis_x,
            co_axis_y,
            label="co-passing magnetic-axis line",
            linestyle=(0, (6, 2, 1, 2)),
            **line_options,
        )
        axis.plot(
            counter_axis_x,
            counter_axis_y,
            label="counter-passing magnetic-axis line",
            linestyle=(0, (2, 2, 1, 2)),
            **line_options,
        )
        trapped_valid = (
            np.isfinite(boundary.lambda_trapped_passing)
            & (boundary.pphi_values > -1.0)
            & (boundary.pphi_values < 0.0)
        )
        trapped_x = np.concatenate(
            [[-1.0], boundary.pphi_values[trapped_valid], [0.0]]
        ) * xscale_value
        trapped_y = np.concatenate(
            [
                [boundary.inboard_loss_lambda_max],
                boundary.lambda_trapped_passing[trapped_valid],
                [boundary.axis_lambda_max],
            ]
        ) * yscale_value
        axis.plot(
            trapped_x,
            trapped_y,
            label="trapped-passing boundary",
            linestyle=(0, (5, 2, 1, 2, 1, 2)),
            **line_options,
        )
        upper_valid = (
            np.isfinite(boundary.lambda_upper)
            & (boundary.pphi_values > -1.0)
            & (boundary.pphi_values < 0.0)
        )
        upper_x = np.concatenate(
            [[-1.0], boundary.pphi_values[upper_valid], [0.0]]
        ) * xscale_value
        upper_y = np.concatenate(
            [
                [boundary.outboard_loss_lambda_max],
                boundary.lambda_upper[upper_valid],
                [boundary.axis_lambda_max],
            ]
        ) * yscale_value
        axis.plot(
            upper_x,
            upper_y,
            label=r"upper $\Lambda$ limit",
            linestyle=(0, (8, 2, 2, 2)),
            **line_options,
        )
        line_xvalues = []
        line_yvalues = []
        for line in axis.lines[first_boundary_line:]:
            line_x = np.asarray(line.get_xdata(), dtype=float)
            line_y = np.asarray(line.get_ydata(), dtype=float)
            finite_line = np.isfinite(line_x) & np.isfinite(line_y)
            line_xvalues.append(line_x[finite_line])
            line_yvalues.append(line_y[finite_line])
        finite_x = np.concatenate([values for values in line_xvalues if values.size])
        finite_y = np.concatenate([values for values in line_yvalues if values.size])

        xspan = plot_xlim[1] - plot_xlim[0]
        yspan = plot_ylim[1] - plot_ylim[0]
        expanded_xlim = (
            min(plot_xlim[0], float(np.min(finite_x)) - 0.01 * xspan),
            max(plot_xlim[1], float(np.max(finite_x)) + 0.01 * xspan),
        )
        expanded_ylim = (
            min(plot_ylim[0], float(np.min(finite_y)) - 0.01 * yspan),
            max(plot_ylim[1], float(np.max(finite_y)) + 0.01 * yspan),
        )
        axis.set_xlim(expanded_xlim)
        axis.set_ylim(expanded_ylim)
        if boundary_legend:
            axis.legend(frameon=False)
        if not overplot:
            figure.tight_layout()
        if outfile is not None:
            figure.savefig(str(outfile))
    return figure, axis
