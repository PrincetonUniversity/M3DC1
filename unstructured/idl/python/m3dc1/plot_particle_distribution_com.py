from __future__ import annotations

from pathlib import Path

from .plot_particle_distribution import plot_particle_distribution


def plot_particle_distribution_com(
    timeslices: int = 0,
    *,
    filename: str | Path | None = None,
    field_filename: str | Path = "C1.h5",
    deltaf: bool = False,
    absolute_value: bool = True,
    sigma: int = 0,
    energy: float | None = None,
    energy_width: float = 1.0,
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
    energy; the default width is 1 keV.
    """
    return plot_particle_distribution(
        timeslices,
        filename=filename,
        field_filename=field_filename,
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
