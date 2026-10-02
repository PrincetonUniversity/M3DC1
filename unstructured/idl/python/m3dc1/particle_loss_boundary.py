from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, proton_mass
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter

from .field_at_point import field_at_point
from .read_field import read_field
from .read_lcfs import read_lcfs
from .read_parameter import read_parameter
from .read_particles import _field_values_at_markers


@dataclass
class ParticleLossBoundary:
    lambda_values: np.ndarray
    co_passing_outboard: np.ndarray
    co_passing_inboard: np.ndarray
    counter_passing_outboard: np.ndarray
    counter_passing_inboard: np.ndarray
    co_axis: np.ndarray
    counter_axis: np.ndarray
    pphi_values: np.ndarray
    lambda_upper: np.ndarray
    lambda_trapped_passing: np.ndarray
    outboard_loss_lambda_max: float
    inboard_loss_lambda_max: float
    axis_lambda_max: float
    psi_axis: float
    psi_edge: float

    @property
    def co_passing(self) -> np.ndarray:
        """Backward-compatible name for the original outboard co branch."""
        return self.co_passing_outboard

    @property
    def counter_passing(self) -> np.ndarray:
        """Backward-compatible name for the original inboard counter branch."""
        return self.counter_passing_inboard

    @property
    def co_loss_lambda_max(self) -> float:
        return self.outboard_loss_lambda_max

    @property
    def counter_loss_lambda_max(self) -> float:
        return self.inboard_loss_lambda_max


def _particle_mass_charge(sps: int, filename: str | Path) -> tuple[float, float]:
    if sps not in (1, 2):
        raise ValueError(f"sps must be 1 or 2, got {sps!r}.")

    if sps == 1:
        mass_number = float(read_parameter("ion_mass", filename=filename))
        charge_number = float(read_parameter("z_ion", filename=filename))
    else:
        mass_number = float(read_parameter("fast_ion_mass", filename=filename))
        charge_number = float(read_parameter("fast_ion_z", filename=filename))
        if mass_number <= 0.0:
            mass_number = float(read_parameter("ion_mass", filename=filename))
        if charge_number == 0.0:
            charge_number = float(read_parameter("z_ion", filename=filename))

    if mass_number <= 0.0 or charge_number == 0.0:
        raise ValueError(f"Valid particle mass and charge are required in {filename}.")
    return mass_number * proton_mass, charge_number * elementary_charge


def _midplane_point(
    *,
    filename: str | Path,
    timeslices: int,
    psi_edge: float,
    axis: np.ndarray,
    points: int,
    field_phi: float,
    outboard: bool,
) -> tuple[float, float]:
    psi = read_field(
        "psi",
        timeslices=int(timeslices),
        filename=filename,
        points=int(points),
        phi=float(field_phi),
        equilibrium=True,
        mks=True,
        return_meta=True,
    )
    axis_r, axis_z = np.asarray(axis, dtype=float).reshape(-1)[:2]
    line_points = max(4 * int(points), 200)
    radial_limit = float(np.max(psi.r) if outboard else np.min(psi.r))
    radius = np.linspace(axis_r, radial_limit, line_points)
    height = np.full(radius.shape, axis_z, dtype=float)
    psi_line = np.asarray(
        field_at_point(psi.data, psi.r, psi.z, radius, height),
        dtype=float,
    )
    valid = np.isfinite(psi_line)
    if psi.mask is not None:
        mask_line = np.asarray(
            field_at_point(psi.mask, psi.r, psi.z, radius, height),
            dtype=float,
        )
        valid &= np.isfinite(mask_line) & (mask_line < 0.5)
    if np.count_nonzero(valid) < 2:
        side = "outboard" if outboard else "inboard"
        raise ValueError(f"Could not sample the {side} magnetic midplane.")

    residual = psi_line - float(psi_edge)
    crossing = np.flatnonzero(
        valid[:-1]
        & valid[1:]
        & (residual[:-1] * residual[1:] <= 0.0)
    )
    if crossing.size:
        index = int(crossing[-1])
        denominator = residual[index + 1] - residual[index]
        fraction = 0.0 if np.isclose(denominator, 0.0) else -residual[index] / denominator
        boundary_r = radius[index] + fraction * (radius[index + 1] - radius[index])
    else:
        valid_indices = np.flatnonzero(valid)
        index = int(valid_indices[np.argmin(np.abs(residual[valid]))])
        boundary_r = radius[index]
    return float(boundary_r), float(axis_z)


def particle_loss_boundary(
    energy: float,
    *,
    lambda_values=None,
    pphi_values=None,
    timeslices: int = 0,
    field_filename: str | Path = "C1.h5",
    sps: int = 2,
    points: int = 200,
    field_points: int = 200,
    field_phi: float = 0.0,
) -> ParticleLossBoundary:
    """Calculate co- and counter-passing loss boundaries in COM space.

    ``energy`` is in keV and ``lambda_values`` contains
    ``Lambda = mu B0 / E``. The returned canonical momentum uses the same
    ``(P_phi-psi_axis)/(psi_axis-psi_edge)`` normalization as the COM plot.
    """
    energy_kev = float(energy)
    if not np.isfinite(energy_kev) or energy_kev <= 0.0:
        raise ValueError(f"energy must be finite and positive, got {energy!r}.")
    midplane_points = int(points)
    if midplane_points < 3:
        raise ValueError("points must be at least 3.")

    if lambda_values is None:
        lambdas = np.linspace(0.0, 1.0, 200)
    else:
        lambdas = np.asarray(lambda_values, dtype=float).reshape(-1)
    if lambdas.size < 2 or np.any(~np.isfinite(lambdas)):
        raise ValueError("lambda_values must contain at least two finite values.")
    if pphi_values is None:
        normalized_pphi = np.linspace(-1.5, 0.5, 200)
    else:
        normalized_pphi = np.asarray(pphi_values, dtype=float).reshape(-1)
    if normalized_pphi.size < 2 or np.any(~np.isfinite(normalized_pphi)):
        raise ValueError("pphi_values must contain at least two finite values.")

    lcfs = read_lcfs(
        filename=field_filename,
        slice=int(timeslices),
        mks=True,
        return_meta=True,
    )
    psi_axis = float(lcfs.flux0)
    psi_edge = float(lcfs.psilim)
    flux_span = psi_axis - psi_edge
    if not np.isfinite(flux_span) or np.isclose(flux_span, 0.0):
        raise ValueError("A finite nonzero psi_axis-psi_edge is required.")

    midplane = [
        _midplane_point(
            filename=field_filename,
            timeslices=int(timeslices),
            psi_edge=psi_edge,
            axis=lcfs.axis,
            points=midplane_points,
            field_phi=float(field_phi),
            outboard=outboard,
        )
        for outboard in (True, False)
    ]
    axis_r, axis_z = np.asarray(lcfs.axis, dtype=float).reshape(-1)[:2]
    outboard_radius = np.linspace(axis_r, midplane[0][0], midplane_points)
    inboard_radius = np.linspace(axis_r, midplane[1][0], midplane_points)
    evaluation_r = np.concatenate([outboard_radius, inboard_radius[1:]])
    evaluation_z = np.full(evaluation_r.shape, axis_z, dtype=float)

    fields = _field_values_at_markers(
        evaluation_r,
        evaluation_z,
        filename=field_filename,
        timeslices=int(timeslices),
        field_points=int(field_points),
        field_phi=float(field_phi),
        equilibrium=True,
        include_psi=True,
    )
    bmag = np.asarray(fields["bmag"], dtype=float)
    current = np.asarray(fields["I"], dtype=float)
    psi_midplane = np.asarray(fields["psi"], dtype=float)
    outboard_index = midplane_points - 1
    inboard_index = 2 * midplane_points - 2
    key_indices = np.asarray([0, outboard_index, inboard_index])
    valid_midplane = (
        np.isfinite(bmag[key_indices])
        & (bmag[key_indices] > 0.0)
        & np.isfinite(current[key_indices])
    )
    if not np.all(valid_midplane):
        raise ValueError("Could not evaluate I/B at the LCFS midplanes and magnetic axis.")

    bzero = float(read_parameter("bzero", filename=field_filename))
    b0_norm = float(read_parameter("b0_norm", filename=field_filename))
    reference_b = abs(bzero) * b0_norm / 1.0e4
    if not np.isfinite(reference_b) or reference_b <= 0.0:
        raise ValueError("A finite positive reference magnetic field is required.")

    mass, charge = _particle_mass_charge(int(sps), field_filename)
    energy_joule = energy_kev * 1.0e3 * elementary_charge
    def _branch(sign: float, index: int, flux: float) -> np.ndarray:
        parallel_fraction = 1.0 - lambdas * bmag[index] / reference_b
        accessible = parallel_fraction >= 0.0
        vparallel = np.sqrt(
            (2.0 * energy_joule / mass) * np.clip(parallel_fraction, 0.0, None)
        )
        orbit_shift = (mass / charge) * vparallel * current[index] / bmag[index]
        normalized = (flux + sign * orbit_shift - psi_axis) / flux_span
        return np.where(accessible, normalized, np.nan)

    outboard_indices = np.arange(midplane_points)
    inboard_indices = np.concatenate(
        [np.asarray([0]), np.arange(midplane_points, 2 * midplane_points - 1)]
    )

    def _stagnation_limit(
        profile_indices: np.ndarray,
        edge_index: int,
        *,
        increases_with_pphi: bool,
    ) -> np.ndarray:
        profile_b = bmag[profile_indices]
        profile_psi = psi_midplane[profile_indices]
        valid_profile = (
            np.isfinite(profile_b)
            & (profile_b > 0.0)
            & np.isfinite(profile_psi)
        )
        limit = np.full(normalized_pphi.shape, np.nan, dtype=float)
        valid_indices = np.flatnonzero(valid_profile)
        if valid_indices.size < 3:
            return limit

        smooth_window = min(31, valid_indices.size)
        if smooth_window % 2 == 0:
            smooth_window -= 1

        def _smooth_profile(values: np.ndarray) -> np.ndarray:
            selected = np.asarray(values[valid_indices], dtype=float)
            if smooth_window < 5:
                return selected
            return savgol_filter(
                selected,
                window_length=smooth_window,
                polyorder=min(3, smooth_window - 1),
                mode="interp",
            )

        smooth_b = _smooth_profile(profile_b)
        smooth_psi = _smooth_profile(profile_psi)
        smooth_b[0] = bmag[0]
        smooth_b[-1] = bmag[edge_index]
        smooth_psi[0] = psi_axis
        smooth_psi[-1] = psi_edge
        profile_pphi = (smooth_psi - psi_axis) / flux_span
        profile_lambda = reference_b / smooth_b
        interior = (
            np.isfinite(profile_pphi)
            & np.isfinite(profile_lambda)
            & (profile_pphi > -1.0)
            & (profile_pphi < 0.0)
        )
        pphi_nodes = np.concatenate([[-1.0], profile_pphi[interior], [0.0]])
        lambda_nodes = np.concatenate(
            [
                [reference_b / bmag[edge_index]],
                profile_lambda[interior],
                [reference_b / bmag[0]],
            ]
        )
        order = np.argsort(pphi_nodes)
        pphi_nodes = pphi_nodes[order]
        lambda_nodes = lambda_nodes[order]
        pphi_nodes, unique_indices = np.unique(pphi_nodes, return_index=True)
        lambda_nodes = lambda_nodes[unique_indices]

        if increases_with_pphi:
            lambda_nodes = np.maximum.accumulate(lambda_nodes)
        else:
            lambda_nodes = np.minimum.accumulate(lambda_nodes)
        stagnation_limit = PchipInterpolator(pphi_nodes, lambda_nodes)
        inside = (normalized_pphi >= -1.0) & (normalized_pphi <= 0.0)
        limit[inside] = stagnation_limit(normalized_pphi[inside])
        return limit

    # On the low-field side, this is the maximum Lambda accessible at a
    # v_parallel=0 stagnation point.
    lambda_upper = _stagnation_limit(
        outboard_indices,
        outboard_index,
        increases_with_pphi=False,
    )
    # On the high-field side, this is where the second midplane crossing has
    # v_parallel=0 and the orbit changes between passing and trapped.
    lambda_trapped_passing = _stagnation_limit(
        inboard_indices,
        inboard_index,
        increases_with_pphi=True,
    )

    return ParticleLossBoundary(
        lambda_values=lambdas,
        co_passing_outboard=_branch(1.0, outboard_index, psi_edge),
        co_passing_inboard=_branch(1.0, inboard_index, psi_edge),
        counter_passing_outboard=_branch(-1.0, outboard_index, psi_edge),
        counter_passing_inboard=_branch(-1.0, inboard_index, psi_edge),
        co_axis=_branch(1.0, 0, psi_axis),
        counter_axis=_branch(-1.0, 0, psi_axis),
        pphi_values=normalized_pphi,
        lambda_upper=lambda_upper,
        lambda_trapped_passing=lambda_trapped_passing,
        outboard_loss_lambda_max=reference_b / bmag[outboard_index],
        inboard_loss_lambda_max=reference_b / bmag[inboard_index],
        axis_lambda_max=reference_b / bmag[0],
        psi_axis=psi_axis,
        psi_edge=psi_edge,
    )
