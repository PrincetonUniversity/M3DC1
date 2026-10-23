from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .read_lcfs import read_lcfs
from .read_parameter import read_parameter


def _find_trace_particle(executable: str | Path | None) -> str:
    requested = (
        "trace_particle"
        if executable is None
        else str(Path(executable).expanduser())
    )
    resolved = shutil.which(requested)
    if resolved is None:
        if executable is None:
            raise FileNotFoundError(
                "trace_particle is not in PATH. Add it to PATH or pass "
                "trace_particle_executable='/path/to/trace_particle'."
            )
        raise FileNotFoundError(
            "trace_particle was not found or is not executable: "
            f"{requested}"
        )
    return resolved


def _scan_points(points: int | Sequence[int]) -> tuple[int, int]:
    values = np.asarray(points, dtype=int).reshape(-1)
    if values.size == 1:
        nx = nlambda = int(values[0])
    elif values.size == 2:
        nx, nlambda = map(int, values)
    else:
        raise ValueError("jacobian_points must be an integer or two integers.")
    if nx < 2 or nlambda < 2:
        raise ValueError("Each jacobian_points value must be at least 2.")
    return nx, nlambda


def _species_mass_charge(
    sps: int,
    field_filename: str | Path,
) -> tuple[float, float]:
    if int(sps) != sps or int(sps) not in (1, 2):
        raise ValueError("Jacobian correction requires sps=1 or sps=2.")
    if int(sps) == 1:
        mass = float(read_parameter("ion_mass", filename=field_filename))
        charge = float(read_parameter("z_ion", filename=field_filename))
    else:
        mass = float(read_parameter("fast_ion_mass", filename=field_filename))
        charge = float(read_parameter("fast_ion_z", filename=field_filename))
        if mass <= 0.0:
            mass = float(read_parameter("ion_mass", filename=field_filename))
        if charge == 0.0:
            charge = float(read_parameter("z_ion", filename=field_filename))
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("A finite positive particle mass is required.")
    if not np.isfinite(charge) or charge == 0.0:
        raise ValueError("A finite nonzero particle charge is required.")
    return mass, charge


def _read_jacobian_table(path: Path) -> np.ndarray:
    try:
        table = np.genfromtxt(path, names=True, dtype=None, encoding=None)
    except (OSError, ValueError) as error:
        raise ValueError(
            f"Could not read trace_particle Jacobian output {path}."
        ) from error
    table = np.atleast_1d(table)
    required = {
        "pphi_over_q_Wb",
        "energy_keV",
        "lambda",
        "sigma",
        "orbit_type",
        "jacobian_relative_keV_s_per_T",
        "status",
    }
    names = set(table.dtype.names or ())
    missing = sorted(required - names)
    if missing:
        raise ValueError(
            f"{path} is missing trace_particle column(s): {', '.join(missing)}."
        )
    return table


def _combine_jacobian_rows(table: np.ndarray, sigma: int):
    pphi = np.asarray(table["pphi_over_q_Wb"], dtype=float)
    lambdas = np.asarray(table["lambda"], dtype=float)
    jacobian = np.asarray(
        table["jacobian_relative_keV_s_per_T"], dtype=float
    )
    orbit_type = np.asarray(table["orbit_type"], dtype=str)
    status = np.asarray(table["status"], dtype=str)
    row_sigma = np.asarray(table["sigma"], dtype=int)
    finite = (
        np.isfinite(pphi)
        & np.isfinite(lambdas)
        & np.isfinite(jacobian)
        & (jacobian > 0.0)
        & (status == "complete")
    )
    if sigma:
        finite &= row_sigma == sigma

    pphi_nodes = np.unique(pphi[np.isfinite(pphi)])
    lambda_nodes = np.unique(lambdas[np.isfinite(lambdas)])
    values = np.full((pphi_nodes.size, lambda_nodes.size), np.nan)
    for i, pphi_value in enumerate(pphi_nodes):
        for j, lambda_value in enumerate(lambda_nodes):
            selected = (
                finite
                & np.isclose(pphi, pphi_value, rtol=1.0e-12, atol=0.0)
                & np.isclose(lambdas, lambda_value, rtol=1.0e-12, atol=0.0)
            )
            if not np.any(selected):
                continue
            if sigma:
                values[i, j] = float(np.sum(jacobian[selected]))
                continue

            passing = selected & (orbit_type == "passing")
            trapped = selected & (orbit_type == "trapped")
            combined = float(np.sum(jacobian[passing]))
            if np.any(trapped):
                # The two sigma values are different starting phases of one
                # trapped orbit, so include their phase volume only once.
                combined += float(np.mean(jacobian[trapped]))
            if combined > 0.0:
                values[i, j] = combined
    return pphi_nodes, lambda_nodes, values


def particle_com_jacobian(
    xgrid,
    lambda_grid,
    *,
    energy: float,
    timeslices: int = 0,
    field_filename: str | Path = "C1.h5",
    sps: int,
    sigma: int = 0,
    xscale: float = 1.0,
    yscale: float = 1.0,
    jacobian_points: int | Sequence[int] = 15,
    trace_particle_executable: str | Path | None = None,
    trace_processes: int = 1,
    trace_mpi_executable: str | Path = "mpirun",
    trace_dt: float = 2.0e-9,
    trace_steps: int = 50_000,
    jacobian_file: str | Path | None = None,
    trace_extra_args: Sequence[str | int | float] | None = None,
) -> np.ndarray:
    """Calculate and interpolate the relative COM phase-space Jacobian."""
    energy_value = float(energy)
    if not np.isfinite(energy_value) or energy_value <= 0.0:
        raise ValueError("A finite positive energy is required for the COM Jacobian.")
    sigma_value = int(sigma)
    if sigma_value not in (-1, 0, 1) or sigma_value != sigma:
        raise ValueError("sigma must be -1, 0, or 1.")
    xscale_value = float(xscale)
    yscale_value = float(yscale)
    if (
        not np.isfinite(xscale_value)
        or np.isclose(xscale_value, 0.0)
        or not np.isfinite(yscale_value)
        or np.isclose(yscale_value, 0.0)
    ):
        raise ValueError("xscale and yscale must be finite and nonzero.")

    xgrid = np.asarray(xgrid, dtype=float)
    lambda_grid = np.asarray(lambda_grid, dtype=float)
    if xgrid.shape != lambda_grid.shape:
        raise ValueError("xgrid and lambda_grid must have matching shapes.")
    physical_x = xgrid / xscale_value
    physical_lambda = lambda_grid / yscale_value
    finite_grid = np.isfinite(physical_x) & np.isfinite(physical_lambda)
    if not np.any(finite_grid):
        raise ValueError("The requested COM grid has no finite coordinates.")

    lambda_min = max(0.0, float(np.min(physical_lambda[finite_grid])))
    lambda_max = float(np.max(physical_lambda[finite_grid]))
    if lambda_max <= lambda_min:
        raise ValueError("The physical Lambda range must have nonzero positive span.")

    field_path = Path(field_filename).expanduser().resolve()
    if not field_path.is_file():
        raise FileNotFoundError(f"M3D-C1 field file does not exist: {field_path}")
    lcfs = read_lcfs(
        filename=field_path,
        slice=int(timeslices),
        mks=True,
        return_meta=True,
    )
    psi0 = float(lcfs.flux0)
    psi_edge = float(lcfs.psilim)
    flux_span = psi0 - psi_edge
    if not np.isfinite(flux_span) or np.isclose(flux_span, 0.0):
        raise ValueError("A finite nonzero psi0-psi_edge is required.")
    pphi_grid = psi0 + physical_x[finite_grid] * flux_span
    pphi_min = float(np.min(pphi_grid))
    pphi_max = float(np.max(pphi_grid))
    if pphi_max <= pphi_min:
        raise ValueError("The physical P_phi/q range must have nonzero span.")

    nx, nlambda = _scan_points(jacobian_points)
    mass, charge = _species_mass_charge(sps, field_path)
    process_count = int(trace_processes)
    if process_count < 1:
        raise ValueError("trace_processes must be at least 1.")

    def generate(output_path: Path) -> None:
        executable = _find_trace_particle(trace_particle_executable)
        command = [
            executable,
            "-m3dc1",
            str(field_path),
            "--timeslice",
            str(int(timeslices)),
            "--equilibrium",
            "1",
            "--perturbed",
            "0",
            "--energy",
            str(energy_value),
            "--lambda",
            str(lambda_min),
            str(lambda_max),
            str(nlambda),
            "--pphi",
            str(pphi_min),
            str(pphi_max),
            str(nx),
            "--sigma",
            str(sigma_value),
            "--mass",
            str(mass),
            "--charge",
            str(charge),
            "--dt",
            str(float(trace_dt)),
            "--steps",
            str(int(trace_steps)),
            "--output",
            str(output_path),
            "-qout",
            "1",
            "-pout",
            "0",
        ]
        if trace_extra_args is not None:
            command.extend(str(value) for value in trace_extra_args)
        if process_count > 1:
            mpi = shutil.which(str(Path(trace_mpi_executable).expanduser()))
            if mpi is None:
                raise FileNotFoundError(
                    f"MPI executable was not found: {trace_mpi_executable}"
                )
            command = [mpi, "-n", str(process_count), *command]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print("Running trace_particle:", " ".join(command))
        subprocess.run(command, cwd=output_path.parent, check=True)

    if jacobian_file is None:
        with tempfile.TemporaryDirectory(prefix="m3dc1_com_jacobian_") as directory:
            path = Path(directory) / "particle_com_jacobian.out"
            generate(path)
            table = _read_jacobian_table(path)
    else:
        path = Path(jacobian_file).expanduser().resolve()
        if not path.is_file():
            generate(path)
        table = _read_jacobian_table(path)

    table_energy = np.asarray(table["energy_keV"], dtype=float)
    energy_rows = np.isclose(
        table_energy, energy_value, rtol=1.0e-10, atol=1.0e-12
    )
    if not np.any(energy_rows):
        raise ValueError(
            f"{path} does not contain the requested energy {energy_value} keV."
        )
    table = table[energy_rows]

    pphi_nodes, lambda_nodes, jacobian_values = _combine_jacobian_rows(
        table, sigma_value
    )
    if pphi_nodes.size < 2 or lambda_nodes.size < 2:
        raise ValueError("The trace_particle output does not span a 2D COM grid.")
    if np.count_nonzero(np.isfinite(jacobian_values)) < 4:
        raise ValueError(
            "Fewer than four completed trace_particle orbits have finite Jacobians."
        )

    x_nodes = (pphi_nodes - psi0) / flux_span
    x_order = np.argsort(x_nodes)
    lambda_order = np.argsort(lambda_nodes)
    interpolator = RegularGridInterpolator(
        (x_nodes[x_order], lambda_nodes[lambda_order]),
        jacobian_values[np.ix_(x_order, lambda_order)],
        bounds_error=False,
        fill_value=np.nan,
    )
    query = np.column_stack([physical_x.ravel(), physical_lambda.ravel()])
    for column, nodes in enumerate((x_nodes[x_order], lambda_nodes[lambda_order])):
        span = float(nodes[-1] - nodes[0])
        tolerance = 1.0e-10 * max(1.0, abs(span))
        near_bounds = (
            (query[:, column] >= nodes[0] - tolerance)
            & (query[:, column] <= nodes[-1] + tolerance)
        )
        query[near_bounds, column] = np.clip(
            query[near_bounds, column], nodes[0], nodes[-1]
        )
    return np.asarray(interpolator(query), dtype=float).reshape(xgrid.shape)
