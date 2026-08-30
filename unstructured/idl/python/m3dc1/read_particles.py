from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import elementary_charge, proton_mass

from .field_at_point import field_at_point
from .radius_matrix import radius_matrix
from .read_field import read_field
from .read_parameter import read_parameter


# Zero-based columns written by hdf5_write_particles() for vspdims=2.
PARTICLE_COLUMNS = {
    "gid": 0,
    "r": 1,
    "phi": 2,
    "z": 3,
    "weight": 4,
    "species": 5,
    "v_parallel": 6,
    "mu_over_q": 7,
    "f0": 8,
    "r0": 9,
    "phi0": 10,
    "z0": 11,
    "v_parallel0": 12,
    "mu_over_q0": 13,
}
DERIVED_PARTICLE_COLUMNS = (
    "energy",
    "xi",
    "bmag",
    "minor_radius",
    "mu",
    "ppar",
    "pperp",
    "pphi",
    "mub0overe",
)


@dataclass
class ParticleResult:
    data: np.ndarray
    columns: tuple[str | int, ...]
    filename: Path
    time: float | None
    velocity_space_dims: int | None
    marker_count: int
    selected_count: int


def _particle_path(
    timeslices: int,
    filename: str | Path | None,
    field_filename: str | Path,
) -> Path:
    if filename is not None:
        return Path(filename)
    slice_idx = int(timeslices)
    if slice_idx < 0:
        raise ValueError(f"Particle timeslices must be nonnegative, got {timeslices!r}.")
    return Path(field_filename).parent / f"ions_{slice_idx:04d}.h5"


def _species_mass_charge(species: np.ndarray, filename: str | Path) -> tuple[np.ndarray, np.ndarray]:
    ion_mass = float(read_parameter("ion_mass", filename=filename))
    ion_z = float(read_parameter("z_ion", filename=filename))
    if ion_mass <= 0.0 or ion_z == 0.0:
        raise ValueError(f"Valid ion_mass and z_ion attributes are required in {filename}.")

    fast_mass = float(read_parameter("fast_ion_mass", filename=filename))
    fast_z = float(read_parameter("fast_ion_z", filename=filename))
    if fast_mass <= 0.0:
        fast_mass = ion_mass
    if fast_z == 0.0:
        fast_z = ion_z

    species_int = np.rint(species).astype(int)
    if np.any(~np.isclose(species, species_int)) or np.any((species_int < 1) | (species_int > 2)):
        values = np.unique(species)
        raise ValueError(f"Only particle species 1 and 2 are supported, got {values}.")

    mass = np.where(species_int == 1, ion_mass, fast_mass) * proton_mass
    charge = np.where(species_int == 1, ion_z, fast_z) * elementary_charge
    return mass, charge


def _field_values_at_markers(
    r: np.ndarray,
    z: np.ndarray,
    *,
    filename: str | Path,
    timeslices: int,
    field_points: int,
    field_phi: float,
    equilibrium: bool = False,
    include_psi: bool = False,
) -> dict[str, np.ndarray]:
    if field_points < 2:
        raise ValueError("field_points must be at least 2.")
    read_kwargs = {
        "timeslices": timeslices,
        "filename": filename,
        "points": field_points,
        "phi": field_phi,
        "mks": True,
        "return_meta": True,
        "equilibrium": equilibrium,
    }
    psi_r = read_field("psi", operation=2, **read_kwargs)
    psi_z = read_field("psi", operation=3, **read_kwargs)
    fp_r = read_field("fp", operation=2, **read_kwargs)
    fp_z = read_field("fp", operation=3, **read_kwargs)
    current = read_field("I", operation=1, **read_kwargs)

    itor = int(read_parameter("itor", filename=filename))
    if itor == 1:
        poloidal_radius = radius_matrix(
            psi_r.r,
            psi_r.z,
            filename=filename,
            mks=True,
        )
    else:
        poloidal_radius = np.ones_like(np.asarray(psi_r.data), dtype=float)
    toroidal_radius = radius_matrix(
        current.r,
        current.z,
        filename=filename,
        mks=True,
    )

    bx = -np.asarray(psi_z.data) / poloidal_radius - np.asarray(fp_r.data)
    by = np.asarray(current.data) / toroidal_radius
    bz = np.asarray(psi_r.data) / poloidal_radius - np.asarray(fp_z.data)
    bmag = np.sqrt(bx**2 + by**2 + bz**2)
    values = {
        "bmag": np.asarray(field_at_point(bmag, psi_r.r, psi_r.z, r, z), dtype=float),
        "I": np.asarray(
            field_at_point(current.data, current.r, current.z, r, z),
            dtype=float,
        ),
    }
    if include_psi:
        psi = read_field("psi", operation=1, **read_kwargs)
        values["psi"] = np.asarray(
            field_at_point(psi.data, psi.r, psi.z, r, z),
            dtype=float,
        )
    return values


def _minor_radius_at_markers(
    r: np.ndarray,
    z: np.ndarray,
    *,
    filename: str | Path,
    timeslices: int,
    field_points: int,
    field_phi: float,
) -> np.ndarray:
    psi_norm = read_field(
        "psi_norm",
        timeslices=timeslices,
        filename=filename,
        points=field_points,
        equilibrium=True,
        phi=field_phi,
        mks=True,
        return_meta=True,
    )
    values = np.asarray(
        field_at_point(psi_norm.data, psi_norm.r, psi_norm.z, r, z),
        dtype=float,
    )
    return np.sqrt(np.clip(values, 0.0, None))


def _normalize_columns(columns) -> tuple[tuple[str | int, ...], set[int], set[str]]:
    if columns is None:
        requested = tuple(PARTICLE_COLUMNS)
    else:
        raw = [columns] if isinstance(columns, (str, int, np.integer)) else list(columns)
        requested = tuple(raw)
    if not requested:
        raise ValueError("At least one particle column must be requested.")

    normalized = []
    raw_indices = set()
    derived = set()
    for column in requested:
        if isinstance(column, str):
            key = column.strip().lower()
            if key in PARTICLE_COLUMNS:
                normalized.append(key)
                raw_indices.add(PARTICLE_COLUMNS[key])
            elif key in DERIVED_PARTICLE_COLUMNS:
                normalized.append(key)
                derived.add(key)
            else:
                available = ", ".join((*PARTICLE_COLUMNS, *DERIVED_PARTICLE_COLUMNS))
                raise KeyError(f"Unknown particle column {column!r}; available names: {available}.")
        else:
            index = int(column)
            if index < 0:
                raise ValueError(f"Particle column indices must be nonnegative, got {index}.")
            normalized.append(index)
            raw_indices.add(index)
    return tuple(normalized), raw_indices, derived


def _read_particles_once(
    timeslices: int = 0,
    *,
    filename: str | Path | None = None,
    columns=None,
    dataset: str = "particles/data",
    field_filename: str | Path = "C1.h5",
    sps: int | None = None,
    minor_radius: float | None = None,
    minor_radius_width: float = 0.05,
    energy: float | None = None,
    energy_width: float = 1.0,
    candidate_particles: int | None = None,
    field_points: int = 200,
    field_phi: float = 0.0,
) -> ParticleResult:
    requested, raw_indices, derived = _normalize_columns(columns)

    need_energy = energy is not None or "energy" in derived
    need_minor_radius = minor_radius is not None or "minor_radius" in derived
    need_bmag = need_energy or bool(derived & {"xi", "bmag", "pperp", "mub0overe"})
    need_pphi = "pphi" in derived
    need_mass_charge = need_energy or bool(
        derived & {"xi", "mu", "ppar", "pperp", "pphi", "mub0overe"}
    )
    if need_minor_radius or need_bmag or need_pphi:
        raw_indices.update((PARTICLE_COLUMNS["r"], PARTICLE_COLUMNS["z"]))
    if need_mass_charge:
        raw_indices.add(PARTICLE_COLUMNS["species"])
    if need_energy or derived & {"xi", "ppar", "pphi", "mub0overe"}:
        raw_indices.add(PARTICLE_COLUMNS["v_parallel"])
    if need_energy or derived & {"xi", "mu", "pperp", "mub0overe"}:
        raw_indices.add(PARTICLE_COLUMNS["mu_over_q"])
    if sps is not None:
        raw_indices.add(PARTICLE_COLUMNS["species"])

    slice_idx = int(timeslices)
    path = _particle_path(slice_idx, filename, field_filename)
    unique_columns = tuple(sorted(raw_indices))
    with h5py.File(path, "r") as h5f:
        if dataset not in h5f:
            raise KeyError(f"Dataset {dataset!r} was not found in {path}.")
        source = h5f[dataset]
        if source.ndim != 2:
            raise ValueError(f"Dataset {dataset!r} in {path} must be 2D, got shape {source.shape}.")

        marker_count, column_count = source.shape
        invalid = [column for column in unique_columns if column >= column_count]
        if invalid:
            raise ValueError(
                f"Particle column(s) {invalid} are outside the available range "
                f"0..{column_count - 1} in {path}."
            )

        if candidate_particles is None or int(candidate_particles) >= marker_count:
            rows = slice(None)
        else:
            count = int(candidate_particles)
            if count < 2:
                raise ValueError("candidate_particles must be at least 2 or None.")
            rows = np.linspace(0, marker_count - 1, count, dtype=np.int64)

        if isinstance(rows, slice):
            selected = np.asarray(source[:, unique_columns])
        else:
            selected = np.asarray(source[rows, :])[:, unique_columns]
        time = float(h5f.attrs["time"]) if "time" in h5f.attrs else None
        velocity_dims = (
            int(h5f.attrs["velocity space dims"])
            if "velocity space dims" in h5f.attrs
            else None
        )

    positions = {column: i for i, column in enumerate(unique_columns)}
    raw_data = {column: selected[:, positions[column]] for column in unique_columns}
    mask = np.ones(selected.shape[0], dtype=bool)

    if sps is not None:
        species_value = int(sps)
        if species_value not in (1, 2) or species_value != sps:
            raise ValueError(f"sps must be 1, 2, or None, got {sps!r}.")
        mask &= np.isclose(raw_data[PARTICLE_COLUMNS["species"]], species_value)

    marker_radius = None
    if need_minor_radius:
        marker_radius = _minor_radius_at_markers(
            raw_data[PARTICLE_COLUMNS["r"]],
            raw_data[PARTICLE_COLUMNS["z"]],
            filename=field_filename,
            timeslices=slice_idx,
            field_points=int(field_points),
            field_phi=float(field_phi),
        )

    if minor_radius is not None:
        radius_value = float(minor_radius)
        radius_width = float(minor_radius_width)
        if not np.isfinite(radius_value) or not 0.0 <= radius_value <= 1.0:
            raise ValueError(f"minor_radius must be between 0 and 1, got {minor_radius!r}.")
        if not np.isfinite(radius_width) or radius_width <= 0.0:
            raise ValueError(f"minor_radius_width must be positive, got {minor_radius_width!r}.")
        mask &= (
            np.isfinite(marker_radius)
            & (marker_radius <= 1.0)
            & (np.abs(marker_radius - radius_value) <= radius_width)
        )

    if not np.any(mask):
        return ParticleResult(
            data=np.empty((0, len(requested)), dtype=float),
            columns=requested,
            filename=path,
            time=time,
            velocity_space_dims=velocity_dims,
            marker_count=int(marker_count),
            selected_count=0,
        )

    derived_data = {}
    if marker_radius is not None:
        derived_data["minor_radius"] = marker_radius[mask]

    mass = None
    charge = None
    if need_mass_charge:
        species = raw_data[PARTICLE_COLUMNS["species"]][mask]
        mass, charge = _species_mass_charge(species, field_filename)

    if "mu" in derived:
        mu_over_q = raw_data[PARTICLE_COLUMNS["mu_over_q"]][mask]
        derived_data["mu"] = charge * mu_over_q

    if "ppar" in derived:
        vparallel = raw_data[PARTICLE_COLUMNS["v_parallel"]][mask]
        derived_data["ppar"] = (mass / proton_mass) * vparallel

    if need_pphi:
        if int(read_parameter("itor", filename=field_filename)) != 1:
            raise ValueError("pphi is only defined for toroidal geometry with itor=1.")
        r = raw_data[PARTICLE_COLUMNS["r"]][mask]
        z = raw_data[PARTICLE_COLUMNS["z"]][mask]
        vparallel = raw_data[PARTICLE_COLUMNS["v_parallel"]][mask]
        equilibrium_field = _field_values_at_markers(
            r,
            z,
            filename=field_filename,
            timeslices=slice_idx,
            field_points=int(field_points),
            field_phi=float(field_phi),
            equilibrium=True,
            include_psi=True,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            derived_data["pphi"] = equilibrium_field["psi"] + (
                (mass / charge)
                * vparallel
                * equilibrium_field["I"]
                / equilibrium_field["bmag"]
            )

    if need_bmag:
        r = raw_data[PARTICLE_COLUMNS["r"]][mask]
        z = raw_data[PARTICLE_COLUMNS["z"]][mask]
        field_values = _field_values_at_markers(
            r,
            z,
            filename=field_filename,
            timeslices=slice_idx,
            field_points=int(field_points),
            field_phi=float(field_phi),
        )
        bmag = field_values["bmag"]
        derived_data["bmag"] = bmag

        if need_energy or derived & {"xi", "pperp", "mub0overe"}:
            mu_over_q = raw_data[PARTICLE_COLUMNS["mu_over_q"]][mask]
            vperp_squared = 2.0 * (charge / mass) * mu_over_q * bmag
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                if "pperp" in derived:
                    derived_data["pperp"] = (mass / proton_mass) * np.sqrt(vperp_squared)
                if need_energy or derived & {"xi", "mub0overe"}:
                    vparallel = raw_data[PARTICLE_COLUMNS["v_parallel"]][mask]
                    speed_squared = vparallel**2 + vperp_squared
                    energy_joule = 0.5 * mass * speed_squared
                    derived_data["energy"] = energy_joule / (
                        1.0e3 * elementary_charge
                    )
                    derived_data["xi"] = vparallel / np.sqrt(speed_squared)
                    if "mub0overe" in derived:
                        bzero = float(read_parameter("bzero", filename=field_filename))
                        b0_norm = float(read_parameter("b0_norm", filename=field_filename))
                        reference_b = abs(bzero) * b0_norm / 1.0e4
                        if not np.isfinite(reference_b) or reference_b <= 0.0:
                            raise ValueError(
                                "Valid nonzero bzero and b0_norm attributes are required "
                                "to calculate muB0overE."
                            )
                        derived_data["mub0overe"] = (
                            charge * mu_over_q * reference_b / energy_joule
                        )

    filtered_mask = np.ones(np.count_nonzero(mask), dtype=bool)
    if energy is not None:
        energy_value = float(energy)
        width_value = float(energy_width)
        if not np.isfinite(energy_value) or energy_value < 0.0:
            raise ValueError(f"energy must be a finite nonnegative value, got {energy!r}.")
        if not np.isfinite(width_value) or width_value <= 0.0:
            raise ValueError(
                f"energy_width must be a finite positive value, got {energy_width!r}."
            )
        marker_energy = derived_data["energy"]
        filtered_mask &= np.isfinite(marker_energy)
        filtered_mask &= np.abs(marker_energy - energy_value) <= width_value

    output = []
    for column in requested:
        if isinstance(column, str) and column in DERIVED_PARTICLE_COLUMNS:
            output.append(derived_data[column][filtered_mask])
        else:
            index = PARTICLE_COLUMNS[column] if isinstance(column, str) else column
            output.append(raw_data[index][mask][filtered_mask])
    data = np.column_stack(output)

    result = ParticleResult(
        data=np.asarray(data),
        columns=requested,
        filename=path,
        time=time,
        velocity_space_dims=velocity_dims,
        marker_count=int(marker_count),
        selected_count=int(np.count_nonzero(filtered_mask)),
    )
    return result


def read_particles(
    timeslices: int = 0,
    *,
    filename: str | Path | None = None,
    columns=None,
    dataset: str = "particles/data",
    field_filename: str | Path = "C1.h5",
    sps: int | None = None,
    minor_radius: float | None = None,
    minor_radius_width: float = 0.05,
    energy: float | None = None,
    energy_width: float = 1.0,
    max_particles: int | None = None,
    field_points: int = 200,
    field_phi: float = 0.0,
    return_meta: bool = False,
):
    """Read raw or physically derived M3D-C1 particle columns.

    ``timeslices`` selects ``ions_NNNN.h5`` beside ``field_filename``. The
    keyword-only ``filename`` can override that path for nonstandard names.
    With ``columns=None``, all named raw HDF5 columns are returned without
    reading the field file. Raw columns may be selected by name or zero-based
    index. Derived columns ``energy``, ``xi``, ``bmag``, ``minor_radius``,
    ``mu``, ``ppar``, ``pperp``, ``pphi``, and ``muB0overE`` trigger only
    their required physical calculations. Names are case-insensitive. Physical
    quantities use SI units, except energy is in keV, ``ppar`` and ``pperp``
    are in ion_mass m/s, and ``muB0overE`` is dimensionless. ``pphi`` is
    canonical toroidal momentum per charge in Wb. ``minor_radius`` and
    ``energy`` are applied before limiting the result to ``max_particles``.
    """
    if max_particles is None:
        target_count = None
    else:
        target_count = int(max_particles)
        if target_count < 2 or target_count != max_particles:
            raise ValueError("max_particles must be an integer of at least 2 or None.")

    has_population_filter = minor_radius is not None or energy is not None
    candidate_count = None if has_population_filter else target_count
    result = _read_particles_once(
        timeslices,
        filename=filename,
        columns=columns,
        dataset=dataset,
        field_filename=field_filename,
        sps=sps,
        minor_radius=minor_radius,
        minor_radius_width=minor_radius_width,
        energy=energy,
        energy_width=energy_width,
        candidate_particles=candidate_count,
        field_points=field_points,
        field_phi=field_phi,
    )

    if has_population_filter and result.selected_count == 0:
        filters = []
        if sps is not None:
            filters.append(f"sps={int(sps)}")
        if minor_radius is not None:
            filters.append(f"minor_radius={float(minor_radius):g}")
        if energy is not None:
            filters.append(f"energy={float(energy):g} keV")
        raise ValueError(f"No particles matched {' and '.join(filters)} in {result.filename}.")

    if target_count is not None and result.selected_count > target_count:
        rows = np.linspace(0, result.selected_count - 1, target_count, dtype=np.int64)
        result.data = result.data[rows]
        result.selected_count = target_count

    return result if return_meta else result.data
