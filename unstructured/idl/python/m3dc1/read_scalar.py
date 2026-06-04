from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np

from .convert_units import convert_units
from .dimensions import dimensions
from .hdf5_file_test import hdf5_file_test
from .parse_units import parse_units
from .read_parameter import read_parameter


@dataclass
class ScalarResult:
    data: np.ndarray
    title: str
    symbol: str
    units: str
    time: np.ndarray


def _cumtrapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y, dtype=float)
    if y.size <= 1:
        return out
    dx = np.diff(x)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dx)
    return out


def _lc_get(dct: dict[str, np.ndarray], key: str) -> np.ndarray | None:
    k = key.lower()
    for name, val in dct.items():
        if name.lower() == k:
            return val
    return None


def _load_scalars(filename: str | Path) -> dict[str, np.ndarray]:
    with h5py.File(str(filename), "r") as h5:
        grp = h5.get("scalars")
        if not isinstance(grp, h5py.Group):
            return {}
        out: dict[str, np.ndarray] = {}
        for k, v in grp.items():
            if isinstance(v, h5py.Dataset):
                out[k] = np.asarray(v[()])
        return out


def _load_pellets(filename: str | Path) -> dict[str, np.ndarray]:
    with h5py.File(str(filename), "r") as h5:
        grp = h5.get("pellet")
        if not isinstance(grp, h5py.Group):
            grp = h5.get("pellets")
        if not isinstance(grp, h5py.Group):
            return {}
        out: dict[str, np.ndarray] = {}
        for k, v in grp.items():
            if isinstance(v, h5py.Dataset):
                out[k] = np.asarray(v[()])
        return out


def _pick_pellet(data: np.ndarray, ipellet: int) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    if arr.ndim < 2:
        return arr
    if ipellet == -1:
        return arr
    return np.asarray(arr[int(ipellet), :], dtype=float)


def _integrate_series(data: np.ndarray, time: np.ndarray, ipellet: int) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    t = np.asarray(time, dtype=float).reshape(-1)
    if arr.ndim == 1:
        return _cumtrapz(arr, t)
    if ipellet == -1:
        out = np.zeros_like(arr, dtype=float)
        for i in range(arr.shape[0]):
            out[i, :] = _cumtrapz(arr[i, :], t[: arr.shape[1]])
        return out
    return _cumtrapz(arr.reshape(-1), t[: arr.size])


def _scalar_data_core(
    scalarname: str,
    scalars: dict[str, np.ndarray],
    filename: str | Path,
    *,
    ipellet: int = 0,
) -> tuple[np.ndarray, str, str, np.ndarray]:
    name = scalarname.lower()
    itor = int(read_parameter("itor", filename=filename))
    rzero = float(read_parameter("rzero", filename=filename))
    version = float(read_parameter("version", filename=filename))
    threed = int(read_parameter("3d", filename=filename))
    pellets = _load_pellets(filename) if version >= 31 else {}
    time = np.asarray(_lc_get(scalars, "time"), dtype=float)
    if time.size == 0:
        raise KeyError("Missing scalar 'time'.")

    d = dimensions()
    title = ""
    symbol = scalarname

    if name in {"toroidal current", "it"}:
        data = np.asarray(_lc_get(scalars, "toroidal_current"), dtype=float)
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Toroidal Current", "$I_p$", dimensions(j0=1, l0=2)
    elif name in {"plasma current", "ip"}:
        data = np.asarray(_lc_get(scalars, "toroidal_current_p"), dtype=float)
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Plasma Current", "$I_p$", dimensions(j0=1, l0=2)
    elif name in {"wall current", "iw"}:
        data = np.asarray(_lc_get(scalars, "toroidal_current_w"), dtype=float)
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Wall Current", "$I_w$", dimensions(j0=1, l0=2)
    elif name in {"total current", "itot"}:
        data = np.asarray(_lc_get(scalars, "toroidal_current_w"), dtype=float) + np.asarray(
            _lc_get(scalars, "toroidal_current"), dtype=float
        )
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Total Current", "$I_{tot}$", dimensions(j0=1, l0=2)
    elif name in {"halo current", "ih"}:
        data = np.asarray(_lc_get(scalars, "toroidal_current"), dtype=float) - np.asarray(
            _lc_get(scalars, "toroidal_current_p"), dtype=float
        )
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Toroidal Halo Current", "$I_{halo}$", dimensions(j0=1, l0=2)
    elif name in {"bootstrap current", "jbs"}:
        data = np.asarray(_lc_get(scalars, "bootstrap_current"), dtype=float)
        if itor == 0 and version < 36:
            data = data / rzero
        title, symbol, d = "Bootstrap Current", "$I_{bs}$", dimensions(j0=1, l0=2)
    elif name == "volume":
        data = np.asarray(_lc_get(scalars, "volume_p"), dtype=float)
        title, symbol, d = "Plasma Volume", "V", dimensions(l0=3)
    elif name == "toroidal flux":
        data = np.asarray(_lc_get(scalars, "toroidal_flux_p"), dtype=float)
        title, symbol, d = "Toroidal Flux", "Flux", dimensions(b0=1, l0=2)
    elif name == "reconnected flux":
        rflux = _lc_get(scalars, "reconnected_flux")
        if rflux is None:
            raise KeyError("Scalar reconnected_flux not present in scalars group.")
        data = np.abs(np.asarray(rflux, dtype=float))
        title, symbol, d = "Reconnected Flux", "psi", dimensions(b0=1, l0=1 + itor)
    elif name in {"time step", "dt"}:
        data = np.abs(np.asarray(_lc_get(scalars, "dt"), dtype=float))
        title, symbol, d = "Time Step", "dt", dimensions(t0=1)
    elif name == "psimin":
        data = np.asarray(_lc_get(scalars, "psimin"), dtype=float)
        title, symbol, d = "Psimin", "$psi_{0}$", dimensions(b0=1, l0=2)
    elif name in {"psibound", "psilim"}:
        data = np.asarray(_lc_get(scalars, "psi_lcfs"), dtype=float)
        title, symbol, d = "Psilim", "$psi_{b}$", dimensions(b0=1, l0=2)
    elif name in {"loop voltage", "vl"}:
        data = np.asarray(_lc_get(scalars, "loop_voltage"), dtype=float)
        title, symbol, d = "Loop Voltage", "$V_{L}$", dimensions(potential=1)
    elif name in {"pellet rate", "pelr"}:
        if version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_rate"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_rate")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet Rate", "$V_L$", dimensions(n0=1, l0=3, t0=-1)
    elif name in {"pellet rate d2", "pelrd2"}:
        if version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_rate_D2"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_rate_D2")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet Rate D2", "$V_L$", dimensions(n0=1, l0=3, t0=-1)
    elif name in {"pellet ablation rate", "pelablr"}:
        if version >= 26:
            raise KeyError("Pellet ablation rate is not present in this M3D-C1 version.")
        data = np.asarray(_lc_get(scalars, "pellet_ablrate"), dtype=float)
        title, symbol, d = "Pellet Ablation Rate", "$V_L$", dimensions(n0=1, t0=-1)
    elif name in {"pellet var", "pelvar"}:
        if version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_var"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_var")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet Var", "$V_L$", dimensions(l0=1)
    elif name in {"pellet radius", "pelrad"}:
        if version < 26:
            data = np.asarray(_lc_get(scalars, "r_p2"), dtype=float)
        elif version < 31:
            data = np.asarray(_lc_get(scalars, "r_p"), dtype=float)
        else:
            raw = _lc_get(pellets, "r_p")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet Radius", "$V_L$", dimensions(l0=1)
    elif name in {"pellet r position", "pelrpos"}:
        if version < 26:
            data = np.asarray(_lc_get(scalars, "pellet_x"), dtype=float)
        elif version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_r"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_r")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet R position", "$V_L$", dimensions(l0=1)
    elif name in {"pellet phi position", "pelphipos"}:
        if version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_phi"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_phi")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet phi position", "$V_L$", dimensions(l0=1)
    elif name in {"pellet z position", "pelzpos"}:
        if version < 31:
            data = np.asarray(_lc_get(scalars, "pellet_z"), dtype=float)
        else:
            raw = _lc_get(pellets, "pellet_z")
            if raw is None:
                raise KeyError("Pellet data not present in this file.")
            data = _pick_pellet(raw, ipellet)
        title, symbol, d = "Pellet Z position", "$V_L$", dimensions(l0=1)
    elif name == "beta":
        gamma = float(read_parameter("gam", filename=filename))
        data = (gamma - 1.0) * np.asarray(_lc_get(scalars, "E_P"), dtype=float) / np.maximum(
            np.asarray(_lc_get(scalars, "E_MP"), dtype=float) + np.asarray(_lc_get(scalars, "E_MT"), dtype=float),
            np.finfo(float).tiny,
        )
        title, symbol, d = ("Plasma Beta" if version >= 26 else "Global Beta"), "$\\beta$", dimensions()
    elif name in {"poloidal beta", "bp"}:
        gamma = float(read_parameter("gam", filename=filename))
        data = 2.0 * (gamma - 1.0) * np.asarray(_lc_get(scalars, "E_P"), dtype=float) / np.maximum(
            np.asarray(_lc_get(scalars, "toroidal_current"), dtype=float) ** 2,
            np.finfo(float).tiny,
        )
        title, symbol, d = ("Plasma Poloidal Beta" if version >= 26 else "Global Poloidal Beta"), "$\\beta_p$", dimensions()
    elif name in {"toroidal beta", "bt"}:
        bzero = float(read_parameter("bzero", filename=filename))
        xmag = _lc_get(scalars, "xmag")
        r0 = np.asarray(xmag if xmag is not None else np.ones_like(time) * rzero, dtype=float)
        bt0 = np.abs(bzero * rzero / np.maximum(r0, np.finfo(float).tiny))
        data = 2.0 * np.asarray(_lc_get(scalars, "Ave_P"), dtype=float) / np.maximum(bt0**2, np.finfo(float).tiny)
        title, symbol, d = "Toroidal Beta", "$\\beta_t$", dimensions()
    elif name in {"normal beta", "bn"}:
        bzero = float(read_parameter("bzero", filename=filename))
        xmag = _lc_get(scalars, "xmag")
        r0 = np.asarray(xmag if xmag is not None else np.ones_like(time) * rzero, dtype=float)
        bt0 = np.abs(bzero * rzero / np.maximum(r0, np.finfo(float).tiny))
        beta_t = 2.0 * np.asarray(_lc_get(scalars, "Ave_P"), dtype=float) / np.maximum(bt0**2, np.finfo(float).tiny)
        volume = np.asarray(_lc_get(scalars, "volume_p"), dtype=float)
        a = np.sqrt(np.maximum(volume / np.maximum(2.0 * np.pi**2 * np.maximum(r0, np.finfo(float).tiny), np.finfo(float).tiny), 0.0))
        ip = convert_units(
            np.asarray(_lc_get(scalars, "toroidal_current_p"), dtype=float),
            dimensions(j0=1, l0=2),
            filename=filename,
            cgs=False,
            mks=True,
        ) / 1.0e6
        data = 100.0 * beta_t * np.abs(bt0 * a / np.maximum(ip, np.finfo(float).tiny))
        title, symbol, d = "Normal Beta", "$\\beta_n$", dimensions()
    elif name in {"kinetic energy", "ke"}:
        data = (
            np.asarray(_lc_get(scalars, "E_KP"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_KT"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_K3"), dtype=float)
        )
        title, symbol, d = "Kinetic Energy", "KE", dimensions(energy=1)
    elif name in {"magnetic energy", "me"}:
        data = np.asarray(_lc_get(scalars, "E_MP"), dtype=float) + np.asarray(
            _lc_get(scalars, "E_MT"), dtype=float
        )
        title, symbol, d = "Magnetic Energy", "ME", dimensions(energy=1)
    elif name in {"poloidal magnetic energy", "wm"}:
        data = np.asarray(_lc_get(scalars, "E_MP"), dtype=float)
        title, symbol, d = "Poloidal Magnetic Energy", "$W_{m}$", dimensions(energy=1)
    elif name in {"thermal energy", "p"}:
        data = np.asarray(_lc_get(scalars, "E_P"), dtype=float)
        title, symbol, d = "Thermal Energy", "TE", dimensions(energy=1)
    elif name in {"electron thermal energy", "pe"}:
        data = np.asarray(_lc_get(scalars, "E_PE"), dtype=float)
        title, symbol, d = "Electron Thermal Energy", "$W_{e}$", dimensions(energy=1)
    elif name in {"ion thermal energy", "pi"}:
        data = np.asarray(_lc_get(scalars, "E_P"), dtype=float) - np.asarray(
            _lc_get(scalars, "E_PE"), dtype=float
        )
        title, symbol, d = "Ion Thermal Energy", "$W_{i}$", dimensions(energy=1)
    elif name == "energy":
        data = (
            np.asarray(_lc_get(scalars, "E_P"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_MP"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_MT"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_KP"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_KT"), dtype=float)
            + np.asarray(_lc_get(scalars, "E_K3"), dtype=float)
        )
        title, symbol, d = "Total Energy", "E", dimensions(energy=1)
    elif name in {"particles", "n"}:
        data = np.asarray(_lc_get(scalars, "particle_number"), dtype=float)
        title, symbol, d = "Particle Number", "N", dimensions(n0=1, l0=3)
    elif name in {"electrons", "ne"}:
        en = _lc_get(scalars, "electron_number")
        if en is not None:
            data = np.asarray(en, dtype=float)
        else:
            zeff = float(read_parameter("zeff", filename=filename))
            data = np.asarray(_lc_get(scalars, "particle_number"), dtype=float) * zeff
        title, symbol, d = "Electron Number", "$N_{e}$", dimensions(n0=1, l0=3)
    elif name == "angular momentum":
        data = np.asarray(_lc_get(scalars, "angular_momentum"), dtype=float)
        title, symbol, d = "Angular Momentum", "$L_{p}$", dimensions(energy=1, t0=-1)
    elif name == "time":
        data = time
        title, symbol, d = "Time", "t", dimensions(t0=1)
    elif name in {"circulation", "vorticity"}:
        data = np.asarray(_lc_get(scalars, "circulation"), dtype=float)
        title, symbol, d = "Circulation", "Circulation", dimensions(v0=1, l0=1)
    elif name == "parallel viscous heating":
        data = np.asarray(_lc_get(scalars, "parallel_viscous_heating"), dtype=float)
        title, symbol, d = "Parallel Viscous Heating", "-int dV P : Gu", dimensions(energy=1, t0=-1)
    elif name == "bwb2":
        amupar = float(read_parameter("amupar", filename=filename))
        data = 4.0 * np.asarray(_lc_get(scalars, "parallel_viscous_heating"), dtype=float) / np.maximum(
            3.0 * amupar, np.finfo(float).tiny
        )
        title, symbol, d = "", "(b.W.b)^2", dimensions(t0=-2, l0=3)
    elif name == "flux":
        data = -2.0 * np.pi * (
            np.asarray(_lc_get(scalars, "psi_lcfs"), dtype=float)
            - np.asarray(_lc_get(scalars, "psimin"), dtype=float)
        )
        if itor == 1:
            data = data * rzero
        title, symbol, d = "Flux", "Phi", dimensions(b0=1, l0=2)
    elif name == "li":
        data = -4.0 * np.pi * (
            np.asarray(_lc_get(scalars, "psi_lcfs"), dtype=float)
            - np.asarray(_lc_get(scalars, "psimin"), dtype=float)
        ) / np.asarray(_lc_get(scalars, "toroidal_current_p"), dtype=float) / rzero
        if itor == 0 and version < 36:
            data = data * rzero
        title, symbol, d = "Normalized Internal Inductance", "$l_{i}$", dimensions()
    elif name == "li3":
        data = (
            4.0
            * np.asarray(_lc_get(scalars, "w_m"), dtype=float)
            / np.asarray(_lc_get(scalars, "toroidal_current_p"), dtype=float) ** 2
            / rzero
        )
        if itor == 0 and version < 36:
            data = data * rzero**2
        title, symbol, d = "Normalized Internal Inductance", "$l_i(3)$", dimensions()
    elif name == "xmag":
        data = np.asarray(_lc_get(scalars, "xmag"), dtype=float)
        title, symbol, d = "R-Coordinate of Magnetic Axis", "$R_{0}$", dimensions(l0=1)
    elif name == "zmag":
        data = np.asarray(_lc_get(scalars, "zmag"), dtype=float)
        title, symbol, d = "Z-Coordinate of Magnetic Axis", "$Z_{0}$", dimensions(l0=1)
    elif name == "runaways":
        data = np.asarray(_lc_get(scalars, "runaways"), dtype=float)
        title, symbol, d = "Runaway Current", "$N_{RE}$", dimensions(j0=1)
    elif name in {"iz", "m_iz"}:
        data = np.asarray(_lc_get(scalars, "M_IZ"), dtype=float) / np.asarray(
            _lc_get(scalars, "toroidal_current_p"), dtype=float
        )
        title, symbol, d = "Plasma Current Centroid", "$Z_{I}$", dimensions(l0=1)
    elif name in {"radiation", "prad"}:
        data = -np.asarray(_lc_get(scalars, "radiation"), dtype=float)
        title, symbol, d = "Radiated Power", "$P_{rad}$", dimensions(p0=1, l0=3, t0=-1)
    elif name in {"line_rad", "pline"}:
        data = -np.asarray(_lc_get(scalars, "line_rad"), dtype=float)
        title, symbol, d = "Line Radiation Power", "$P_{rad}$", dimensions(p0=1, l0=3, t0=-1)
    elif name in {"brem_rad", "pbrem"}:
        data = -np.asarray(_lc_get(scalars, "brem_rad"), dtype=float)
        title, symbol, d = "Bremsstrahlung Radiation Power", "$P_{rad}$", dimensions(
            p0=1, l0=3, t0=-1
        )
    elif name in {"ion_loss", "pion"}:
        data = -np.asarray(_lc_get(scalars, "ion_loss"), dtype=float)
        title, symbol, d = "Ionization Power", "$P_{rad}$", dimensions(p0=1, l0=3, t0=-1)
    elif name in {"reck_rad", "preck"}:
        data = -np.asarray(_lc_get(scalars, "reck_rad"), dtype=float)
        title, symbol, d = "Recombination Radiation Power (Kinetic)", "$P_{rad}$", dimensions(
            p0=1, l0=3, t0=-1
        )
    elif name in {"recp_rad", "precp"}:
        data = -np.asarray(_lc_get(scalars, "recp_rad"), dtype=float)
        title, symbol, d = "Recombination Radiation Power (Potential)", "$P_{rad}$", dimensions(
            p0=1, l0=3, t0=-1
        )
    elif name == "temax":
        data = np.asarray(_lc_get(scalars, "temax"), dtype=float)
        title, symbol, d = "Maximum Te", "$max(T_{e})$", dimensions(temperature=1)
    elif name == "pohm":
        data = -(
            np.asarray(_lc_get(scalars, "e_mpd"), dtype=float)
            + np.asarray(_lc_get(scalars, "e_mtd"), dtype=float)
        )
        title, symbol, d = "Ohmic Heating", "$P_{ohm}$", dimensions(p0=1, t0=-1, l0=3)
    elif name == "ave_p":
        data = np.asarray(_lc_get(scalars, "ave_p"), dtype=float)
        title, symbol, d = "Average Pressure", "<p>", dimensions(p0=1)
    elif name == "helicity":
        data = np.asarray(_lc_get(scalars, "helicity"), dtype=float)
        title, symbol, d = "Magnetic Helicity", "H", dimensions(b0=2, l0=4)
    elif name == "pinj":
        data = np.asarray(_lc_get(scalars, "power_injected"), dtype=float)
        title, symbol, d = "Power Injected", "$P_{inj}$", dimensions(p0=1, t0=-1, l0=3)
    elif name == "kprad_n0":
        data = np.asarray(_lc_get(scalars, "kprad_n0"), dtype=float)
        title, symbol, d = "Neutral Impurities", "$n_{I0}$", dimensions(n0=1, l0=3)
    elif name == "kprad_n":
        data = np.asarray(_lc_get(scalars, "kprad_n"), dtype=float)
        title, symbol, d = "Total Impurities", "$n_{I}$", dimensions(n0=1, l0=3)
    elif name == "kprad_ion_frac":
        data = 1.0 - np.asarray(_lc_get(scalars, "kprad_n0"), dtype=float) / np.asarray(
            _lc_get(scalars, "kprad_n"), dtype=float
        )
        title, symbol, d = "Impurity Ionization Fraction", "Impurity Ionization Fraction", dimensions()
    elif name == "zeff":
        data = np.asarray(_lc_get(scalars, "electron_number"), dtype=float) / (
            np.asarray(_lc_get(scalars, "particle_number"), dtype=float)
            + np.asarray(_lc_get(scalars, "kprad_n"), dtype=float)
        )
        title, symbol, d = "Free electrons per nucleus", "$Z_{eff}$", dimensions()
    else:
        raw = _lc_get(scalars, scalarname)
        if raw is None:
            if version >= 31:
                praw = _lc_get(pellets, scalarname)
                if praw is not None:
                    data = _pick_pellet(praw, ipellet)
                    title, symbol, d = "", scalarname, dimensions()
                    if scalarname.lower().startswith("flux_thermal"):
                        symbol = "$C_t$"
                        d = dimensions(p0=1, l0=3, t0=-1)
                    return data, title, symbol, d
            raise KeyError(f"Scalar {scalarname} not recognized.")
        data = np.asarray(raw, dtype=float)
        title, symbol, d = "", scalarname, dimensions()
        if scalarname.lower().startswith("wall_force"):
            d = dimensions(p0=1, l0=2)
            if version < 27 and threed == 0:
                data = data * (2.0 * np.pi)
        if scalarname.lower().startswith("flux_thermal"):
            symbol = "$C_t$"
            d = dimensions(p0=1, l0=3, t0=-1)

    return data, title, symbol, d


def read_scalar(
    scalarname: str,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    final: bool = False,
    integrate: bool = False,
    growth: bool = False,
    ipellet: int = 0,
    cgs: bool = False,
    mks: bool = False,
    return_meta: bool = False,
    growth_rate: bool | None = None,
):
    """
    Python port of read_scalar.pro (core behavior).
    """
    if not scalarname:
        raise ValueError("Error: no scalar name provided")
    if growth_rate is not None:
        growth = bool(growth_rate)

    if isinstance(filename, (list, tuple)):
        vals = []
        for f in filename:
            x = read_scalar(
                scalarname,
                filename=f,
                final=True,
                integrate=integrate,
                growth=growth,
                ipellet=ipellet,
                cgs=cgs,
                mks=mks,
                return_meta=False,
            )
            vals.append(float(np.asarray(x).reshape(-1)[-1]))
        return np.asarray(vals, dtype=float)

    if not hdf5_file_test(filename):
        return 0.0

    scalars = _load_scalars(filename)
    if not scalars:
        return 0.0

    data, title, symbol, d = _scalar_data_core(scalarname, scalars, filename=filename, ipellet=ipellet)
    time = np.asarray(_lc_get(scalars, "time"), dtype=float)

    if integrate:
        data = _integrate_series(np.asarray(data, dtype=float), time, ipellet=ipellet)
        d = d + dimensions(t0=1)

    data = convert_units(np.asarray(data, dtype=float), d, filename=filename, cgs=cgs, mks=mks)
    time_out = convert_units(time, dimensions(t0=1), filename=filename, cgs=cgs, mks=mks)
    units = parse_units(d, cgs=cgs, mks=mks)

    if growth:
        tiny = np.finfo(float).tiny
        if np.asarray(data).ndim == 1:
            n = min(len(np.asarray(time_out).reshape(-1)), len(np.asarray(data).reshape(-1)))
            data = np.gradient(
                np.log(np.maximum(np.abs(np.asarray(data).reshape(-1)[:n]), tiny)),
                np.asarray(time_out).reshape(-1)[:n],
            )
            time_out = np.asarray(time_out).reshape(-1)[:n]
        else:
            data_arr = np.asarray(data, dtype=float)
            t_arr = np.asarray(time_out).reshape(-1)
            n = min(len(t_arr), data_arr.shape[1])
            out = np.zeros((data_arr.shape[0], n), dtype=float)
            for i in range(data_arr.shape[0]):
                out[i, :] = np.gradient(
                    np.log(np.maximum(np.abs(data_arr[i, :n]), tiny)),
                    t_arr[:n],
                )
            data = out
            time_out = t_arr[:n]
        symbol = "gamma"
        units = parse_units(dimensions(t0=-1), cgs=cgs, mks=mks)

    if final:
        arr = np.asarray(data)
        if ipellet == -1 and arr.ndim >= 2:
            out = np.asarray(arr[:, -1], dtype=float)
        else:
            out = np.asarray(arr).reshape(-1)[-1]
        if return_meta:
            return ScalarResult(
                data=np.asarray(out, dtype=float).reshape(-1),
                title=title,
                symbol=symbol,
                units=units,
                time=np.asarray([time_out.reshape(-1)[-1]], dtype=float),
            )
        if np.asarray(out).size == 1:
            return float(np.asarray(out).reshape(-1)[0])
        return np.asarray(out, dtype=float)

    if ipellet == -1:
        n = np.asarray(data).shape[-1] if np.asarray(data).ndim >= 2 else np.asarray(data).size
        time_out = np.asarray(time_out).reshape(-1)[:n]
        if np.asarray(data).ndim >= 2:
            data = np.asarray(data)[:, : time_out.size]
        else:
            data = np.asarray(data).reshape(-1)[: time_out.size]
    else:
        if np.asarray(data).ndim == 1:
            n = min(np.asarray(data).size, np.asarray(time_out).size)
            data = np.asarray(data).reshape(-1)[:n]
            time_out = np.asarray(time_out).reshape(-1)[:n]

    if return_meta:
        return ScalarResult(data=np.asarray(data), title=title, symbol=symbol, units=units, time=time_out)
    return np.asarray(data)
