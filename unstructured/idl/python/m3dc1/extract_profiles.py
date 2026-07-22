from __future__ import annotations

import fnmatch
import tarfile
from pathlib import Path
from typing import Callable

import numpy as np


ProfileExpr = Callable[[list[float]], tuple[float, float]]


def _as_float(token: str) -> float:
    return float(str(token).replace("D", "E").replace("d", "e"))


def _read_numeric_rows(path: Path, *, skip: int = 0, d_exponent: bool = False) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh):
            if lineno < int(skip):
                continue
            text = line.strip()
            if not text:
                continue
            parts = text.split()
            if d_exponent:
                row = [_as_float(v) for v in parts]
            else:
                row = [float(v) for v in parts]
            rows.append(row)
    return rows


def _write_pairs(
    path: Path,
    rows: list[list[float]],
    expr: ProfileExpr,
    *,
    printf: bool = False,
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for row in rows:
            x, y = expr(row)
            if printf:
                fh.write(f"{x:10g}\t{y:10g}\n")
            else:
                fh.write(f"{x:g} {y:g}\n")


def _write_tail(path: Path, sources: list[Path], *, skip: int = 1) -> None:
    with path.open("w", encoding="utf-8") as out:
        for src in sources:
            with src.open("r", encoding="utf-8", errors="replace") as fh:
                for lineno, line in enumerate(fh):
                    if lineno >= int(skip):
                        out.write(line)


def _copy_file(src: Path, dst: Path) -> None:
    dst.write_bytes(src.read_bytes())


def _extract_block(path: Path, start: str) -> list[str]:
    out: list[str] = []
    active = False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if active:
                if "psinorm" in line:
                    break
                out.append(line)
                continue
            if start in line:
                active = True
    return out


def _write_block(path: Path, source: Path, start: str) -> None:
    path.write_text("".join(_extract_block(source, start)), encoding="utf-8")


def _print_c1input(*settings: str) -> None:
    print("In C1input set")
    for setting in settings:
        print(f" {setting}")


def _usage() -> None:
    print("Usage: ./extract_profiles <profile_file>")
    print("  where <profile_file> is one of: ")
    print("  * a m3dc1_profiles_*.txt from Shafer")
    print("  * a .tgz file from Osborne's phython tools")
    print("  * a p-eqdsk file")
    print("  * a ntvin.dat file")
    print("  * a NSTX_profiles.dat file")
    print("  * a pdbne*.dat, pdbte*.dat, or pdbomgeb*.dat file")


def extract_profiles(
    filename: str | Path,
    *,
    output_dir: str | Path = ".",
    profiles_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Extract profile files using the same conventions as extract_profiles.sh.

    Parameters
    ----------
    filename:
        Input profile file or tarball.
    output_dir:
        Directory where ``profile_*`` files are written. Defaults to the current
        working directory, matching the shell script.
    profiles_dir:
        Directory used to unpack Osborne tarballs. Defaults to
        ``output_dir / "profiles"``.

    Returns
    -------
    dict
        Mapping of logical profile names to written paths.
    """
    src = Path(filename)
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    def out(name: str) -> Path:
        path = outdir / name
        written[name] = path
        return path

    if not src.exists():
        _usage()
        raise FileNotFoundError(src)

    base = src.name

    if fnmatch.fnmatch(base, "m3dc1_profiles_*.txt"):
        print("Reading Shafer profiles.")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_te"), rows, lambda r: (r[0], r[1]))
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0], r[2] * 0.1))
        _write_pairs(out("profile_omega"), rows, lambda r: (r[0], r[3]))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega_ExB = 1")
        return written

    if fnmatch.fnmatch(base, "*_ipec_profs.dat"):
        print("Reading IPEC profiles.")
        rows = _read_numeric_rows(src, skip=6)
        _write_pairs(out("profile_te"), rows, lambda r: (r[0], r[4] * 1e-3))
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0], r[2] * 1e-20))
        _write_pairs(out("profile_omega"), rows, lambda r: (r[0], r[5]))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega_ExB = 1")
        return written

    if fnmatch.fnmatch(base, "*.tgz"):
        print("Reading tarball from Osborne's tools.")
        unpack_dir = Path(profiles_dir) if profiles_dir is not None else outdir / "profiles"
        unpack_dir.mkdir(parents=True, exist_ok=True)
        print("Unpacking tarball")
        with tarfile.open(src, "r:gz") as tf:
            tf.extractall(unpack_dir)
        print("Extracting profiles")
        _write_tail(out("profile_ne"), sorted(unpack_dir.glob("netanh*psi_*.dat")), skip=1)
        _write_tail(out("profile_te"), sorted(unpack_dir.glob("tetanh*psi_*.dat")), skip=1)
        _write_tail(out("profile_omega.ExB"), sorted(unpack_dir.glob("omgebspl*psi_*.dat")), skip=1)
        _write_tail(out("profile_omega.ion"), sorted(unpack_dir.glob("ommvbspl*psi_*.dat")), skip=1)
        _copy_file(outdir / "profile_omega.ExB", out("profile_omega"))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega_ExB = 1")
        return written

    if fnmatch.fnmatch(base, "*_ntvin.dat"):
        print("Reading NTVIN file")
        rows = _read_numeric_rows(src, skip=10)
        _write_pairs(out("profile_omega"), rows, lambda r: (r[1], r[2] * 1e-3))
        _write_pairs(out("profile_ne"), rows, lambda r: (r[1], r[4] * 1e-20))
        _write_pairs(out("profile_te"), rows, lambda r: (r[1], r[8] * 1e-3))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega = 1")
        return written

    if fnmatch.fnmatch(base, "k.*"):
        print("Reading k. file")
        rows = _read_numeric_rows(src, skip=10)
        _write_pairs(out("profile_omega"), rows, lambda r: (r[0], r[7] * 1e-3))
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0], r[3] * 1e-20))
        _write_pairs(out("profile_te"), rows, lambda r: (r[0], r[5] * 1e-3))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega_ExB = 1")
        return written

    if fnmatch.fnmatch(base, "NSTX*_profiles.dat"):
        print("Reading NSTX profiles file")
        rows = _read_numeric_rows(src, skip=6)
        _write_pairs(out("profile_omega"), rows, lambda r: (r[0] * r[0], r[5]))
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0] * r[0], r[1] * 1e-20))
        _write_pairs(out("profile_te"), rows, lambda r: (r[0] * r[0], r[2] * 1e-3))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega = 1")
        return written

    if fnmatch.fnmatch(base, "pdbne*.dat"):
        print("Reading pdbne file")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_ne_rho_0"), rows, lambda r: (r[0], r[1] * 1e-6))
        _print_c1input("iread_ne = 4")
        return written

    if fnmatch.fnmatch(base, "pdbte*.dat"):
        print("Reading pdbte file")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_te_rho_3"), rows, lambda r: (r[0], r[1]))
        _print_c1input("iread_te = 4")
        return written

    if fnmatch.fnmatch(base, "pdbomgeb*.dat"):
        print("Reading pdbomgeb file")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_omega_rho_0"), rows, lambda r: (r[0], r[1] * 1e3))
        _print_c1input("iread_omega_ExB = 4")
        return written

    if fnmatch.fnmatch(base, "p*.*"):
        print("Reading p-eqdsk file")
        _write_block(out("profile_ne"), src, "psinorm ne")
        _write_block(out("profile_te"), src, "psinorm te")
        _write_block(out("profile_omega.ExB"), src, "psinorm omgeb")
        _write_block(out("profile_omega.C"), src, "psinorm omgvb")
        _write_block(out("profile_omega.ion"), src, "psinorm ommvb")
        _write_block(out("profile_omega.electron"), src, "psinorm omevb")
        _copy_file(outdir / "profile_omega.ExB", out("profile_omega"))
        _print_c1input("iread_ne = 1", "iread_te = 1", "iread_omega_ExB = 1")
        return written

    if fnmatch.fnmatch(base, "neprof_*.asc"):
        print("Reading neprof.asc file")
        rows = _read_numeric_rows(src, d_exponent=True)
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0] * r[0], r[1] * 1e-20))
        _print_c1input("iread_ne = 1")
        return written

    if fnmatch.fnmatch(base, "Teprof_*.asc"):
        print("Reading Teprof.asc file")
        rows = _read_numeric_rows(src, d_exponent=True)
        _write_pairs(out("profile_te"), rows, lambda r: (r[0] * r[0], r[1] * 1e-3))
        _print_c1input("iread_te = 1")
        return written

    if fnmatch.fnmatch(base, "vtprof_*.asc"):
        print("Reading vtprof.asc file")
        rows = _read_numeric_rows(src, d_exponent=True)
        _write_pairs(out("profile_vphi"), rows, lambda r: (r[0] * r[0], r[1]))
        _print_c1input("iread_omega = 3")
        return written

    if fnmatch.fnmatch(base, "*_profiles.dat"):
        print("Reading AUG profile file")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_ne"), rows, lambda r: (r[0] * r[0], r[7] * 1e-4), printf=True)
        _write_pairs(out("profile_te"), rows, lambda r: (r[0] * r[0], r[9] * 1e-3), printf=True)
        _print_c1input("iread_te = 1", "iread_ne = 1")
        return written

    if fnmatch.fnmatch(base, "*_vtor_*.dat"):
        print("Reading AUG rotation profile file")
        rows = _read_numeric_rows(src, skip=1)
        _write_pairs(out("profile_omega"), rows, lambda r: (r[2], r[3] * 1e-3 / r[0]), printf=True)
        _print_c1input("iread_omega = 1")
        return written

    _usage()
    raise ValueError(f"Unsupported profile file: {src}")


def convert_p(filename: str | Path = "profile_p") -> np.ndarray:
    """Convert pressure columns in ``profile_p`` from kPa to Pa in place."""
    path = Path(filename)
    profile = np.loadtxt(path, dtype=float)
    if profile.ndim == 1:
        profile = profile.reshape(1, -1)
    if profile.shape[1] < 3:
        raise ValueError("profile_p must have at least three columns.")
    profile[:, 1] = profile[:, 1] * 1000.0
    profile[:, 2] = profile[:, 2] * 1000.0
    print(profile)
    with path.open("w", encoding="utf-8") as pf:
        for ln in profile[:, :3]:
            pf.write(" {0:8.6f}   {1:>9.3f}   {2:>9.3f}\n".format(*ln))
    return profile


def smooth_profile(psin, values, psirange):
    psirange = np.asarray(psirange)
    from scipy.interpolate import interp1d
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt

    from . import fpylib as fpyl

    profile = np.column_stack([psin, values])
    profile_old = np.copy(profile)

    if len(psirange.shape) == 1:
        psirange = np.asarray([psirange])

    psi_replace = []
    psi_replace_ind = []
    for i, data in enumerate(profile):
        for interval in psirange:
            if data[0] >= interval[0] and data[0] <= interval[1]:
                psi_replace.append(data[0])
                psi_replace_ind.append(i)

    profile_removed = np.delete(profile, psi_replace_ind, 0)

    f = interp1d(profile_removed[:, 0], profile_removed[:, 1], kind="cubic")
    new_values = f(psi_replace)

    for i in range(len(psi_replace_ind)):
        profile[psi_replace_ind[i], 1] = new_values[i]

    fig = plt.figure(constrained_layout=True, figsize=(10, 5))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
    f2_ax1 = fig.add_subplot(spec2[0, 0])
    f2_ax2 = fig.add_subplot(spec2[0, 1])

    f2_ax1.plot(profile_removed[:, 0], profile_removed[:, 1], lw=0, marker=".", label="original (points removed)")
    f2_ax1.plot(profile_old[:, 0], profile_old[:, 1], ls="--", lw=2, label="original")
    f2_ax1.plot(profile[:, 0], profile[:, 1], lw=2, label="smoothed profile")
    f2_ax1.set_xlabel(r"$\psi_N$")
    f2_ax1.set_ylabel("Rotation profile")
    f2_ax1.grid()

    f2_ax2.plot(profile[:, 0], fpyl.deriv(profile[:, 1], profile[:, 0]))
    f2_ax2.set_xlabel(r"$\psi_N$")
    f2_ax2.set_ylabel("Rotation profile derivative")
    f2_ax2.grid()

    return profile
