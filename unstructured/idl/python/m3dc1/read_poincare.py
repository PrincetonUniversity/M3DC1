from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass
class PoincareResult:
    data: list[np.ndarray]
    files: list[Path]
    r: np.ndarray
    z: np.ndarray
    phi: np.ndarray
    value: np.ndarray


def _natural_key(path: Path) -> tuple[str, int, str]:
    match = re.search(r"(\d+)$", path.name)
    if match is None:
        return path.name, -1, path.name
    return path.name[: match.start()], int(match.group(1)), path.name


def _resolve_files(files=None, *, directory: str | Path = ".", pattern: str = "out*") -> list[Path]:
    if files is None:
        paths = list(Path(directory).glob(pattern))
    elif isinstance(files, (str, Path)):
        p = Path(files)
        if any(ch in str(p) for ch in "*?[]"):
            paths = list(Path(directory).glob(str(p)))
        else:
            paths = [p]
    else:
        paths = [Path(p) for p in files]
    return sorted((p for p in paths if p.is_file()), key=_natural_key)


def read_poincare(
    files=None,
    *,
    directory: str | Path = ".",
    pattern: str = "out*",
    rcol: int = 1,
    zcol: int = 2,
    phicol: int | None = 0,
    valuecol: int | None = None,
    max_files: int | None = None,
    skip_empty: bool = True,
    return_meta: bool = False,
):
    """Read Poincare point files written as numeric ``out*`` text tables."""
    paths = _resolve_files(files, directory=directory, pattern=pattern)
    if max_files is not None:
        paths = paths[: max(0, int(max_files))]

    arrays: list[np.ndarray] = []
    used: list[Path] = []
    required = [int(rcol), int(zcol)]
    if phicol is not None:
        required.append(int(phicol))
    if valuecol is not None:
        required.append(int(valuecol))
    required_col = max(required)
    for path in paths:
        if skip_empty and path.stat().st_size == 0:
            continue
        try:
            arr = np.loadtxt(str(path), ndmin=2)
        except ValueError:
            if skip_empty:
                continue
            raise
        if arr.size == 0 or arr.shape[0] == 0:
            if skip_empty:
                continue
            raise ValueError(f"{path} does not contain numeric Poincare data.")
        if arr.ndim != 2 or arr.shape[1] <= required_col:
            raise ValueError(f"{path} has shape {arr.shape}; column {required_col} is required.")
        arrays.append(np.asarray(arr, dtype=float))
        used.append(path)

    if arrays:
        r = np.concatenate([a[:, int(rcol)] for a in arrays])
        z = np.concatenate([a[:, int(zcol)] for a in arrays])
        phi = np.concatenate([a[:, int(phicol)] for a in arrays]) if phicol is not None else np.asarray([], dtype=float)
        value = np.concatenate([a[:, int(valuecol)] for a in arrays]) if valuecol is not None else np.asarray([], dtype=float)
    else:
        r = np.asarray([], dtype=float)
        z = np.asarray([], dtype=float)
        phi = np.asarray([], dtype=float)
        value = np.asarray([], dtype=float)

    result = PoincareResult(data=arrays, files=used, r=r, z=z, phi=phi, value=value)
    if return_meta:
        return result
    return np.column_stack([r, z]) if r.size else np.empty((0, 2), dtype=float)
