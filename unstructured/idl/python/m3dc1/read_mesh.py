from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .hdf5_file_test import hdf5_file_test
from .read_parameter import read_parameter
from .time_name import time_name


@dataclass
class MeshData:
    elements: np.ndarray
    nelms: int
    nplanes: int
    period: float


def _read_dataset_like(obj) -> np.ndarray:
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[()])
    if isinstance(obj, h5py.Group):
        for key in ("_data", "data", "value"):
            if key in obj and isinstance(obj[key], h5py.Dataset):
                return np.asarray(obj[key][()])
    raise KeyError("Unsupported mesh object layout.")


def read_mesh(filename: str | Path = "C1.h5", slice: int = 0) -> MeshData:
    """Return mesh data at a given time slice."""
    if not hdf5_file_test(filename):
        raise ValueError(f"Invalid HDF5 file: {filename}")

    group_name = time_name(int(slice))
    with h5py.File(str(filename), "r") as h5:
        if group_name not in h5:
            if "equilibrium" in h5:
                group_name = "equilibrium"
            else:
                raise KeyError(f"Missing group '{group_name}'.")

        g = h5[group_name]
        if "mesh" not in g or not isinstance(g["mesh"], h5py.Group):
            raise KeyError(f"Missing mesh group in '{group_name}'.")

        mg = g["mesh"]
        if "elements" not in mg:
            raise KeyError("Missing mesh/elements dataset.")

        elements = _read_dataset_like(mg["elements"])
        elements = np.asarray(elements, dtype=float)

        if elements.ndim != 2:
            raise ValueError(f"Unexpected mesh/elements shape: {elements.shape}")

        # IDL code uses [n_coeff, n_elem]. Some files store [n_elem, n_coeff].
        if elements.shape[0] > elements.shape[1]:
            elements = elements.T

        nelms = int(elements.shape[1])

        if "nplanes" in mg:
            nplanes = int(np.asarray(_read_dataset_like(mg["nplanes"])).reshape(-1)[0])
        else:
            nplanes = int(read_parameter("nplanes", filename=filename)) if int(read_parameter("3d", filename=filename)) == 1 else 1

        if "period" in mg:
            period = float(np.asarray(_read_dataset_like(mg["period"])).reshape(-1)[0])
        else:
            itor = int(read_parameter("itor", filename=filename))
            if itor == 1:
                period = float(2.0 * np.pi)
            else:
                rzero = float(read_parameter("rzero", filename=filename))
                period = float(2.0 * np.pi * rzero)

    return MeshData(elements=elements, nelms=nelms, nplanes=max(nplanes, 1), period=period)
