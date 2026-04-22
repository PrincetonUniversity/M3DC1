from __future__ import annotations

import numpy as np


def power_spectrum(signal: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Port of IDL power_spectrum.pro output layout.

    Returns:
      phi: power spectrum array
      frequency: corresponding frequency array
    """
    if signal.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    smax = np.max(np.abs(signal))
    s = signal / smax if smax > 0 else signal.copy()
    fftv = np.fft.fft(s)
    fftv = fftv * smax

    n = fftv.size
    phi_old = np.abs(fftv) ** 2
    phi = np.zeros(n, dtype=float)
    frequency = np.zeros(n, dtype=float)

    mid = n // 2
    right = int(n / 2.0 - 0.5)
    left = n - right - 1

    phi[mid] = phi_old[0]
    frequency[mid] = 0.0

    for i in range(1, left + 1):
        phi[left - i] = phi_old[n - i]
        frequency[left - i] = -i / float(t)
    for i in range(1, right + 1):
        phi[mid + i] = phi_old[i]
        frequency[mid + i] = i / float(t)

    return phi, frequency
