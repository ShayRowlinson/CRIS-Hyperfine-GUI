"""Doppler correction from the lab frame to the rest frame of the ions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .constants import AMU_EV


def beta(voltage: ArrayLike, exact_mass: float) -> np.ndarray:
    """Ion velocity v/c after acceleration through ``voltage`` (V) for a singly charged ion of ``exact_mass`` (amu)."""
    rest_energy = exact_mass * AMU_EV
    voltage = np.asarray(voltage, dtype=float)
    return np.sqrt(1 - rest_energy**2 / (voltage + rest_energy) ** 2)


def to_rest_frame(wavenumber: ArrayLike, voltage: ArrayLike, exact_mass: float, harmonic: int = 1) -> np.ndarray:
    """Rest-frame wavenumber for collinear geometry.

    The lab wavenumber is multiplied by ``harmonic`` first, then shifted by
    ``(1 - beta) / sqrt(1 - beta^2)``.
    """
    b = beta(voltage, exact_mass)
    return np.asarray(wavenumber, dtype=float) * harmonic * (1 - b) / np.sqrt(1 - b * b)
