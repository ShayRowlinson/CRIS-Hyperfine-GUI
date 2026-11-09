"""Combining repeated measurements (centroids of several scans, moments of several fits)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass
class WeightedMean:
    mean: float
    error: float  # error on the mean, 1 / sqrt(sum of weights)
    scatter: float  # weighted standard deviation of the points (what v3 plotted as +-1 sigma)
    redchi: float  # reduced chi^2 of the points about the mean
    n: int

    @property
    def birge_error(self) -> float:
        """Error on the mean inflated by sqrt(redchi) when the points scatter more than their errors allow."""
        return self.error * max(1.0, np.sqrt(self.redchi))


def weighted_mean(values: ArrayLike, errors: ArrayLike) -> WeightedMean:
    """Inverse-variance weighted mean, ignoring points with missing or non-positive errors."""
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    ok = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    v, e = values[ok], errors[ok]
    if len(v) == 0:
        raise ValueError("No points with a valid error")
    w = 1 / e**2
    mean = float(np.sum(w * v) / np.sum(w))
    scatter = float(np.sqrt(np.average((v - mean) ** 2, weights=w)))
    redchi = float(np.sum(w * (v - mean) ** 2) / (len(v) - 1)) if len(v) > 1 else float("nan")
    return WeightedMean(mean, float(np.sqrt(1 / np.sum(w))), scatter, redchi, len(v))
