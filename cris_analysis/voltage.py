"""Beam energy: ISCOOL voltage filtering and CEC voltage calibration."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

from .io import device_path, read_device


def read_iscool(folders: Iterable[str | Path], filename: str = "iscool2") -> tuple[np.ndarray, np.ndarray]:
    """Read ISCOOL readings from one or more scan folders. Returns (timestamps, raw voltages), time-sorted."""
    frames = [read_device(f, filename) for f in folders if device_path(f, filename).exists()]
    if not frames:
        raise FileNotFoundError(f"No {filename}_ds.csv found")
    raw = pd.concat(frames).sort_values("timestamp")
    return raw["timestamp"].to_numpy(), raw["iscool_voltage"].to_numpy()


def apply_iscool(
    data: pd.DataFrame,
    method: str,
    multiplier: float,
    folders: Iterable[str | Path] = (),
    filename: str = "iscool2",
    savgol_window: int = 51,
    savgol_polyorder: int = 3,
) -> pd.DataFrame:
    """Replace the ISCOOL readout on every event with a filtered value in volts.

    ``avg``: the mean readout over all events of the scan (v3 default).
    ``spline``: linear interpolation of the raw readings in time.
    ``savgol``: Savitzky-Golay smoothing of the raw readings, then interpolation.

    ``folders`` are the scan folders whose ISCOOL files are used for ``spline`` and
    ``savgol``; passing neighbouring scans as well helps when readings are sparse.
    """
    data = data.copy()
    positive = data["iscool_voltage"] > 0
    method = method.lower()
    if method == "avg":
        data.loc[positive, "iscool_voltage"] = np.average(data["iscool_voltage"]) * multiplier
        return data
    if method not in ("spline", "spline_iscool", "savgol"):
        raise ValueError(f"Unknown ISCOOL method '{method}'. Use 'avg', 'spline' or 'savgol'.")

    times, volts = read_iscool(folders, filename)
    volts = volts * multiplier
    if method == "savgol":
        window = min(savgol_window, len(volts) - (1 - len(volts) % 2))
        volts = savgol_filter(volts, window_length=window, polyorder=min(savgol_polyorder, window - 1))
    spline = interp1d(times, volts, bounds_error=False, fill_value="extrapolate")
    data.loc[positive, "iscool_voltage"] = spline(data.loc[positive, "timestamp_copy"])
    return data


def calibrate_voltage(data: pd.DataFrame, voltage_scanning: bool, cec_calibration: float = 1.0) -> pd.DataFrame:
    """Set the ``voltage`` column to the total acceleration voltage seen by the ions.

    Laser scanning: voltage = ISCOOL. Voltage scanning: voltage = ISCOOL - CEC readout,
    with the CEC readout multiplied by 1000 (voltage divider) and ``cec_calibration``.
    """
    data = data.copy()
    if voltage_scanning:
        if "cec_measured_voltage" not in data.columns:
            raise KeyError("Voltage scanning needs cec_voltage_ds.csv in the scan folder")
        data["cec_measured_voltage"] = data["cec_measured_voltage"] * 1000 * cec_calibration
        data["voltage"] = data["iscool_voltage"] - data["cec_measured_voltage"]
    else:
        data["voltage"] = data["iscool_voltage"]
    return data
