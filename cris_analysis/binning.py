"""Binning events into a spectrum, and the full raw-folder-to-spectrum pipeline."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.stats import chi2

from . import __version__
from .config import BinningSettings
from .constants import DEVICE_COLUMNS, MHZ_PER_INVERSE_CM
from .doppler import to_rest_frame
from .gating import (
    advanced_noise_cut,
    apply_reference_correction,
    apply_tof_gate,
    find_tof_gate,
    gate_wavenumber,
    simple_noise_cut,
)
from .io import load_scan
from .voltage import apply_iscool, calibrate_voltage

#: v3 converted MHz to cm^-1 with 3.34e-5 instead of 1/29979.2458 = 3.33564e-5,
#: making bins 0.13% wider than requested. Kept as an option to reproduce v3 output.
V3_CM_PER_MHZ = 0.0000334


def poisson_upper(counts: ArrayLike, alpha: float = 0.32) -> np.ndarray:
    """Upper end of the central (1 - alpha) Poisson confidence interval."""
    return chi2.ppf(1 - alpha / 2, 2 * np.asarray(counts) + 2) / 2


def poisson_lower(counts: ArrayLike, alpha: float = 0.32) -> np.ndarray:
    """Lower end of the central (1 - alpha) Poisson confidence interval (0 for zero counts)."""
    return np.nan_to_num(chi2.ppf(alpha / 2, 2 * np.asarray(counts)) / 2)


def _summarise(data: pd.DataFrame, groups, x_column: str, bunch_column: str) -> pd.DataFrame:
    df = pd.DataFrame()
    df[["x", "xerr"]] = groups[x_column].agg(["mean", "std"])
    df["xerr"] = df["xerr"].fillna(0)
    empty = groups["delta_t"].agg(lambda dt: dt.le(0).sum())
    df["y"] = groups["timestamp_copy"].count() - empty
    df["bunches"] = groups[bunch_column].nunique()
    df["yerr"] = poisson_upper(df["y"]) - df["y"]
    df["yerr_low"] = df["y"] - poisson_lower(df["y"])
    df["median_scan_time"] = np.median(data["timestamp_copy"])
    return df


def _to_relative_MHz(df: pd.DataFrame, transition_wavenumber: float) -> pd.DataFrame:
    df["x"] = (df["x"] - transition_wavenumber) * MHZ_PER_INVERSE_CM
    df["xerr"] = df["xerr"] * MHZ_PER_INVERSE_CM
    return df


def bin_by_frequency(
    data: pd.DataFrame,
    channel: str,
    transition_wavenumber: float,
    bin_size_MHz: float,
    legacy_bin_width: bool = False,
) -> pd.DataFrame:
    """Bin rest-frame wavenumbers (already Doppler corrected in ``channel``) in steps of ``bin_size_MHz``."""
    width = bin_size_MHz * (V3_CM_PER_MHZ if legacy_bin_width else 1 / MHZ_PER_INVERSE_CM)
    wn = data[channel]
    edges = np.arange(wn.min(), wn.max() + width, width)
    data = data.assign(digit_index=np.digitize(wn, edges))
    df = _summarise(data, data.groupby("digit_index"), channel, "bunch_no")
    return _to_relative_MHz(df, transition_wavenumber).sort_values(by=["x"])


def bin_by_voltage(
    data: pd.DataFrame, channel: str, transition_wavenumber: float, bin_size_voltage: float
) -> pd.DataFrame:
    """Bin voltage-scanning data in steps of the scanning voltage, then average rest-frame wavenumber per bin."""
    voltage = data["voltage"].to_numpy()
    edges = np.arange(voltage.min() - bin_size_voltage / 2, voltage.max() + bin_size_voltage / 2, bin_size_voltage)
    data = data.assign(digit_index=np.digitize(voltage, edges))
    df = _summarise(data, data.groupby("digit_index"), channel, "timestamp_copy")
    return _to_relative_MHz(df, transition_wavenumber).sort_values(by=["x"])


@dataclass
class BinnedScan:
    """A binned spectrum plus everything needed to reproduce it."""

    spectrum: pd.DataFrame
    events: pd.DataFrame
    settings: BinningSettings
    folders: list[Path]
    tof_gate: tuple[float, float]
    metadata: dict = field(default_factory=dict)

    def save(self, path: str | Path) -> None:
        """Write the spectrum as v3-compatible CSV, with a ``.json`` file of metadata beside it."""
        from .output import write_spectrum

        write_spectrum(self.spectrum, path, self.metadata)


def prepare_events(folder: str | Path, settings: BinningSettings) -> tuple[pd.DataFrame, tuple[float, float]]:
    """Run every selection step on one scan folder and return Doppler-ready events and the TOF gate used."""
    raw = load_scan(folder, settings.devices)
    data = raw.data

    if settings.noise_cut == "advanced":
        data = advanced_noise_cut(data, settings.noise_threshold, settings.signal_channel, settings.noise_channel)
    elif settings.noise_cut == "simple":
        data = simple_noise_cut(data, settings.signal_channel)
    elif settings.noise_cut != "none":
        raise ValueError(f"Unknown noise_cut '{settings.noise_cut}'")

    gate = settings.tof_gate or find_tof_gate(data, settings.tof_range, settings.tof_bin_size, settings.tof_min_width)
    data = apply_tof_gate(data, gate)

    data = apply_iscool(
        data,
        settings.iscool_method,
        settings.iscool_multiplier,
        folders=[folder],
        filename=settings.iscool_file,
        savgol_window=settings.savgol_window,
        savgol_polyorder=settings.savgol_polyorder,
    )
    data = calibrate_voltage(data, settings.voltage_scanning, settings.cec_calibration)

    if settings.apply_reference_correction:
        ref = settings.diode_wavenumber if settings.reference == "diode" else settings.hene_wavenumber
        data = apply_reference_correction(data, settings.reference_channel, ref, DEVICE_COLUMNS["wavemeter"][2:])
    data = gate_wavenumber(data, settings.wavenumber_channel, settings.wavenumber_bounds)
    return data, gate


def bin_scans(
    folders: str | Path | Iterable[str | Path], settings: BinningSettings, legacy_bin_width: bool = False
) -> BinnedScan:
    """Bin one scan folder, or several scans of the same isotope combined, into a spectrum.

    Several folders are combined the way v3 did: events from every scan are pooled and
    binned together with one set of settings.
    """
    if isinstance(folders, (str, Path)):
        folders = [folders]
    folders = [Path(f) for f in folders]

    parts, gates = [], []
    for folder in folders:
        events, gate = prepare_events(folder, settings)
        parts.append(events)
        gates.append(gate)
    events = pd.concat(parts)

    channel = settings.wavenumber_channel
    events[channel] = to_rest_frame(events[channel], events["voltage"], settings.exact_mass, settings.harmonic)

    use_voltage = settings.voltage_scanning and "cec_measured_voltage" in events.columns
    if use_voltage:
        spectrum = bin_by_voltage(events, channel, settings.transition_wavenumber, settings.bin_size_voltage)
    else:
        spectrum = bin_by_frequency(
            events, channel, settings.transition_wavenumber, settings.bin_size_MHz, legacy_bin_width
        )

    metadata = {
        "cris_analysis_version": __version__,
        "scan_folders": [str(f) for f in folders],
        "tof_gates_us": [list(g) for g in gates],
        "binned_by": "voltage" if use_voltage else "frequency",
        "legacy_bin_width": legacy_bin_width,
        "settings": settings.to_dict(),
    }
    return BinnedScan(spectrum, events, settings, folders, gates[0], metadata)
