"""Event selection: noise rejection, time-of-flight gate and wavenumber gate.

Rows outside a gate are not deleted. They are turned into empty-bunch markers
(channel -1, delta_t -0.005 us) so the bunch still counts towards the number of
bunches in its bin, exactly as in v3.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .constants import EMPTY_BUNCH_CHANNEL, EMPTY_BUNCH_DT


def _keep_signal_and_empty(data: pd.DataFrame, signal_channel: int) -> pd.DataFrame:
    keep = (data["channel"] == EMPTY_BUNCH_CHANNEL) | (data["channel"] == signal_channel)
    return data[keep]


def simple_noise_cut(data: pd.DataFrame, signal_channel: int = 1) -> pd.DataFrame:
    """Keep only signal events and empty-bunch markers."""
    return _keep_signal_and_empty(data, signal_channel)


def advanced_noise_cut(
    data: pd.DataFrame, threshold: float, signal_channel: int = 1, noise_channel: int = 3
) -> pd.DataFrame:
    """Reject signal events that coincide in time with a noise-channel event.

    For every bunch with both signal and noise events, each row whose delta_t lies
    within ``threshold`` (microseconds) of a noise event in the same bunch is flagged.

    * If only some rows of that bunch are flagged, the flagged rows are discarded and
      the rest of the bunch is kept.
    * If every row is flagged, the bunch is kept as an empty bunch.

    Noise-channel events everywhere become empty-bunch markers. Everything that is not
    a signal event or an empty-bunch marker is then dropped.

    This is v3's ``AdvCutNoise`` without the Python loop. v3 discarded rows by timestamp
    label, and every row of a bunch shares its timestamp, so v3 threw away the *whole*
    bunch (good events included, and the bunch itself from the bunch count) whenever
    one event coincided with noise. v3 also relied on chained assignment, which stops
    working under pandas 3.
    """
    data = data.copy()
    data["_row"] = np.arange(len(data))
    is_signal = data["channel"] == signal_channel
    is_noise = data["channel"] == noise_channel
    noise_flag = np.zeros(len(data), dtype=bool)
    discard = np.zeros(len(data), dtype=bool)

    both = np.intersect1d(data.loc[is_signal, "timestamp_copy"].unique(), data.loc[is_noise, "timestamp_copy"].unique())
    if len(both):
        rows = data.loc[data["timestamp_copy"].isin(both), ["timestamp_copy", "delta_t", "_row"]]
        noise = data.loc[is_noise & data["timestamp_copy"].isin(both), ["timestamp_copy", "delta_t"]]
        pairs = rows.merge(noise, on="timestamp_copy", suffixes=("", "_noise"))
        pairs["close"] = (pairs["delta_t"] - pairs["delta_t_noise"]).abs() < threshold
        flagged = pairs.groupby("_row")["close"].any()
        per_row = rows.set_index("_row").join(flagged.rename("flag"))
        per_bunch = per_row.groupby("timestamp_copy")["flag"].agg(["any", "all"])
        partial = per_bunch.index[per_bunch["any"] & ~per_bunch["all"]]
        flag_rows = per_row.index[per_row["flag"]].to_numpy()
        partial_rows = per_row.index[per_row["flag"] & per_row["timestamp_copy"].isin(partial)].to_numpy()
        discard[partial_rows] = True
        noise_flag[np.setdiff1d(flag_rows, partial_rows)] = True

    noise_flag |= is_noise.to_numpy()
    noise_flag &= ~discard

    data.loc[discard, "channel"] = 5  # neither signal nor empty: dropped below
    data.loc[noise_flag, "channel"] = EMPTY_BUNCH_CHANNEL
    data.loc[noise_flag, "delta_t"] = EMPTY_BUNCH_DT
    data = data.drop(columns="_row")
    return _keep_signal_and_empty(data, signal_channel)


def tof_spectrum(
    delta_t: pd.Series | np.ndarray, tof_range: tuple[float, float] = (0.0, 1000.0), bin_size: float = 0.1
) -> tuple[np.ndarray, np.ndarray]:
    """Histogram of delta_t. Returns (bin centres, counts)."""
    bins = np.arange(tof_range[0], tof_range[1] + bin_size, bin_size)
    counts, edges = np.histogram(delta_t, bins=bins)
    return edges[:-1] + np.diff(edges) / 2, counts


def find_tof_gate(
    data: pd.DataFrame,
    tof_range: tuple[float, float] = (0.0, 1000.0),
    bin_size: float = 0.1,
    min_width: float = 1.0,
    fraction: float = 0.05,
) -> tuple[float, float]:
    """Find a TOF gate around the highest peak, where counts drop below ``fraction`` of the peak."""
    centres, counts = tof_spectrum(data["delta_t"], tof_range, bin_size)
    peak = int(np.argmax(counts))
    level = fraction * counts.max()
    left = np.nonzero(counts[:peak] < level)[0][-1]
    right = np.nonzero(counts[peak:] < level)[0][0] + peak
    if centres[right] - centres[left] < min_width:
        raise RuntimeError(
            f"Automatic TOF gate is only {centres[right] - centres[left]:.2f} us wide "
            f"(minimum {min_width} us). Set the gate manually."
        )
    return float(centres[left]), float(centres[right])


def apply_tof_gate(data: pd.DataFrame, gate: tuple[float, float]) -> pd.DataFrame:
    """Turn events outside ``gate`` (microseconds) into empty-bunch markers."""
    lo, hi = gate
    data = data.copy()
    outside = (data["delta_t"] < lo) | (data["delta_t"] > hi)
    data.loc[outside, "delta_t"] = EMPTY_BUNCH_DT
    data.loc[outside, "channel"] = EMPTY_BUNCH_CHANNEL
    return data.loc[(data["delta_t"] < 0) | ((data["delta_t"] > lo) & (data["delta_t"] < hi))]


def apply_reference_correction(
    data: pd.DataFrame, reference_channel: str, reference_wavenumber: float, wavenumber_columns: list[str]
) -> pd.DataFrame:
    """Shift every wavemeter channel by the drift of the reference laser (diode or HeNe)."""
    data = data.copy()
    correction = data[reference_channel] - reference_wavenumber
    if np.all(np.abs(correction) > 100):
        warnings.warn(
            f"Reference correction is over 100 cm^-1; check the reference wavenumber ({reference_wavenumber}).",
            stacklevel=2,
        )
    for name in wavenumber_columns:
        if name != reference_channel and name in data.columns:
            data[name] = data[name] - correction
    return data


def gate_wavenumber(data: pd.DataFrame, channel: str, bounds: tuple[float, float]) -> pd.DataFrame:
    """Keep rows whose wavenumber lies strictly inside ``bounds``."""
    if channel == "wavenumber_mixed":
        data = data.copy()
        data[channel] = data["wavenumber_4"] + data["wavenumber_3"]
    if channel not in data.columns:
        raise KeyError(f"Wavenumber channel '{channel}' is not in the data")
    return data[(data[channel] > bounds[0]) & (data[channel] < bounds[1])]
