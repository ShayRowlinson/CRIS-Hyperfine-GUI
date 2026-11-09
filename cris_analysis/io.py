"""Reading raw CRIS scan folders.

A scan folder holds one ``<device>_ds.csv`` per device (tagger, wavemeter, iscool2, ...),
each a semicolon-separated table with no header.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import DEVICE_COLUMNS, EVENT_COLUMNS


@dataclass
class RawScan:
    """All device readings of one scan merged into one time-ordered table."""

    data: pd.DataFrame
    folder: Path
    loaded_devices: dict[str, bool] = field(default_factory=dict)


def device_path(folder: str | Path, device: str) -> Path:
    return Path(folder) / f"{device}_ds.csv"


def read_device(folder: str | Path, device: str) -> pd.DataFrame:
    """Read one device file into a DataFrame with named columns."""
    return pd.read_csv(device_path(folder, device), sep=";", header=None, names=DEVICE_COLUMNS[device])


def load_scan(folder: str | Path, devices: tuple[str, ...] | None = None) -> RawScan:
    """Read every device file in ``folder`` and merge them on timestamp.

    Device readings (wavemeter, ISCOOL, ...) are forward- then back-filled onto the
    tagger rows, so each tagger event carries the latest reading of every device.
    Rows that are not tagger events are then dropped. ``delta_t`` is converted from
    tagger ticks to microseconds.
    """
    folder = Path(folder)
    devices = tuple(DEVICE_COLUMNS) if devices is None else devices

    frames = []
    loaded: dict[str, bool] = {}
    for device in devices:
        path = device_path(folder, device)
        if not path.exists():
            loaded[device] = False
            continue
        frame = read_device(folder, device).drop(columns=["offset"])
        frames.append(frame)
        loaded[device] = True

    if not loaded.get("tagger"):
        raise FileNotFoundError(f"No tagger_ds.csv found in {folder}")

    data = pd.concat(frames)
    data["voltage"] = np.zeros(len(data["timestamp"]))
    data["timestamp_copy"] = data["timestamp"]
    data = data.set_index(["timestamp"])
    data = data.sort_index()

    # Fill device readings onto every row, but never the event columns.
    for device in devices:
        for column in DEVICE_COLUMNS[device]:
            if column in EVENT_COLUMNS or column not in data.columns:
                continue
            data[column] = data[column].ffill().bfill()

    data = data.dropna()
    data["delta_t"] = data["delta_t"] / 2000  # tagger ticks -> microseconds
    return RawScan(data=data, folder=folder, loaded_devices=loaded)
