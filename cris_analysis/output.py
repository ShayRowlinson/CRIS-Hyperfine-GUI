"""Reading and writing spectra and results.

Spectra are written in the v3 CSV layout (semicolon separated, ``digit_index;x;xerr;y;
yerr;bunches;median_scan_time``) so the v3 GUI can still read them. Everything about how
the file was made goes into a JSON file next to it with the same name.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SPECTRUM_COLUMNS = ["x", "xerr", "y", "yerr", "bunches", "median_scan_time"]


def metadata_path(path: str | Path) -> Path:
    return Path(path).with_suffix(".json")


def write_spectrum(spectrum: pd.DataFrame, path: str | Path, metadata: dict | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    spectrum[SPECTRUM_COLUMNS].to_csv(path, sep=";")
    if metadata is not None:
        with open(metadata_path(path), "w") as fh:
            json.dump(metadata, fh, indent=2, default=str)


def read_spectrum(path: str | Path) -> pd.DataFrame:
    """Read a binned spectrum written by v3 or v4 (semicolon or comma separated)."""
    for sep in (";", ","):
        df = pd.read_csv(path, sep=sep)
        if {"x", "y", "yerr", "bunches"} <= set(df.columns):
            return df
    raise ValueError(f"{path} needs columns x, y, yerr and bunches")


def read_metadata(path: str | Path) -> dict | None:
    meta = metadata_path(path)
    if not meta.exists():
        return None
    with open(meta) as fh:
        return json.load(fh)
