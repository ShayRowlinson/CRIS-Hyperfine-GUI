"""Settings for binning a CRIS scan.

Every number that v3 hard-coded in ``binning_v3.py`` or ``Dopplershift_analysis.py`` lives
here, with the v3 value as the default, so results match v3 unless you change something.
Settings can be built in Python or loaded from a TOML file::

    [scan]
    mass_number = 115
    exact_mass = 114.903878776
    transition_wavenumber = 40506.421
    harmonic = 3
    bin_size_MHz = 10
    tof_gate = [15.5, 19.5]
"""

from __future__ import annotations

import dataclasses
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class BinningSettings:
    """Everything needed to turn one raw scan folder into a binned spectrum."""

    # Isotope and transition
    mass_number: int
    exact_mass: float  # amu
    transition_wavenumber: float  # cm^-1, reference that x is measured from
    harmonic: int = 4  # laser frequency multiplier

    # Binning
    bin_size_MHz: float = 30.0
    bin_size_voltage: float = 4.0  # V, only used for voltage scanning
    voltage_scanning: bool = False

    # Time-of-flight gate in microseconds. None means find the gate automatically.
    tof_gate: tuple[float, float] | None = None
    tof_range: tuple[float, float] = (0.0, 1000.0)  # histogram range for the automatic gate
    tof_bin_size: float = 0.1  # microseconds
    tof_min_width: float = 1.0  # microseconds; v3 compared 10 *bins* of 0.1 us, i.e. 1 us

    # Tagger channels and noise rejection
    signal_channel: int = 1
    noise_channel: int = 3
    noise_cut: str = "advanced"  # "advanced", "simple" or "none"
    noise_threshold: float = 0.15  # microseconds

    # Wavemeter
    wavenumber_channel: str = "wavenumber_1"
    wavenumber_bounds: tuple[float, float] = (0.0, 1e15)
    apply_reference_correction: bool = False  # v3 GUI never applied it
    reference: str = "diode"  # "diode" or "hene"
    reference_channel: str = "wavenumber_2"
    diode_wavenumber: float = 12578.84915
    hene_wavenumber: float = 15798.0

    # Voltages
    iscool_method: str = "avg"  # "avg", "spline" or "savgol"
    iscool_multiplier: float = 10001.645  # voltage divider
    iscool_file: str = "iscool2"
    savgol_window: int = 51
    savgol_polyorder: int = 3
    cec_calibration: float = 1.005030

    # Which device files to read. None reads every known device that exists.
    devices: tuple[str, ...] | None = None

    extra: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> BinningSettings:
        known = {f.name for f in dataclasses.fields(cls)}
        unknown = set(values) - known
        if unknown:
            raise ValueError(f"Unknown binning settings: {', '.join(sorted(unknown))}")
        values = dict(values)
        for key in ("tof_gate", "tof_range", "wavenumber_bounds", "devices"):
            if values.get(key) is not None:
                values[key] = tuple(values[key])
        return cls(**values)

    @classmethod
    def from_toml(cls, path: str | Path, section: str = "scan") -> BinningSettings:
        with open(path, "rb") as fh:
            data = tomllib.load(fh)
        return cls.from_dict(data[section] if section else data)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)
