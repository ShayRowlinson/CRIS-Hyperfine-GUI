"""cris_analysis: binning, fitting and analysis of CRIS hyperfine-structure scans.

The physics behind the CRIS Hyperfine GUI, usable from scripts and notebooks::

    from cris_analysis import BinningSettings, bin_scans

    settings = BinningSettings(mass_number=115, exact_mass=114.903878776,
                               transition_wavenumber=40506.421, harmonic=3,
                               bin_size_MHz=10, tof_gate=(15.5, 19.5))
    scan = bin_scans("data/scan_3603", settings)
    scan.save("results/In/115/3603.csv")
"""

__version__ = "4.0.0.dev0"

from .binning import BinnedScan, bin_scans, poisson_lower, poisson_upper
from .centroids import weighted_mean
from .config import BinningSettings
from .moments import compute_moments, parse_spin
from .output import read_spectrum, write_spectrum

__all__ = [
    "BinnedScan",
    "BinningSettings",
    "bin_scans",
    "compute_moments",
    "parse_spin",
    "poisson_lower",
    "poisson_upper",
    "read_spectrum",
    "weighted_mean",
    "write_spectrum",
]
