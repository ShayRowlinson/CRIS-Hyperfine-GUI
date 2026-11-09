"""v4 must reproduce v3 on In-115 scan 3603 (needs the real data, see tests/README.md)."""

import json

import numpy as np
import pandas as pd

from cris_analysis import BinningSettings, bin_scans
from cris_analysis.fitting import HFSComponent, fit_spectrum

SETTINGS = {
    "mass_number": 115,
    "exact_mass": 114.903878776,
    "transition_wavenumber": 40506.421,
    "harmonic": 3,
    "bin_size_MHz": 10,
    "tof_gate": (15.5, 19.5),
}


def test_binning_matches_v3(data_dir):
    scan = bin_scans(data_dir / "scan_3603", BinningSettings(**SETTINGS), legacy_bin_width=True)
    v3 = pd.read_csv(data_dir / "v3_binned_3603.csv", sep=";").set_index("digit_index")
    spectrum = scan.spectrum
    assert list(spectrum.index) == list(v3.index)
    for column in ("y", "bunches"):
        np.testing.assert_array_equal(spectrum[column].to_numpy(), v3[column].to_numpy())
    for column in ("x", "xerr", "yerr"):
        np.testing.assert_allclose(spectrum[column].to_numpy(), v3[column].to_numpy(), atol=1e-5)


def test_fit_matches_v3(data_dir):
    spectrum = pd.read_csv(data_dir / "v3_binned_3603.csv", sep=";")
    start = HFSComponent(I=4.5, J_l=1.5, J_u=0.5, A_l=242, A_u=-300, B_l=450, centroid=-500, scale=5)
    result = fit_spectrum(spectrum, [start], background=0.02)
    v3 = json.loads((data_dir / "v3_fit_3603.json").read_text())
    params = result.parameters["hfs1"]
    for name in ("centroid", "Al", "Au", "Bl", "FWHMG", "FWHML", "scale"):
        value, error = v3[name]
        assert abs(params[name].value - value) < 0.05 * error, name
        assert abs(params[name].uncertainty - error) < 0.01 * error, name
    assert abs(result.chisqr - v3["chisqr"]) / v3["chisqr"] < 1e-5
