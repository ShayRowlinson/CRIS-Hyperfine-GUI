import numpy as np
import pandas as pd
import pytest
from uncertainties import ufloat

from cris_analysis import BinningSettings, bin_scans, compute_moments, parse_spin, weighted_mean
from cris_analysis.binning import poisson_lower, poisson_upper
from cris_analysis.doppler import beta, to_rest_frame
from cris_analysis.fitting import HFSComponent, _build_model
from cris_analysis.gating import advanced_noise_cut, apply_tof_gate
from cris_analysis.output import read_metadata, read_spectrum
from cris_analysis.voltage import apply_iscool


def events(rows):
    df = pd.DataFrame(rows, columns=["ts", "channel", "delta_t"])
    df["timestamp_copy"] = df["ts"]
    df["bunch_no"] = df["ts"]
    return df.set_index("ts")


def test_doppler_matches_classical_beta():
    b = beta(30000.0, 115.0)
    assert b == pytest.approx(np.sqrt(2 * 30000 / (115 * 931.494e6)), rel=1e-4)
    assert to_rest_frame(10000.0, 0.0, 115.0, harmonic=3) == pytest.approx(30000.0)


def test_poisson_intervals():
    assert poisson_upper(0) == pytest.approx(1.8326, abs=1e-3)  # 68% upper limit for zero counts
    assert poisson_lower(0) == 0
    assert poisson_lower(100) < 100 < poisson_upper(100)


def test_noise_cut_keeps_good_events_in_a_noisy_bunch():
    data = events(
        [
            (1.0, 1, 17.0),
            (1.0, 1, 30.0),
            (1.0, 3, 17.05),  # one coincident event: drop only that one
            (2.0, 1, 18.0),
            (2.0, 3, 18.1),  # every row coincident: becomes an empty bunch
            (3.0, 1, 16.0),  # clean bunch
        ]
    )
    out = advanced_noise_cut(data, threshold=0.15)
    assert out.loc[[1.0], "delta_t"].tolist() == [30.0]  # v3 dropped the whole bunch here
    assert (out.loc[[2.0], "channel"] == -1).all()
    assert out.loc[[3.0], "delta_t"].tolist() == [16.0]


def test_tof_gate_turns_outside_events_into_empty_bunches():
    out = apply_tof_gate(events([(1.0, 1, 17.0), (2.0, 1, 40.0)]), (15, 20))
    assert out["channel"].tolist() == [1, -1]


def test_iscool_spline_and_savgol_run(synthetic_scan):
    """v3 crashed on these methods (it used an attribute that was never set)."""
    data = events([(1000.5, 1, 17.0), (1001.5, 1, 17.0)])
    data["iscool_voltage"] = 3.0
    for method in ("spline", "savgol"):
        out = apply_iscool(data, method, 10000.0, folders=[synthetic_scan], savgol_window=3, savgol_polyorder=1)
        assert out["iscool_voltage"].tolist() == pytest.approx([30000.0, 30000.0])


def test_bin_synthetic_scan_end_to_end(synthetic_scan, tmp_path):
    settings = BinningSettings(
        mass_number=100,
        exact_mass=100.0,
        transition_wavenumber=30000.0,
        harmonic=1,
        bin_size_MHz=10,
        tof_gate=(15, 20),
        iscool_multiplier=10000.0,
    )
    scan = bin_scans(synthetic_scan, settings)
    assert scan.spectrum["y"].sum() == 100
    assert scan.spectrum["bunches"].sum() == 200
    assert len(scan.spectrum) == 4

    path = tmp_path / "out" / "1.csv"
    scan.save(path)
    assert len(read_spectrum(path)) == 4
    assert read_metadata(path)["settings"]["transition_wavenumber"] == 30000.0


def test_scans_do_not_share_state(synthetic_scan, tmp_path):
    """v3 kept loaded-device flags on the class, so one scan could affect the next."""
    import shutil

    from cris_analysis.io import load_scan

    other = tmp_path / "scan_2"
    shutil.copytree(synthetic_scan, other)
    (other / "powermeter_1_ds.csv").write_text("1000.0;0;0.0002\n")
    assert load_scan(other).loaded_devices["powermeter_1"] is True
    assert load_scan(synthetic_scan).loaded_devices["powermeter_1"] is False


def test_settings_from_toml(tmp_path):
    path = tmp_path / "in115.toml"
    path.write_text(
        "[scan]\nmass_number = 115\nexact_mass = 114.9\ntransition_wavenumber = 40506.421\n"
        "harmonic = 3\ntof_gate = [15.5, 19.5]\n"
    )
    s = BinningSettings.from_toml(path)
    assert s.tof_gate == (15.5, 19.5) and s.harmonic == 3
    with pytest.raises(ValueError):
        BinningSettings.from_dict({"mass_number": 1, "exact_mass": 1, "transition_wavenumber": 1, "typo": 2})


def test_mu_computed_when_reference_B_is_zero():
    """v3 refused to compute mu if any reference B was zero (e.g. an I = 1/2 reference)."""
    ref = {"I": "1/2", "mu": ufloat(1.0, 0.01), "Q": 0, "A_l": 1000.0, "A_u": 100.0, "B_l": 0.0, "B_u": 0.0}
    scan = {"I": "9/2", "A_l": 200.0, "A_u": 20.0, "B_l": 50.0, "B_u": 0.0}
    m = compute_moments(scan, ref)
    assert m.mu_l.n == pytest.approx(1.0 * 200 * 4.5 / (1000 * 0.5))
    assert m.Q_l is None and m.Q_u is None


def test_parse_spin():
    assert parse_spin("9/2") == 4.5 and parse_spin(3) == 3.0


def test_weighted_mean():
    wm = weighted_mean([1.0, 3.0], [1.0, 1.0])
    assert wm.mean == 2.0 and wm.error == pytest.approx(1 / np.sqrt(2))
    assert wm.redchi == pytest.approx(2.0)


def test_B_fixed_for_J_half():
    model = _build_model(HFSComponent(I=4.5, J_l=1.5, J_u=0.5, B_u=10.0), "hfs1")
    assert model.params["Bu"].vary is False and model.params["Bu"].value == 0
    assert model.params["Bl"].vary is True


def test_two_components_with_same_spin_get_unique_names():
    from cris_analysis.fitting import fit_spectrum

    x = np.linspace(-2000, 2000, 200)
    spectrum = pd.DataFrame({"x": x, "y": np.full_like(x, 10.0), "yerr": np.full_like(x, 3.0), "bunches": 100.0})
    comps = [HFSComponent(I=4.5, J_l=1.5, J_u=0.5), HFSComponent(I=4.5, J_l=1.5, J_u=0.5, centroid=500)]
    result = fit_spectrum(spectrum, comps, background=0.1)
    assert set(result.parameters) == {"hfs1", "hfs2", "bkg"}
