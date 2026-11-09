# Tests

`pytest` runs everything. Most tests build tiny synthetic scans and run anywhere, including GitHub Actions.

The regression tests in `test_regression_3603.py` check that v4 reproduces v3 on real data (In-115, scan 3603).
That data is not committed. To run them, put these in `tests/data/` (or point `CRIS_TEST_DATA` at a folder holding them):

- `scan_3603/` – the raw scan folder (`tagger_ds.csv`, `wavemeter_ds.csv`, `iscool2_ds.csv`, ...)
- `v3_binned_3603.csv` – v3 binning of that scan (transition 40506.421 cm^-1, harmonic 3, 10 MHz bins, TOF gate 15.5–19.5 us)
- `v3_fit_3603.json` – v3 fit of that spectrum (I = 9/2, J 3/2 -> 1/2, starting values in the test)

Without them those tests are skipped.
