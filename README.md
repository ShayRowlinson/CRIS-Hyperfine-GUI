# CRIS-Hyperfine-GUI
Package for quick binning of CRIS data and analysis of hyperfine structure, centroid shifts between scans and extraction of nuclear magnetic dipole and electric quadrupole moments.

Version 2 now with scan combination and isomeric fitting.
Version 3 with elementary machine learning parameter estimation

GUI is configured to be ideal for a 1080p monitor. If you're using a higher resolution, increasing the figure sizes likely a good idea

Please extract the zip file into the folder where the python files are stored.

## v4 analysis package (in progress)

The physics is moving out of the GUI into `cris_analysis`, a Python package you can use from scripts and notebooks. The v3 GUI files are unchanged and still work.

```bash
pip install -e .            # add [dev] for pytest and ruff
```

```python
from cris_analysis import BinningSettings, bin_scans
from cris_analysis.fitting import HFSComponent, fit_spectrum

settings = BinningSettings.from_toml("configs/In115_example.toml")
scan = bin_scans("path/to/scan_3603", settings)
scan.save("results/In/115/3603.csv")      # v3-format CSV plus 3603.json with every setting used

result = fit_spectrum(scan.spectrum, [HFSComponent(I=4.5, J_l=1.5, J_u=0.5, A_l=242, A_u=-300, B_l=450, centroid=-500, scale=5)])
print(result.table(), result.redchi)
```

v4 reproduces v3 binning bin for bin, and v3 fits to well within their errors, on In-115 scan 3603 (see `tests/`). Changes from v3:

- Every hard-coded setting (transition wavenumber, harmonic, ISCOOL divider, CEC calibration, noise threshold, ...) is in `BinningSettings`, loadable from a TOML file.
- Fixed: the ISCOOL `spline` and `savgol` methods crashed; the advanced noise cut dropped whole bunches (good events and the bunch count) when one event coincided with noise; device flags leaked between scans; mu could not be computed when a reference B was zero.
- Fitting: B (and A for J = 0) is fixed automatically when J is too small for it to exist, so uncertainties are no longer lost; no hidden second background; components with the same spin no longer clash.
- Bins are now exactly the requested width (v3 used 3.34e-5 cm^-1/MHz, 0.13% too wide). Pass `legacy_bin_width=True` to reproduce v3.
- scipy is pinned below 1.18 until SATLAS2 supports it.

-----------------------------------------------------------------------------------------------------------------------------------

Ensure SATLAS2 is installed!

Full functionality of GUI requires the following:

-----------------------------------------------------------------------------------------------------------------------------------

Set up a CSV file in

{program location}/Elements/'[Element symbol].csv'

e.g C:/CRIS-Hyperfine-GUI-3.0/Elements/Sb.csv

with headings: 

Mass, ExactMass, A_l, A_u, B_l, B_u, I, centroid,

Example Sb.csv left as a guide



Full functionality is possible with only Mass and ExactMass columns. Everything else can be inputted into GUI manually.

-----------------------------------------------------------------------------------------------------------------------------------

For binning and fitting:

Binning requires you to go select the scan folder where things like wavemeter_ds.csv and tagger_ds.csv are saved for that scan

Compatible Scan folder set up contains:

cec_voltage_ds, diodes_ds, iscool2_ds, powermeter_1_ds, powermeter_3_ds, tagger_ds, wavemeter_ds, wavemeter_pdl_ds

These file names are given in Dopplershift_analysis if you wish to change the file names the program looks for.


The binning process formats data ready for fitting. Format for fitting looks like:

x, xerr, y, yerr, bunches

COLUMN NAMES ARE SENSITIVE but column order doesn't matter


Fiting is to y/bunches vs x

x is relative frequency, y is counts.

If your data is already normalised to counts per bunch: simply set bunches to 1 for each data point.


If the fit or show estimate buttons don't do anything: Check your I and Js before looking into data or code.

If error pops up that uncertainties couldn't be estimated: your starting parameters are likely too far away from the final fit

-----------------------------------------------------------------------------------------------------------------------------------

For moment extraction:


Reference isotope values are required for the moment calculation.
Program calls upon the I values from the Elements csv so ensure this column is filled.


Set up a CSV file in

{program location}/References/'[Element symbol]_[Mass Number].csv'

e.g C:/CRIS-Hyperfine-GUI-3.0/References/Sb_123.csv

with headings:

A_l, A_u, B_l, B_u, mu, Q, I,

Please use mu as the header: pandas doesn't like μ

Example Sb_123.csv left as guide

-----------------------------------------------------------------------------------------------------------------------------------

For any questions: contact me at shayrowlinson@gmail.com or shay.rowlinson@student.manchester.ac.uk
