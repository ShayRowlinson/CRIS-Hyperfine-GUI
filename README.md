# CRIS-Hyperfine-GUI
Package for quick binning of CRIS data and analysis of hyperfine structure, centroid shifts between scans and extraction of nuclear magnetic dipole and electric quadrupole moments.

GUI is configured to be ideal for a 1080p monitor. If you're using a higher resolution, increasing the figure size in fitting_v1.py will make the fitting frame more clear.

Ensure SATLAS2 is installed!!!

Full functionality of GUI requires the following:

-----------------------------------------------------------------------------------------------------------------------------------

For binning and fitting:

Binning requires you to go select the scan folder where things like wavemeter_ds.csv and tagger_ds.csv are saved for that scan

Set up a CSV file in
{program location}/Elements/'[Element symbol].csv'
with headings: 
Mass, ExactMass, A_l, A_u, B_l, B_u, I, centroid,

Example Sb.csv left as a guide

Functionality is possible with only Mass, ExactMass,

If fit or show estimate buttons don't do anything: Check your I and Js. 
J values must be integer or half integer and fulfill the triangle relation

-----------------------------------------------------------------------------------------------------------------------------------

For moment extraction:

Need reference values for moment calculation
Also calls upon the I values from the Elements csv

Set up a CSV file in
{program location}/References/'[Element symbol]_[Mass Number].csv'
with headings:
A_l, A_u, B_l, B_u, mu, Q, I,
Please use mu as the header: pandas doesn't like μ

Example Sb_123.csv left as guide

-----------------------------------------------------------------------------------------------------------------------------------


