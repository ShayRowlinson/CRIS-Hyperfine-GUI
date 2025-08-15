# CRIS-Hyperfine-GUI
Package for quick binning of CRIS data and analysis of hyperfine structure, centroid shifts between scans and extraction of nuclear magnetic dipole and electric quadrupole moments.

Version 2 now with scan combination and isomeric fitting.
Version 3 with elementary machine learning parameter estimation

GUI is configured to be ideal for a 1080p monitor. If you're using a higher resolution, increasing the figure sizes likely a good idea

Please extract the zip file into the folder where the python files are stored.

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
