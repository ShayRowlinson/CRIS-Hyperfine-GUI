"""Physical constants and unit conversions used throughout the package."""

import scipy.constants as _csts

#: Speed of light expressed as MHz per cm^-1 (c in cm/s divided by 1e6).
MHZ_PER_INVERSE_CM = 29979.2458

#: Atomic mass unit in eV/c^2.
AMU_EV = _csts.physical_constants["atomic mass unit-electron volt relationship"][0]

#: Column names written by the CRIS DAQ for each device file ``<device>_ds.csv``.
DEVICE_COLUMNS = {
    "tagger": ["timestamp", "offset", "bunch_no", "events_per_bunch", "channel", "delta_t"],
    "wavemeter": [
        "timestamp",
        "offset",
        "wavenumber_1",
        "wavenumber_2",
        "wavenumber_3",
        "wavenumber_4",
        "wavenumber_mixed",
    ],
    "wavemeter_pdl": ["timestamp", "offset", "wavenumber_pdl_1"],
    "iscool2": ["timestamp", "offset", "iscool_voltage"],
    "diodes": ["timestamp", "offset", "photodiode_1", "photodiode_2", "photodiode_3", "photodiode_4"],
    "powermeter_1": ["timestamp", "offset", "powermeter_1"],
    "powermeter_3": ["timestamp", "offset", "powermeter_3"],
    "cec_voltage": ["timestamp", "offset", "cec_measured_voltage", "cec_set_voltage"],
}

#: Tagger columns that must never be forward- or back-filled (that would invent events).
EVENT_COLUMNS = ("bunch_no", "events_per_bunch", "delta_t")

#: delta_t value (in microseconds) used to mark a row as an empty bunch.
EMPTY_BUNCH_DT = -0.005

#: Tagger channel used to mark empty bunches.
EMPTY_BUNCH_CHANNEL = -1
