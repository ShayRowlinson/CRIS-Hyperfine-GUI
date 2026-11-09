import os
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(os.environ.get("CRIS_TEST_DATA", Path(__file__).parent / "data"))


@pytest.fixture
def data_dir():
    """Folder with real beamtime data (scan_3603 and v3 baselines). Tests skip without it."""
    if not (DATA_DIR / "scan_3603").exists():
        pytest.skip(f"Real scan data not found in {DATA_DIR} (see tests/README.md)")
    return DATA_DIR


def write_device(folder: Path, device: str, rows) -> None:
    np.savetxt(folder / f"{device}_ds.csv", np.asarray(rows, dtype=float), delimiter=";")


@pytest.fixture
def synthetic_scan(tmp_path):
    """A tiny scan folder: 200 bunches, laser stepped in 4 wavenumber steps, one event every other bunch."""
    folder = tmp_path / "scan_1"
    folder.mkdir()
    tagger, wavemeter, iscool = [], [], []
    for b in range(200):
        t = 1000.0 + b * 0.01
        if b % 50 == 0:
            wavemeter.append([t - 0.001, 0, 10000.0 + (b // 50) * 0.001, 12578.84915, 0, 0, 0])
            iscool.append([t - 0.001, 0, 3.0])
        if b % 2 == 0:
            tagger.append([t, 0, b, 1, 1, 17.0 * 2000])  # signal event at 17 us
        else:
            tagger.append([t, 0, b, 0, -1, -1])  # empty bunch
    write_device(folder, "tagger", tagger)
    write_device(folder, "wavemeter", wavemeter)
    write_device(folder, "iscool2", iscool)
    return folder
