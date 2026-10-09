"""cf.fastcc reproduces DIALS's own rij and neff (run in the DIALS env on saved trials)."""

from cf.result_io import glob_results


import glob

import numpy as np
import pytest

from cf.fastcc import dense_from_cosym, weighted_cc_matrix
from cf.paths import DATA

FILES = [f for f in sorted(glob_results(str(DATA / "trials" / "*" / "*.npz"))) if "eps1" in np.load(f).files][:6]


@pytest.mark.skipif(not FILES, reason="no saved trial arrays")
@pytest.mark.parametrize("path", FILES)
def test_matches_dials(path):
    d = np.load(path)
    n = len(d["names"])
    y, v, m = dense_from_cosym(d["hkl"], d["i"], d["sig"], d["ds"], n, d["eps1"])
    r, ne, _ = weighted_cc_matrix(y, v, m)
    rij = d["rij"].copy()
    np.fill_diagonal(rij, 0)
    assert np.allclose(r, rij, atol=1e-6)
    # DIALS rescales the neff weights so that their sum equals the number of non-zero
    # entries; compare up to that scale.
    w = d["wij"]
    nz = w > 0
    scale = np.count_nonzero(ne) / ne.sum()
    assert np.allclose(ne[nz] * scale, w[nz], rtol=1e-6)


@pytest.mark.skipif(not FILES, reason="no saved trial arrays")
def test_amplitude_matrix_is_pearson_of_amplitudes():
    from cf.fastcc import trial_arrays

    d = np.load(FILES[0])
    a = trial_arrays(d, "amplitude")
    y, m = a["y"], a["m"]
    for i, j in [(0, 1), (1, 2), (0, 3)]:
        mm = m[i] & m[j]
        if mm.sum() < 5:
            continue
        expected = np.corrcoef(y[i, mm], y[j, mm])[0, 1]
        assert abs(a["rij"][i, j] - expected) < 1e-9
        assert a["wij"][i, j] == mm.sum()
