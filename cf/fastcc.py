"""Pairwise sigma-weighted correlation coefficients and Kish effective sample sizes, as
dials.algorithms.scaling.scaling_library.ExtendedDatasetStatistics.weighted_cchalf
computes them for the cosym rij matrix, vectorised over a dense crystal x reflection
layout. Validated against DIALS's own rij on the trial data (tests/test_fastcc.py)."""

from __future__ import annotations

import numpy as np


def dense_from_cosym(hkl, i, sig, ds, n, eps1=None):
    """Dense intensity, variance and observed-mask matrices (n crystals x unique hkl).

    `eps1` (bool per row) keeps only reflections with epsilon == 1 in the Patterson group,
    the only ones DIALS's cosym matcher pairs."""
    if eps1 is not None:
        hkl, i, sig, ds = hkl[eps1], i[eps1], sig[eps1], ds[eps1]
    key = (hkl[:, 0].astype(np.int64) * 4096 + hkl[:, 1]) * 4096 + hkl[:, 2]
    uniq, col = np.unique(key, return_inverse=True)
    y = np.zeros((n, len(uniq)))
    v = np.zeros((n, len(uniq)))
    m = np.zeros((n, len(uniq)), bool)
    y[ds, col] = i
    v[ds, col] = sig**2
    m[ds, col] = True
    return y, v, m


def trial_arrays(d, representation="intensity"):
    """Names, rij, wij and the dense arrays of a saved trial, aligned with the rows of DIALS's
    matrix. cosym drops a data set with fewer than its minimum number of reflections without
    renumbering the others; the matrix rows are then the data sets still present in `ds`, in
    order, and the crystal dropped is in no row (it keeps whatever OPTICS gave it).

    representation="amplitude": the same reflections as amplitudes, sqrt(max(I, 0)), with unit
    weights, and the matrix recomputed from them (Pearson's r over the shared reflections; wij
    the number of shared reflections where r is defined); y, v are then the amplitudes and unit
    variances, so every noise model computed from them refers to this matrix."""
    names = [str(x) for x in d["names"]]
    rij = np.array(d["rij"], dtype=np.float64)
    wij = np.array(d["wij"], dtype=np.float64)
    ds = np.asarray(d["ds"])
    rows = np.unique(ds)
    if rij.shape[0] != len(names):
        assert len(rows) == rij.shape[0], (len(rows), rij.shape, len(names))
        names = [names[i] for i in rows]
        ds = np.searchsorted(rows, ds)
    y, v, m = dense_from_cosym(d["hkl"], d["i"], d["sig"], ds, len(names), d["eps1"])
    if representation == "amplitude":
        y = np.where(m, np.sqrt(np.maximum(y, 0.0)), 0.0)
        v = m.astype(np.float64)
        rij, neff, _ = weighted_cc_matrix(y, v, m)
        wij = np.where(rij != 0, neff, 0.0)
    elif representation != "intensity":
        raise ValueError(representation)
    return dict(names=names, rij=rij, wij=wij, y=y, v=v, m=m)


def weighted_cc_matrix(y, v, m, min_pairs: int = 3):
    """Return (rij, neff, npairs) with rij = 0 and neff = 0 where incalculable."""
    n = y.shape[0]
    r = np.zeros((n, n))
    ne = np.zeros((n, n))
    npairs = np.zeros((n, n), int)
    for i in range(n - 1):
        mm = m[i][None, :] & m[i + 1 :]
        cnt = mm.sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(mm, 1.0 / (v[i][None, :] + v[i + 1 :]), 0.0)
            sw = w.sum(1)
            nw = w / sw[:, None]
            xbar = (nw * y[i][None, :]).sum(1)
            ybar = (nw * y[i + 1 :]).sum(1)
            dx = np.where(mm, y[i][None, :] - xbar[:, None], 0.0)
            dy = np.where(mm, y[i + 1 :] - ybar[:, None], 0.0)
            sxy = (nw * dx * dy).sum(1)
            sx = (nw * dx * dx).sum(1)
            syy = (nw * dy * dy).sum(1)
            cc = sxy / np.sqrt(sx * syy)
            neff = 1.0 / (nw**2).sum(1)
        ok = (cnt >= min_pairs) & (sx > 0) & (syy > 0) & np.isfinite(cc)
        r[i, i + 1 :] = np.where(ok, cc, 0.0)
        ne[i, i + 1 :] = np.where(ok, neff, 0.0)
        npairs[i, i + 1 :] = cnt
    return r + r.T, ne + ne.T, npairs + npairs.T
