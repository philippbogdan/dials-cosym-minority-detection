"""Counting-noise distribution of DIALS's pairwise correlations by a bootstrap over
reflections: every replicate redraws the unique reflections (Poisson weights, the same for all
crystals, so pairs keep their shared reflections) and recomputes every sigma-weighted CC
exactly as DIALS does (cf.fastcc). Crystal-to-crystal differences are not resampled, so the
replicates measure counting noise only.

Stored per trial (results/boot/<design>/<shard>.npz): the B replicate matrices minus the
observed matrix (float32)."""

from __future__ import annotations

from cf.result_io import glob_results

import glob
import hashlib
import os
import sys

import numpy as np

from cf.fastcc import dense_from_cosym, trial_arrays
from cf.paths import DATA, RESULTS, base_design, representation


def weighted_cc_matrix_w(y, v, m, wt, min_pairs=3):
    """As cf.fastcc.weighted_cc_matrix, with an extra weight per reflection (wt >= 0)."""
    n = y.shape[0]
    r = np.zeros((n, n))
    for i in range(n - 1):
        mm = m[i][None, :] & m[i + 1 :] & (wt[None, :] > 0)
        cnt = (m[i][None, :] & m[i + 1 :] & (wt[None, :] > 0)).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(mm, wt[None, :] / (v[i][None, :] + v[i + 1 :]), 0.0)
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
        ok = (cnt >= min_pairs) & (sx > 0) & (syy > 0) & np.isfinite(cc)
        r[i, i + 1 :] = np.where(ok, cc, 0.0)
    return r + r.T


def replicates(d, B=100, seed=0, rep="intensity"):
    a = trial_arrays(d, rep)
    n = len(a["names"])
    y, v, m = a["y"], a["v"], a["m"]
    rng = np.random.default_rng(seed)
    out = np.empty((B, n, n), dtype=np.float32)
    for b in range(B):
        wt = rng.poisson(1.0, y.shape[1]).astype(float)
        out[b] = weighted_cc_matrix_w(y, v, m, wt)
    return out


def parametric(d, B=100, seed=0):
    """Null-world replicates (measurement errors redrawn around the consensus, reflections
    fixed), centred: the counting-noise model of the frozen procedure v2."""
    import importlib.util

    from cf.paths import ROOT

    path = ROOT / "cf" / "frozen" / "groups_v6_1.py"
    spec = importlib.util.spec_from_file_location("pg_boot", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pg_boot"] = mod
    spec.loader.exec_module(mod)
    a = trial_arrays(d)
    return mod.parametric_replicates(a["y"], a["v"], a["m"], n_bootstrap=B, seed=seed)


def neff_matrix(y, v, m):
    """Kish effective number of reflections of every pair's sigma-weighted CC (weights
    1 / (v_i + v_j)), the quantity DIALS keeps in wij before rescaling it."""
    n = y.shape[0]
    out = np.zeros((n, n))
    for i in range(n - 1):
        mm = m[i][None, :] & m[i + 1 :]
        w = np.where(mm, 1.0 / (v[i][None, :] + v[i + 1 :]), 0.0)
        s1, s2 = w.sum(1), (w * w).sum(1)
        out[i, i + 1 :] = np.where(s2 > 0, s1**2 / np.where(s2 > 0, s2, 1.0), 0.0)
    return out + out.T


def analytic_replicates(r, neff, B=100, seed=0, min_neff=4.0):
    """Counting-noise replicates of a CC matrix from the large-sample error of Fisher's z,
    se(z_ij) = 1 / sqrt(neff_ij - 3), independent over pairs; pairs with neff <= min_neff get
    no replicate noise (and so are left out by the detector). Returns replicate - observed."""
    n = r.shape[0]
    r = np.array(r, dtype=np.float64)
    np.fill_diagonal(r, 0.0)
    ok = (neff > min_neff) & (r != 0)
    np.fill_diagonal(ok, False)
    sez = np.where(ok, 1.0 / np.sqrt(np.where(ok, neff - 3.0, 1.0)), 0.0)
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(seed)
    out = np.zeros((B, n, n), dtype=np.float32)
    for b in range(B):
        eps = np.zeros((n, n))
        eps[iu] = rng.standard_normal(len(iu[0])) * sez[iu]
        eps = eps + eps.T
        out[b] = np.where(ok, np.tanh(z + eps) - r, 0.0)
    return out


def analytic(d, B=100, seed=0, rep="intensity"):
    a = trial_arrays(d, rep)
    return analytic_replicates(a["rij"], neff_matrix(a["y"], a["v"], a["m"]), B, seed)


SUBDIR = {"reflections": "boot", "parametric": "boot_param", "analytic": "boot_analytic"}


def noise_dir(design, model):
    """Where the replicates of a design live: results/<subdir>/<design> for cryo; for rt (about
    1 MB per trial) next to the trial arrays, DATA/<subdir>/<design>, which is not copied back."""
    return (DATA if design.startswith("rt") else RESULTS) / SUBDIR[model] / design


def main(design, shard, nshards, B=100, model="reflections"):
    files = sorted(glob_results(str(DATA / "trials" / base_design(design) / "*.npz")))[shard::nshards]
    rep = representation(design)
    if rep != "intensity" and model == "parametric":
        raise ValueError("the parametric model is defined for intensities only")
    out = {}
    for f in files:
        tid = os.path.basename(f)[:-4]
        d = np.load(f)
        if "eps1" not in d.files:
            continue
        seed = int.from_bytes(hashlib.sha256(tid.encode()).digest()[:4], "little")
        if model == "parametric":
            out[tid] = parametric(d, B=B, seed=seed).astype(np.float32)
        elif model == "analytic":
            out[tid] = analytic(d, B=B, seed=seed, rep=rep)
        else:
            r_obs = trial_arrays(d, rep)["rij"]
            np.fill_diagonal(r_obs, 0.0)
            out[tid] = (replicates(d, B=B, seed=seed, rep=rep) - r_obs[None, :, :]).astype(np.float32)
    od = noise_dir(design, model)
    od.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(od / f"shard_{shard:03d}.npz", **out)
    print(len(out), "trials")


def load(design, model="reflections", ids=None):
    """{trial id: replicates - observed}; only the trials in `ids` if given."""
    res = {}
    for f in sorted(glob_results(str(noise_dir(design, model) / "shard_*.npz"))):
        d = np.load(f)
        for k in d.files:
            if ids is None or k in ids:
                res[k] = d[k]
    return res


if __name__ == "__main__":
    main(sys.argv[1], int(os.environ.get("CF_SHARD", 0)), int(os.environ.get("CF_NSHARDS", 1)),
         model=sys.argv[2] if len(sys.argv) > 2 else "reflections")
