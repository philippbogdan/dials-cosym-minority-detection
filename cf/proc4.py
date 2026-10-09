"""Procedure candidate P4: group scan on DIALS's pairwise correlations with a counting-noise
null from a bootstrap over reflections.

1. One population: rank-1 fit rho of rij (weights wij); residuals e = r - rho rho^T.
2. Counting noise: bootstrap replicates r_b (reflections redrawn, cf.bootstrap), each refitted:
   e_b; SE_ij = sd_b(e_b,ij); noise replicates z*_b = (e_b - e) / SE.
3. Scan on z = e / SE: greedy growth from every pair, best gap per size; standardised per size
   by the noise replicates' mean and SD; T = max over sizes; p from a Gumbel fitted to the
   replicates' maxima.
4. Features of the chosen group G against the rest, for the population criterion:
   gap_raw (correlation units), within_raw, rel_gap = gap_raw / within_raw, sigma_het (excess
   spread of the rest's residuals over counting noise), R = gap_raw / sigma_het.
"""

from __future__ import annotations

import numpy as np

from cf.oracle import rank1_fit
from cf.scan import grow_all


def resid(r, w):
    w = np.array(w, float)
    np.fill_diagonal(w, 0.0)
    rho = rank1_fit(r, w)
    e = np.array(r, float) - np.outer(rho, rho)
    np.fill_diagonal(e, 0.0)
    return e


def noise_model(r, w, deltas):
    valid = np.array(w) > 0
    np.fill_diagonal(valid, False)
    e = resid(r, w)
    eb = np.empty_like(deltas, dtype=float)
    ok = np.ones(r.shape, bool)
    for b, d in enumerate(deltas):
        rb = r + d
        np.fill_diagonal(rb, 0.0)
        ok &= ~((rb == 0) & valid)  # pair dropped in a replicate: unreliable
        eb[b] = resid(rb, w)
    valid &= ok
    se = eb.std(0)
    valid &= se > 0
    noise = eb - e[None]
    return e, se, noise, valid


def per_size_best(z, valid, kmax):
    zz = np.where(valid, z, 0.0)
    order, within, outside = grow_all(zz, kmax)
    g = within - outside
    return g, order


def test(r, w, deltas, kmax=None):
    n = r.shape[0]
    kmax = n // 2 if kmax is None else min(kmax, n - 2)
    if n < 4 or kmax < 2:
        return None
    e, se, noise, valid = noise_model(r, w, deltas)
    if valid.sum() < 6:
        return None
    z = np.where(valid, e / np.where(valid, se, 1.0), 0.0)
    g, order = per_size_best(z, valid, kmax)
    obs = np.nanmax(g[:, 2:], axis=0)
    B = len(noise)
    sims = np.empty((B, len(obs)))
    for b in range(B):
        zb = np.where(valid, noise[b] / np.where(valid, se, 1.0), 0.0)
        sims[b] = np.nanmax(per_size_best(zb, valid, kmax)[0][:, 2:], axis=0)
    mu, sd = sims.mean(0), sims.std(0) + 1e-12
    t_obs = (obs - mu) / sd
    t_sim = ((sims - mu) / sd).max(1)
    beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
    loc = t_sim.mean() - 0.5772156649 * beta
    T = float(np.max(t_obs))
    p = float(-np.expm1(-np.exp(-(T - loc) / beta)))
    k = int(np.argmax(t_obs)) + 2
    path = int(np.nanargmax(g[:, k]))
    members = np.sort(order[path, :k])
    rest = np.setdiff1d(np.arange(n), members)
    ew = np.where(valid, e, np.nan)
    within = float(np.nanmean(ew[np.ix_(members, members)][np.triu_indices(k, 1)]))
    outs = np.nanmean(ew[np.ix_(members, rest)], axis=0)
    gap = float(within - np.nanmax(outs)) if len(rest) else within
    iu = np.triu_indices(len(rest), 1)
    er = ew[np.ix_(rest, rest)][iu]
    sr = np.where(valid, se, np.nan)[np.ix_(rest, rest)][iu]
    okr = np.isfinite(er) & np.isfinite(sr)
    if okr.sum() >= 3:
        mad = np.median(np.abs(er[okr] - np.median(er[okr]))) * 1.4826
        het2 = mad**2 - np.median(sr[okr] ** 2)
    else:
        het2 = np.nan
    sigma_het = float(np.sqrt(het2)) if het2 > 0 else 0.0
    return dict(p=p, T=T, k=k, members=members, gap=gap, within=within,
                rel_gap=gap / within if within > 0 else -np.inf,
                sigma_het=sigma_het, R=gap / sigma_het if sigma_het > 0 else np.inf,
                se_median=float(np.median(se[valid])))
