"""Procedure candidate P3: groups that are significant against counting noise and more
distinct than the structure inside the populations they are separated from.

On DIALS's own pairwise correlations rij (weights wij):
1. residuals e from one population (rank-1 fit), standardised by their robust scale;
2. scan: greedy growth from every pair; the best standardised gap over sizes 2..n/2;
   significance against a noise-only null (residuals resampled i.i.d., B times, max gap);
3. if significant, both parts are analysed the same way (recursively), and the split is
   kept only if its gap (in correlation units) exceeds `ratio` times the strongest gap found
   inside any final part (structure inside a population) when that is itself significant.
"""

from __future__ import annotations

import numpy as np

from cf.oracle import rank1_fit
from cf.scan import grow_all


def _resid(r, w):
    w = np.array(w, float)
    np.fill_diagonal(w, 0.0)
    rho = rank1_fit(r, w)
    e = np.array(r, float) - np.outer(rho, rho)
    valid = w > 0
    np.fill_diagonal(valid, False)
    return e, valid


def best_gap(z, valid, kmax):
    zz = np.where(valid, z, 0.0)
    order, within, outside = grow_all(zz, kmax)
    g = within - outside
    k_idx = np.nanargmax(np.nanmax(g[:, 2:], axis=0)) + 2
    path = int(np.nanargmax(g[:, k_idx]))
    return float(g[path, k_idx]), np.sort(order[path, :k_idx])


def scan_test(r, w, B=200, seed=0, null="gauss"):
    """Best group of the set and its p-value against an i.i.d. resampling null."""
    n = r.shape[0]
    kmax = n // 2
    if n < 4:
        return None
    e, valid = _resid(r, w)
    iu = np.triu_indices(n, 1)
    vals = e[iu][valid[iu]]
    if len(vals) < 3:
        return None
    s0 = np.median(np.abs(vals - np.median(vals))) * 1.4826
    if not np.isfinite(s0) or s0 <= 0:
        return None
    z = (e - np.median(vals)) / s0
    zpool = z[iu][valid[iu]]
    t_obs, members = best_gap(z, valid, kmax)
    rng = np.random.default_rng(seed)
    nulls = np.empty(B)
    for b in range(B):
        zs = np.zeros((n, n))
        if null == "resample":
            zs[iu] = rng.choice(zpool, size=len(iu[0]))
        else:
            zs[iu] = rng.standard_normal(len(iu[0]))
        zs = zs + zs.T
        nulls[b] = best_gap(zs, valid, kmax)[0]
    # Gumbel tail for a continuous p-value
    beta = np.sqrt(6.0) * nulls.std() / np.pi + 1e-12
    loc = nulls.mean() - 0.5772156649 * beta
    p = float(-np.expm1(-np.exp(-(t_obs - loc) / beta)))
    # raw gap (correlation units) of the chosen group
    rest = np.setdiff1d(np.arange(n), members)
    ew = np.where(valid, e, np.nan)
    k = len(members)
    within = np.nanmean(ew[np.ix_(members, members)][np.triu_indices(k, 1)])
    outs = np.nanmean(ew[np.ix_(members, rest)], axis=0)
    gap_raw = float(within - np.nanmax(outs)) if len(rest) else float(within)
    return dict(p=p, t=t_obs, members=members, gap_raw=gap_raw, s0=float(s0))


def analyze(r, w, alpha=0.01, alpha_inside=0.05, ratio=2.5, B=200, seed=0, depth=0, null="gauss"):
    """Returns (list of final groups as index arrays, strongest significant inside gap)."""
    n = r.shape[0]
    idx = np.arange(n)
    if n < 4 or depth > 6:
        return [idx], None
    t = scan_test(r, w, B=B, seed=seed + depth, null=null)
    if t is None:
        return [idx], None
    if t["p"] >= alpha:
        inside = t["gap_raw"] if t["p"] < alpha_inside else None
        return [idx], inside
    g = t["members"]
    rest = np.setdiff1d(idx, g)
    parts, inside = [], []
    for part in (g, rest):
        sub, ins = analyze(r[np.ix_(part, part)], w[np.ix_(part, part)], alpha, alpha_inside,
                           ratio, B, seed + 11, depth + 1, null)
        parts.extend(part[s] for s in sub)
        if ins is not None:
            inside.append(ins)
    t_int = max(inside) if inside else None
    if t_int is None or t["gap_raw"] > ratio * t_int:
        return parts, t_int
    # the split is no more distinct than structure inside the parts: one population
    return [idx], max(t["gap_raw"], t_int)


def detect(r, w, **kw):
    return analyze(r, w, **kw)[0]
