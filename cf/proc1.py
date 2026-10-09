"""Procedure candidate P1: calibrated group scan on DIALS's own pairwise correlations.

1. Single-population model: rank-1 fit rho of rij (weights wij, as cosym in 1 dimension);
   residuals e_ij = r_ij - rho_i rho_j, standardised z_ij = e_ij * sqrt(wij) / s0.
2. Scan: from every pair, greedy growth; for each size k the best gap
   g_k = max over paths [ mean z within the group - max over outsiders of mean z to the group ].
3. Null: one population with continuous crystal-to-crystal variation plus counting noise:
   e*_ij = b_i . b_j + noise_ij, with b_i ~ N(0, Theta) (Theta from the leading eigenvalues of
   e) and the noise resampled from the standardised scatter of e around its low-rank part.
   Each size's observed g_k is ranked among B simulations; the smallest per-size p-value is
   calibrated against the same simulations (min-p).
4. If the global p-value is below alpha, the group is reported and both parts are searched
   again (on their own rij), until nothing is significant.
"""

from __future__ import annotations

import numpy as np

from cf.oracle import rank1_fit
from cf.scan import grow_all


def _prep(r, w):
    n = r.shape[0]
    w = np.array(w, float)
    np.fill_diagonal(w, 0.0)
    rho = rank1_fit(r, w)
    e = np.array(r, float) - np.outer(rho, rho)
    valid = w > 0
    e = np.where(valid, e, 0.0)
    np.fill_diagonal(e, 0.0)
    return e, w, valid, rho


def _lowrank(e, valid, n_iter=10, max_dims=6, edge_factor=2.0, noise_var=None, robust=False):
    """Low-rank (heterogeneity) part of e with imputed diagonal; dims above the noise edge."""
    n = e.shape[0]
    a = e.copy()
    low = np.zeros_like(a)
    dims = 0
    for _ in range(n_iter):
        filled = np.where(valid, a, low)
        np.fill_diagonal(filled, np.diag(low))
        ev, evec = np.linalg.eigh(filled)
        idx = np.argsort(ev)[::-1]
        ev, evec = ev[idx], evec[:, idx]
        edge = edge_factor * np.sqrt(n * noise_var) if noise_var is not None else 0.0
        dims = int(min(max_dims, np.sum(ev > edge)))
        low = (evec[:, :dims] * ev[:dims]) @ evec[:, :dims].T
    theta = np.maximum(ev[:dims], 0.0) / n
    if robust and dims:
        # spread of the crystals' loadings by median absolute deviation: a few crystals
        # set apart (a minority) cannot inflate it, continuous variation sets it
        b = evec[:, :dims] * np.sqrt(np.maximum(ev[:dims], 0.0))[None, :]
        mad = np.median(np.abs(b - np.median(b, 0)[None, :]), 0) * 1.4826
        theta = mad**2
    return low, theta, dims


def gaps(z, kmax):
    order, within, outside = grow_all(z, kmax)
    g = within - outside
    best = np.nanmax(g[:, 2:], axis=0)  # per size k = 2..kmax
    return best, order, g


def calibrated_test(r, w, B=200, alpha=0.01, seed=0, kmax=None, robust=True, hetero=True):
    n = r.shape[0]
    kmax = (n // 2) if kmax is None else min(kmax, n - 2)
    if n < 4 or kmax < 2:
        return dict(p=1.0)
    e, w, valid, rho = _prep(r, w)
    iu = np.triu_indices(n, 1)
    vmask = valid[iu]
    sw = np.sqrt(w)
    # noise scale from the scatter around the low-rank part (two passes)
    s0 = np.median(np.abs(e[iu][vmask] * sw[iu][vmask])) / 0.6745
    mean_inv_w = np.mean(1.0 / w[iu][vmask])
    low, theta, dims = _lowrank(e, valid, noise_var=s0**2 * mean_inv_w, robust=robust)
    u = (e - low)[iu][vmask] * sw[iu][vmask]
    s0 = np.median(np.abs(u - np.median(u))) / 0.6745
    low, theta, dims = _lowrank(e, valid, noise_var=s0**2 * mean_inv_w, robust=robust)
    u = (e - low)[iu][vmask] * sw[iu][vmask]
    if not hetero:
        dims = 0
    if not np.isfinite(s0) or s0 <= 0:
        return dict(p=1.0)
    z = np.where(valid, e * sw / s0, 0.0)
    obs, order, g = gaps(z, kmax)
    if not np.any(np.isfinite(obs)):
        return dict(p=1.0)
    rng = np.random.default_rng(seed)
    sims = np.empty((B, len(obs)))
    for b in range(B):
        es = np.zeros((n, n))
        if dims:
            bl = rng.standard_normal((n, dims)) * np.sqrt(theta)[None, :]
            es = bl @ bl.T
        noise = np.zeros((n, n))
        noise[iu] = rng.choice(u, size=len(iu[0])) / np.where(w[iu] > 0, sw[iu], 1.0)
        es = es + noise + noise.T
        np.fill_diagonal(es, 0.0)
        zs = np.where(valid, es * sw / s0, 0.0)
        sims[b] = gaps(zs, kmax)[0]
    # standardise each size by the null's mean and SD; global statistic = max over sizes;
    # its null tail from a Gumbel fitted (by moments) to the simulations' maxima
    mu = sims.mean(0)
    sd = sims.std(0) + 1e-12
    t_obs = np.where(np.isfinite(obs), (obs - mu) / sd, -np.inf)
    t_sim = np.nan_to_num((sims - mu[None, :]) / sd[None, :], nan=-np.inf).max(1)
    t_sim = t_sim[np.isfinite(t_sim)]
    if len(t_sim) < 10:
        return dict(p=1.0)
    beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
    loc = t_sim.mean() - 0.5772156649 * beta
    tmax = float(np.max(t_obs))
    p = float(-np.expm1(-np.exp(-(tmax - loc) / beta)))
    p_emp = (1 + np.sum(t_sim >= tmax)) / (1 + len(t_sim))
    p_k = (1 + (sims >= obs[None, :]).sum(0)) / (1 + B)
    kb = int(np.argmax(t_obs)) + 2
    path = int(np.nanargmax(g[:, kb]))
    members = np.sort(order[path, :kb])
    return dict(p=float(p), p_emp=float(p_emp), k=kb, members=members, p_k=p_k, t=t_obs,
                dims=dims, s0=float(s0))


def detect(r, w, alpha=0.01, B=200, seed=0, min_size=2, depth=0, **kw):
    """Recursive calibrated detection; returns a list of groups (index arrays)."""
    n = r.shape[0]
    idx = np.arange(n)
    if n < 2 * min_size or depth > 5:
        return [idx]
    t = calibrated_test(r, w, B=B, alpha=alpha, seed=seed + depth, **kw)
    if t["p"] >= alpha:
        return [idx]
    g = t["members"]
    rest = np.setdiff1d(idx, g)
    out = []
    for part in (g, rest):
        if len(part) >= 2 * min_size:
            sub = detect(r[np.ix_(part, part)], w[np.ix_(part, part)], alpha, B, seed + 7,
                         min_size, depth + 1, **kw)
            out.extend(part[s] for s in sub)
        else:
            out.append(part)
    return out
