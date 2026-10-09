"""The oracle and decoy measurements of PREREGISTRATION.md (amendment 1, items 3 and 4).

Input is what DIALS computes its pairwise correlations from: each crystal's merged,
normalised, resolution-filtered intensities with sigmas (saved per trial by cf.run_trials
--save-arrays), restricted to the reflections DIALS's matcher pairs (epsilon == 1).

c_ij: crystals i and j are compared with the inverse-variance weighted mean of all other
crystals of the trial; each deviation is divided by its standard deviation; c_ij is the
correlation of the two deviations over the reflections both observe. C(G) is the mean of c_ij
over the pairs of a group: how much its members deviate from the rest in the same way. A
crystal without diffraction deviates a lot but not in the same way as anyone else, so it does
not inflate C (unlike a chi-square of group means).
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from cf.fastcc import dense_from_cosym


def deviation_correlations(y, v, m):
    n, nu = y.shape
    w = np.where(m, 1.0 / np.where(m, v, 1.0), 0.0)
    W = w.sum(0)
    S = (w * y).sum(0)
    c = np.zeros((n, n))
    npairs = np.zeros((n, n), int)
    for i in range(n - 1):
        j = np.arange(i + 1, n)
        Wij = W[None, :] - w[i][None, :] - w[j]
        Sij = S[None, :] - (w[i] * y[i])[None, :] - w[j] * y[j]
        ok = m[i][None, :] & m[j] & (Wij > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ref = np.where(ok, Sij / np.where(ok, Wij, 1.0), 0.0)
            vref = np.where(ok, 1.0 / np.where(ok, Wij, 1.0), 0.0)
            di = np.where(ok, (y[i][None, :] - ref) / np.sqrt(v[i][None, :] + vref), 0.0)
            dj = np.where(ok, (y[j] - ref) / np.sqrt(v[j] + vref), 0.0)
            num = (di * dj).sum(1)
            den = np.sqrt((di * di).sum(1) * (dj * dj).sum(1))
            cc = np.where(den > 0, num / den, 0.0)
        c[i, j] = cc
        npairs[i, j] = ok.sum(1)
    return c + c.T, npairs + npairs.T


def from_npz(d):
    n = len(d["names"])
    y, v, m = dense_from_cosym(d["hkl"], d["i"], d["sig"], d["ds"], n, d["eps1"])
    return deviation_correlations(y, v, m)


def group_means(c, subsets: np.ndarray) -> np.ndarray:
    """C for each row of `subsets` (B x k index array): mean of c over the pairs."""
    k = subsets.shape[1]
    block = c[subsets[:, :, None], subsets[:, None, :]]
    return (block.sum((1, 2)) - np.einsum("bii->b", block)) / (k * (k - 1))


def _all_subsets(pool, k, chunk=200_000):
    it = itertools.combinations(pool, k)
    while True:
        block = list(itertools.islice(it, chunk))
        if not block:
            return
        yield np.array(block, dtype=np.int64)


def oracle_p(c, minority, majority, max_subsets=20_000, seed=0):
    """p-value of C(minority) among size-k subsets of the majority (exact when feasible)."""
    minority = np.asarray(minority)
    majority = np.asarray(majority)
    k = len(minority)
    c_obs = group_means(c, minority[None, :])[0]
    total = math.comb(len(majority), k)
    if total <= max_subsets:
        vals = np.concatenate([group_means(c, s) for s in _all_subsets(majority, k)])
    else:
        rng = np.random.default_rng(seed)
        subs = np.array([rng.choice(majority, k, replace=False) for _ in range(max_subsets)])
        vals = group_means(c, subs)
    p = (1 + np.sum(vals >= c_obs)) / (1 + len(vals))
    return dict(C=float(c_obs), p=float(p), n_null=int(len(vals)), exact=total <= max_subsets,
                null_q95=float(np.quantile(vals, 0.95)))


def decoy_max(c, majority, k, exhaustive_limit=2_000_000, seed=0):
    """Largest C over size-k subsets of the majority: exact when there are at most
    `exhaustive_limit` subsets, otherwise a lower bound from greedy starts + swap search."""
    majority = np.asarray(majority)
    total = math.comb(len(majority), k)
    if total <= exhaustive_limit:
        best = -np.inf
        best_set = None
        for s in _all_subsets(majority, k):
            vals = group_means(c, s)
            j = int(np.argmax(vals))
            if vals[j] > best:
                best, best_set = float(vals[j]), s[j]
        if best_set is None:  # no subset has a calculable pair
            return dict(max=-np.inf, members=[], exact=True)
        return dict(max=best, members=[int(x) for x in best_set], exact=True)
    best, best_set = -np.inf, None
    sub = c[np.ix_(majority, majority)]
    m = len(majority)
    for a in range(m):
        cur = [a]
        while len(cur) < k:
            rest = [x for x in range(m) if x not in cur]
            gain = sub[np.ix_(rest, cur)].sum(1)
            cur.append(rest[int(np.argmax(gain))])
        improved = True
        while improved:
            improved = False
            for pos in range(k):
                others = cur[:pos] + cur[pos + 1 :]
                cand = [x for x in range(m) if x not in cur]
                gains = sub[np.ix_(cand, others)].sum(1)
                old = sub[cur[pos], others].sum()
                j = int(np.argmax(gains))
                if gains[j] > old + 1e-15:
                    cur[pos] = cand[j]
                    improved = True
        val = float(group_means(sub, np.array([cur]))[0])
        if val > best:
            best, best_set = val, [int(majority[x]) for x in cur]
    return dict(max=best, members=best_set, exact=False)


# --- amendment 2: residuals of DIALS's rij from a single population (rank-1 fit) ---------


def rank1_fit(r, w, iters: int = 5000, tol: float = 1e-12):
    """Weighted least-squares rank-1 fit rho (as cosym in one dimension), diagonal excluded.
    Each step moves halfway to every rho_i's least-squares value given the others (the plain
    simultaneous update oscillates when the rho_i differ widely)."""
    r = np.array(r, dtype=float)
    w = np.array(w, dtype=float)
    np.fill_diagonal(w, 0.0)
    wsum = w.sum(1)
    rho = np.sqrt(np.clip(np.where(wsum > 0, (w * r).sum(1) / np.maximum(wsum, 1e-300), 0), 1e-6, None))
    for _ in range(iters):
        num = (w * r * rho[None, :]).sum(1)
        den = (w * rho[None, :] ** 2).sum(1)
        new = 0.5 * (np.where(den > 0, num / np.maximum(den, 1e-300), 0.0) + rho)
        if np.max(np.abs(new - rho)) < tol:
            rho = new
            break
        rho = new
    return rho


def residuals(r, w):
    """e_ij = r_ij - rho_i rho_j, with NaN where DIALS could not calculate r_ij (w == 0)."""
    rho = rank1_fit(r, w)
    e = np.array(r, dtype=float) - np.outer(rho, rho)
    e[np.array(w) <= 0] = np.nan
    np.fill_diagonal(e, np.nan)
    return e, rho


def group_means_nan(e, subsets: np.ndarray) -> np.ndarray:
    """S for each row of `subsets`: mean of e over the pairs, ignoring missing pairs."""
    block = e[subsets[:, :, None], subsets[:, None, :]]
    valid = np.isfinite(block)
    tot = np.where(valid, block, 0.0).sum((1, 2))
    cnt = valid.sum((1, 2))
    return np.where(cnt > 0, tot / np.maximum(cnt, 1), -np.inf)


def oracle_p_resid(e, minority, majority, max_subsets=20_000, seed=0):
    minority = np.asarray(minority)
    majority = np.asarray(majority)
    k = len(minority)
    s_obs = group_means_nan(e, minority[None, :])[0]
    total = math.comb(len(majority), k)
    if total <= max_subsets:
        vals = np.concatenate([group_means_nan(e, s) for s in _all_subsets(majority, k)])
    else:
        rng = np.random.default_rng(seed)
        subs = np.array([rng.choice(majority, k, replace=False) for _ in range(max_subsets)])
        vals = group_means_nan(e, subs)
    p = (1 + np.sum(vals >= s_obs)) / (1 + len(vals))
    return dict(S=float(s_obs), p=float(p), n_null=int(len(vals)), exact=total <= max_subsets,
                null_q95=float(np.quantile(vals, 0.95)))


def decoy_max_resid(e, majority, k, exhaustive_limit=2_000_000):
    """Largest S over size-k subsets of the majority (exact, or a local-search lower bound)."""
    majority = np.asarray(majority)
    total = math.comb(len(majority), k)
    if total <= exhaustive_limit:
        best, best_set = -np.inf, None
        for s in _all_subsets(majority, k):
            vals = group_means_nan(e, s)
            j = int(np.argmax(vals))
            if vals[j] > best:
                best, best_set = float(vals[j]), s[j]
        if best_set is None:  # no subset has a calculable pair
            return dict(max=-np.inf, members=[], exact=True)
        return dict(max=best, members=[int(x) for x in best_set], exact=True)
    sub = np.nan_to_num(e[np.ix_(majority, majority)], nan=0.0)
    res = decoy_max(sub, np.arange(len(majority)), k, exhaustive_limit=0)
    members = [int(majority[x]) for x in res["members"]]
    val = float(group_means_nan(e, np.array([members]))[0])
    return dict(max=val, members=members, exact=False)
