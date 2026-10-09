"""Calibrated coherent-group detection (prototype of the proposed DIALS procedure).

Input: the matrix c of deviation correlations between crystals (cf.oracle) and the number of
reflections each pair shares. A group is a set of at least two crystals that deviate from
the rest in the same way; its coherence C(G) is the mean of c over its pairs.

Search: from every pair, a group is grown greedily (adding the crystal with the largest
summed c to the members), keeping, for each size, the most coherent group found.

Null: one population whose crystals differ continuously (crystal states drawn from a single
Gaussian, as unit-cell or hydration variation would make them) plus counting noise:
    c_ij = mu0 + b_i . b_j + e_ij,  b_i ~ N(0, Theta),  e_ij ~ N(0, kappa / n_ij),
with mu0, Theta (from the leading eigenvalues of c) and kappa (from the remaining scatter)
fitted to the data at hand. A minority inflates Theta, which makes the null conservative.

Test: for each size k the observed best coherence is ranked among those of B simulated null
matrices searched the same way; the smallest per-size p-value is itself calibrated against
the same simulations (min-p), and a group is reported when that global p-value is below
alpha. The search then continues on the remaining crystals and inside the group.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def best_groups(c: np.ndarray, kmax: int | None = None):
    """Greedy growth from every pair. Returns (best C per size k=2..kmax, members)."""
    n = c.shape[0]
    kmax = n - 2 if kmax is None else min(kmax, n - 2)
    best = np.full(kmax + 1, -np.inf)
    members: list = [None] * (kmax + 1)
    iu, ju = np.triu_indices(n, 1)
    for a, b in zip(iu, ju):
        inset = np.zeros(n, bool)
        inset[[a, b]] = True
        tot = c[a] + c[b]  # summed c of every crystal to the current members
        s = c[a, b]  # sum over member pairs
        k = 2
        while True:
            val = s / (k * (k - 1) / 2)
            if val > best[k]:
                best[k] = val
                members[k] = np.flatnonzero(inset)
            if k == kmax:
                break
            cand = np.where(inset, -np.inf, tot)
            x = int(np.argmax(cand))
            s += tot[x]
            inset[x] = True
            tot = tot + c[x]
            k += 1
    return best, members


@dataclass
class NullModel:
    mu0: float
    loadings_var: np.ndarray  # Theta's eigenvalues (per state dimension)
    kappa: float  # noise scale: var(e_ij) = kappa / n_ij


def fit_null(c: np.ndarray, npairs: np.ndarray, max_dims: int = 10) -> NullModel:
    n = c.shape[0]
    iu = np.triu_indices(n, 1)
    mu0 = float(np.median(c[iu]))
    a = c - mu0
    np.fill_diagonal(a, 0.0)
    inv_n = np.where(npairs > 0, 1.0 / np.maximum(npairs, 1), 0.0)
    # noise level from the scatter of c around its low-rank part: start from all of it
    ev, evec = np.linalg.eigh(a)
    order = np.argsort(ev)[::-1]
    ev, evec = ev[order], evec[:, order]
    resid = a.copy()
    kappa = float(np.median(resid[iu] ** 2 / np.maximum(inv_n[iu], 1e-12)))
    dims = 0
    for _ in range(5):
        edge = 2.0 * np.sqrt(n * kappa * np.mean(inv_n[iu]))
        dims = int(min(max_dims, np.sum(ev > edge)))
        low = (evec[:, :dims] * ev[:dims]) @ evec[:, :dims].T
        resid = a - low
        np.fill_diagonal(resid, 0.0)
        kappa = float(np.mean(resid[iu] ** 2 / np.maximum(inv_n[iu], 1e-12)))
    theta = np.maximum(ev[:dims], 0.0) / n
    return NullModel(mu0=mu0, loadings_var=theta, kappa=kappa)


def simulate(model: NullModel, npairs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = npairs.shape[0]
    L = len(model.loadings_var)
    out = np.full((n, n), model.mu0)
    if L:
        b = rng.standard_normal((n, L)) * np.sqrt(model.loadings_var)[None, :]
        out = out + b @ b.T
    e = rng.standard_normal((n, n)) * np.sqrt(model.kappa / np.maximum(npairs, 1))
    e = np.triu(e, 1)
    out = out + e + e.T
    np.fill_diagonal(out, 0.0)
    return out


@dataclass
class Test:
    p_global: float
    k_best: int
    members: np.ndarray
    p_by_size: np.ndarray
    observed: np.ndarray
    model: NullModel = field(repr=False)


def test_once(c, npairs, n_sim=200, seed=0, kmax=None) -> Test:
    n = c.shape[0]
    kmax = n - 2 if kmax is None else min(kmax, n - 2)
    obs, mem = best_groups(c, kmax)
    model = fit_null(c, npairs)
    rng = np.random.default_rng(seed)
    sims = np.array([best_groups(simulate(model, npairs, rng), kmax)[0] for _ in range(n_sim)])
    ks = np.arange(2, kmax + 1)
    # per-size p-values of the observation and of every simulation (rank among the others)
    p_obs = (1 + (sims[:, ks] >= obs[ks][None, :]).sum(0)) / (1 + n_sim)
    ranks = sims[:, ks].argsort(0).argsort(0)  # 0 = smallest
    p_sim = (n_sim - ranks) / (1 + n_sim)  # p of each simulation among the others (approx.)
    minp_obs = p_obs.min()
    minp_sim = p_sim.min(1)
    p_global = (1 + np.sum(minp_sim <= minp_obs)) / (1 + n_sim)
    kb = int(ks[np.argmin(p_obs)])
    return Test(p_global=float(p_global), k_best=kb, members=mem[kb], p_by_size=p_obs,
                observed=obs, model=model)


def detect(c, npairs, alpha=0.01, n_sim=200, seed=0, min_size=2, _depth=0):
    """Recursive detection. Returns a list of groups (arrays of crystal indices)."""
    n = c.shape[0]
    idx = np.arange(n)
    if n < 2 * min_size or _depth > 6:
        return [idx], []
    t = test_once(c, npairs, n_sim=n_sim, seed=seed)
    if t.p_global >= alpha:
        return [idx], [t]
    g = np.asarray(t.members)
    rest = np.setdiff1d(idx, g)
    tests = [t]
    groups = []
    for part in (g, rest):
        if len(part) >= 2 * min_size:
            sub, subtests = detect(c[np.ix_(part, part)], npairs[np.ix_(part, part)],
                                   alpha, n_sim, seed + 1 + _depth, min_size, _depth + 1)
            groups.extend(part[s] for s in sub)
            tests.extend(subtests)
        else:
            groups.append(part)
    return groups, tests


def gap_scan(c: np.ndarray, npairs: np.ndarray, kappa: float, kmax: int | None = None):
    """Greedy growth from every pair; for every group met, the standardized gap
        Z(G) = [C(G) - max_{x not in G} mean_{i in G} c_ix] / (s * sqrt(2/(k(k-1)) + 1/k)),
    with s the counting-noise SD of one c (sqrt(kappa / median n_ij)). Returns the best Z,
    its members, and the best Z per size."""
    n = c.shape[0]
    kmax = n - 2 if kmax is None else min(kmax, n - 2)
    iu = np.triu_indices(n, 1)
    s_e = np.sqrt(kappa / max(1.0, float(np.median(npairs[iu]))))
    best = -np.inf
    best_members = None
    best_k = np.full(kmax + 1, -np.inf)
    for a, b in zip(*iu):
        inset = np.zeros(n, bool)
        inset[[a, b]] = True
        tot = c[a] + c[b]
        ssum = c[a, b]
        k = 2
        while True:
            within = ssum / (k * (k - 1) / 2)
            outside = np.where(inset, -np.inf, tot / k)
            z = (within - outside.max()) / (s_e * np.sqrt(2.0 / (k * (k - 1)) + 1.0 / k))
            if z > best_k[k]:
                best_k[k] = z
            if z > best:
                best, best_members = z, np.flatnonzero(inset)
            if k == kmax:
                break
            x = int(np.argmax(np.where(inset, -np.inf, tot)))
            ssum += tot[x]
            inset[x] = True
            tot = tot + c[x]
            k += 1
    return float(best), best_members, best_k


def gap_test(c, npairs, n_sim=200, seed=0, kmax=None):
    model = fit_null(c, npairs)
    z, mem, zk = gap_scan(c, npairs, model.kappa, kmax)
    rng = np.random.default_rng(seed)
    null = np.empty(n_sim)
    for b in range(n_sim):
        cs = simulate(model, npairs, rng)
        null[b] = gap_scan(cs, npairs, fit_null(cs, npairs).kappa, kmax)[0]
    p = (1 + np.sum(null >= z)) / (1 + n_sim)
    return dict(p=float(p), z=z, members=mem, z_by_size=zk, null=null, model=model)
