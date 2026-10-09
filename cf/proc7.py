"""Prototype (cryo, in-sample): a detector for the sparse regime, where v6's usable-pairs guard leaves
no data set to judge. On the pairs whose correlation has an analytic error (Fisher's z, the noise
v6 uses, at face value), the group with the largest inverse-variance statistic sum(z u) /
sqrt(sum(u^2)) (u = 1 / se) is searched by greedy growth over sizes 3 to n // 2, standardised per
size by the same search in the noise replicates; its look-elsewhere p-value (Gumbel tail of the
maximum over sizes) must be below `alpha`, and the group must pass v6's rules (separation from the
rest's own structure, or coherence beyond min_coherence * 2 / k)."""
from __future__ import annotations

import numpy as np

from cf.limits2 import gumbel_p
from cf.limits3 import Trial, grow_iv
from cf.proc6 import _mean_within, decoy


def sparse_groups(rij, wij, deltas, alpha=0.001, kmin=3, separation_ratio=3.5, min_coherence=0.05):
    n = rij.shape[0]
    kmax = n // 2
    if n < 6 or kmax < kmin:
        return None, None
    T = Trial(rij, wij, deltas, 1.0)
    if T.valid.sum() < 6:
        return None, None
    u2 = T.u ** 2

    def best_by_size(z):
        g = grow_iv(np.where(T.valid, z, 0.0) * T.u, u2, kmax)
        return g

    obs = best_by_size(T.z)
    sims = np.array([best_by_size(x)[:, kmin:].max(0) for x in T.zb])  # [B, sizes]
    sims = np.where(np.isfinite(sims), sims, np.nan)
    mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
    o = obs[:, kmin:]
    t = np.where(np.isfinite(o), (o - mu[None]) / sd[None], -np.inf)
    path, j = np.unravel_index(int(np.nanargmax(t)), t.shape)
    k = j + kmin
    p = gumbel_p(t[path, j], np.nanmax((sims - mu) / sd, 1))
    # the members of the best path at size k: replay greedy growth for that seed
    seeds = np.array(np.triu_indices(n, 1)).T
    x = np.where(T.valid, T.z, 0.0) * T.u
    members = [int(seeds[path][0]), int(seeds[path][1])]
    tot = x[members[0]] + x[members[1]]
    while len(members) < k:
        cand = np.where(np.isin(np.arange(n), members), -np.inf, tot)
        b = int(np.argmax(cand))
        members.append(b)
        tot = tot + x[b]
    members = np.array(sorted(members))
    rest = np.setdiff1d(np.arange(n), members)
    coherence = _mean_within(T.eS, members)
    d = decoy(T.eS, [rest], k)
    ok = p < alpha and (coherence >= min_coherence * 2 / k or (d > 0 and coherence / d >= separation_ratio)
                        or (d <= 0 and coherence > 0))
    log = dict(p=float(p), size=int(k), coherence=float(coherence), decoy=float(d), accepted=bool(ok))
    return (members if ok else None), log


def _local_search(z, u, valid, mu, sd, G, kmin, kmax, max_iter=50):
    """Improve G by single moves (add a crystal, remove one, or swap one in for one out) while the
    inverse-variance statistic standardised by its size's noise maxima rises."""
    n = z.shape[0]
    zz, uu = np.where(valid, z, 0.0), np.where(valid, u, 0.0)

    def t_of(members):
        k = len(members)
        if k < kmin or k > kmax:
            return -np.inf
        blk = np.ix_(members, members)
        iu = np.triu_indices(k, 1)
        den = (uu[blk][iu] ** 2).sum()
        if den <= 0:
            return -np.inf
        x = (zz[blk][iu] * uu[blk][iu]).sum() / np.sqrt(den)
        return (x - mu[k - kmin]) / sd[k - kmin]

    G = sorted(int(i) for i in G)
    best = t_of(np.array(G))
    for _ in range(max_iter):
        cands = []
        out_ = [i for i in range(n) if i not in G]
        cands += [sorted(G + [a]) for a in out_]
        cands += [[x for x in G if x != r] for r in G]
        cands += [sorted([x for x in G if x != r] + [a]) for r in G for a in out_]
        vals = [t_of(np.array(c)) for c in cands]
        j = int(np.argmax(vals))
        if vals[j] <= best + 1e-12:
            break
        best, G = vals[j], cands[j]
    return np.array(G), best


def sparse_groups_local(rij, wij, deltas, alpha=0.001, kmin=3, separation_ratio=3.5, min_coherence=0.05):
    """As sparse_groups, the greedy result improved by local moves, in the data and in every noise
    replicate alike (the look-elsewhere null of the same search)."""
    n = rij.shape[0]
    kmax = n // 2
    if n < 6 or kmax < kmin:
        return None, None
    T = Trial(rij, wij, deltas, 1.0)
    if T.valid.sum() < 6:
        return None, None
    u2 = T.u ** 2
    # per-size standardisation from the greedy maxima of the replicates
    sims = np.array([grow_iv(np.where(T.valid, x, 0.0) * T.u, u2, kmax)[:, kmin:].max(0) for x in T.zb])
    sims = np.where(np.isfinite(sims), sims, np.nan)
    mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
    seeds = np.array(np.triu_indices(n, 1)).T

    def search(z):
        g = grow_iv(np.where(T.valid, z, 0.0) * T.u, u2, kmax)[:, kmin:]
        t = np.where(np.isfinite(g), (g - mu[None]) / sd[None], -np.inf)
        path, j = np.unravel_index(int(np.argmax(t)), t.shape)
        k = j + kmin
        x = np.where(T.valid, z, 0.0) * T.u
        members = [int(seeds[path][0]), int(seeds[path][1])]
        tot = x[members[0]] + x[members[1]]
        while len(members) < k:
            cand = np.where(np.isin(np.arange(n), members), -np.inf, tot)
            b = int(np.argmax(cand))
            members.append(b)
            tot = tot + x[b]
        return _local_search(z, T.u, T.valid, mu, sd, members, kmin, kmax)

    members, t_obs = search(T.z)
    null = np.array([search(x)[1] for x in T.zb])
    p = gumbel_p(t_obs, null[np.isfinite(null)])
    rest = np.setdiff1d(np.arange(n), members)
    k = len(members)
    coherence = _mean_within(T.eS, members)
    d = decoy(T.eS, [rest], k)
    ok = p < alpha and (coherence >= min_coherence * 2 / k or (d > 0 and coherence / d >= separation_ratio)
                        or (d <= 0 and coherence > 0))
    return (members if ok else None), dict(p=float(p), size=int(k), coherence=float(coherence),
                                           decoy=float(d), accepted=bool(ok))
