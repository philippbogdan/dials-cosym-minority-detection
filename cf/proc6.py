"""Procedure candidate P6: noise-calibrated group scan + self-decoy population criterion.

As P5 (cf.proc5) for the scan and the counting-noise test (bootstrap over reflections).
Population criterion (the unsupervised counterpart of the pre-registered decoy): a candidate
group G is a population of its own when its coherence S(G) (mean residual of rij from one
population within G, correlation units) exceeds `ratio` times the strongest structure inside
the populations it is separated from: the most coherent subset of each of those populations,
of size min(|G|, |P| // 2) (exhaustive search when there are at most `exhaustive` subsets,
otherwise greedy growth with swaps). A rejected candidate is merged back into the population
it is closest to, and the search continues in the rest, so a minority is still found when the
strongest structure was a sub-population of the majority.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from cf.proc5 import Data, floor_p, scan


def _mean_within(e, members):
    k = len(members)
    if k < 2:
        return -np.inf
    blk = e[np.ix_(members, members)]
    iu = np.triu_indices(k, 1)
    v = blk[iu]
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else -np.inf


def _affinity(r, a, b):
    """Mean correlation between sets a and b over the geometric mean of their mean internal
    correlations (as dials.algorithms.correlation.groups._affinity)."""
    between = np.nanmean(r[np.ix_(a, b)])
    inside = [np.nanmean(r[np.ix_(x, x)][np.triu_indices(len(x), 1)]) if len(x) > 1 else np.nan
              for x in (a, b)]
    inside = [v for v in inside if np.isfinite(v) and v > 0]
    if not np.isfinite(between) or not inside:
        return -np.inf
    return float(between / np.sqrt(np.prod(inside) if len(inside) == 2 else inside[0] ** 2))


def best_subset(e, pool, k, exhaustive=200_000):
    """Most coherent size-k subset of `pool` (by mean residual within)."""
    pool = np.asarray(pool)
    if k < 2 or len(pool) < k:
        return -np.inf, None
    total = math.comb(len(pool), k)
    ez = np.nan_to_num(e, nan=0.0)
    if total <= exhaustive:
        best, best_set = -np.inf, None
        combos = np.array(list(itertools.combinations(range(len(pool)), k)))
        sub = ez[np.ix_(pool, pool)]
        for start in range(0, len(combos), 50_000):
            c = combos[start : start + 50_000]
            blk = sub[c[:, :, None], c[:, None, :]]
            vals = (blk.sum((1, 2)) - np.einsum("bii->b", blk)) / (k * (k - 1))
            j = int(np.argmax(vals))
            if vals[j] > best:
                best, best_set = float(vals[j]), pool[c[j]]
        return best, best_set
    sub = ez[np.ix_(pool, pool)]
    m = len(pool)
    best, best_set = -np.inf, None
    for a in range(m):
        cur = [a]
        while len(cur) < k:
            rest = [x for x in range(m) if x not in cur]
            cur.append(rest[int(np.argmax(sub[np.ix_(rest, cur)].sum(1)))])
        improved = True
        while improved:
            improved = False
            for pos in range(k):
                others = cur[:pos] + cur[pos + 1 :]
                cand = [x for x in range(m) if x not in cur]
                gains = sub[np.ix_(cand, others)].sum(1)
                if gains.max() > sub[cur[pos], others].sum() + 1e-15:
                    cur[pos] = cand[int(np.argmax(gains))]
                    improved = True
        val = _mean_within(sub, np.array(cur))
        if val > best:
            best, best_set = val, pool[np.array(cur)]
    return best, best_set


def decoy_refit(D, populations, k, exhaustive=200_000):
    """P6b: each population's structure on residuals of its own rank-1 fit."""
    from cf.oracle import rank1_fit

    best = -np.inf
    for p in populations:
        if len(p) < 4:
            continue
        sub = np.ix_(p, p)
        r, w = D.r[sub], D.w[sub]
        rho = rank1_fit(r, w)
        e = r - np.outer(rho, rho)
        e[w <= 0] = np.nan
        np.fill_diagonal(e, np.nan)
        val, _ = best_subset(e, np.arange(len(p)), min(k, len(p) // 2), exhaustive)
        best = max(best, val)
    return best


def decoy_full(e, populations, k, exhaustive=200_000):
    """P6 (first version): on the candidate's residuals, not centred."""
    best = -np.inf
    for p in populations:
        if len(p) < 4:
            continue
        val, _ = best_subset(e, p, min(k, len(p) // 2), exhaustive)
        best = max(best, val)
    return best


def decoy(e, populations, k, exhaustive=200_000):
    """Strongest structure inside the given populations (each of at least 4 crystals): the
    most coherent subset of size min(k, |P| // 2), on the same residuals as the candidate,
    measured from the population's own mean within-pair residual (the part a single
    population fit leaves in every pair of a population that is one side of a split)."""
    best = -np.inf
    for p in populations:
        if len(p) < 4:
            continue
        base = _mean_within(e, p)
        kd = min(k, len(p) // 2)
        val, _ = best_subset(e, p, kd, exhaustive)
        best = max(best, val - base)
    return best


def _gap(e, members, outsiders):
    k = len(members)
    within = _mean_within(e, members)
    if len(outsiders) == 0:
        return within
    outs = np.nanmean(e[np.ix_(members, outsiders)], axis=0)
    return float(within - np.nanmax(outs))


def gap_decoy(e, populations, k):
    """Largest gap of a subset of size <= min(k, |P| // 2) inside each population (greedy
    growth from every pair, the population's other members as outsiders)."""
    from cf.scan import grow_all

    best = -np.inf
    for p in populations:
        if len(p) < 4:
            continue
        kd = min(k, len(p) // 2)
        if kd < 2:
            continue
        sub = np.nan_to_num(e[np.ix_(p, p)], nan=0.0)
        order, within, outside = grow_all(sub, kd)
        g = within[:, 2 : kd + 1] - outside[:, 2 : kd + 1]
        best = max(best, float(np.nanmax(g)))
    return best


def usable(D: Data, min_pairs):
    """Crystals with at least `min_pairs` pairs that have a correlation and counting noise,
    dropping the others until every remaining crystal qualifies."""
    ok = (D.w > 0) & (np.abs(D.d).sum(0) > 0)
    np.fill_diagonal(ok, False)
    keep = np.ones(len(ok), bool)
    while True:
        cnt = (ok & keep[None, :] & keep[:, None]).sum(1)
        drop = keep & (cnt < min_pairs)
        if not drop.any():
            return np.where(keep)[0]
        keep &= ~drop


def analyze_usable(D: Data, min_pairs=3, **kw):
    """analyze() on the crystals with enough usable pairs; the others join the largest group
    (they are not judged)."""
    keep = usable(D, min_pairs)
    n = D.r.shape[0]
    if len(keep) == n:
        return analyze(D, **kw)
    if len(keep) < 2:
        return [np.arange(n)]
    groups = [keep[g] for g in analyze(D.sub(keep), **kw)]
    j = int(np.argmax([len(g) for g in groups]))
    groups[j] = np.sort(np.concatenate([groups[j], np.setdiff1d(np.arange(n), keep)]))
    return groups


def analyze(D: Data, alpha=0.01, ratio=3.0, min_pop=5, depth=0, trace=None, exhaustive=200_000,
            stat="gap", scan_fn=None, criterion="coherence", decoy_mode="centred", abs_min=None,
            calibrate=False, abs_q=0.0, abs_test=False, abs_alpha=None, abs_p=None,
            abs_scale=None, alpha_small=None, small_k=2):
    n = D.r.shape[0]
    idx = np.arange(n)
    if n < 2 + min_pop or depth > 8:
        return [idx]
    s = scan_fn(D) if scan_fn is not None else scan(D, stat=stat, calibrate=calibrate,
                                                    floor_null=abs_test)
    if s is None or s["p"] >= alpha:
        return [idx]
    G = s["members"]
    rest = np.setdiff1d(idx, G)
    if len(rest) < min_pop:
        return [idx]
    e = np.where(s["valid"], s["e"], np.nan)
    # "absolute" (frontier point p6abs) judges the gap; "coherence" and "either" the mean within
    sG = _gap(e, G, rest) if criterion in ("gap", "absolute") else _mean_within(e, G)
    rest_groups = [rest[g] for g in analyze(D.sub(rest), alpha, ratio, min_pop, depth + 1,
                                            None, exhaustive, stat, scan_fn, criterion, decoy_mode, abs_min,
                                            calibrate, abs_q, abs_test, abs_alpha, abs_p, abs_scale, alpha_small, small_k)]
    if criterion == "gap":
        d = gap_decoy(e, rest_groups, len(G))
    elif decoy_mode == "full":
        d = decoy_full(e, rest_groups, len(G), exhaustive)
    elif decoy_mode == "refit":
        d = decoy_refit(D, rest_groups, len(G), exhaustive)
    else:
        d = decoy(e, rest_groups, len(G), exhaustive)
    if criterion == "absolute":  # frontier point: a fixed coherence instead of the ratio
        ok = sG >= abs_min
    elif criterion == "either":  # separated from the rest's structure, or coherent beyond abs_min
        floor = abs_min * (2.0 / len(G) if abs_scale == "inv_k" else 1.0)
        se_sG = 0.0
        if abs_q and s.get("noise") is not None:  # beyond abs_min by abs_q noise SDs of S(G)
            iu = np.triu_indices(len(G), 1)
            blk = s["noise"][:, G][:, :, G][:, iu[0], iu[1]]
            vb = s["valid"][np.ix_(G, G)][iu]
            se_sG = float(np.std(blk[:, vb].mean(1))) if vb.any() else 0.0
        if abs_test:  # more coherent than abs_min beyond counting noise, look-elsewhere included
            floor_ok = floor_p(s, G, floor) < (alpha if abs_alpha is None else abs_alpha)
        else:
            floor_ok = sG - abs_q * se_sG >= floor
        if abs_p is not None:  # the floor route needs stronger evidence against counting noise
            floor_ok = floor_ok and s["p"] < abs_p
        ok = sG > 0 and (d <= 0 or sG > ratio * d or floor_ok)
    else:
        ok = sG > 0 and (d <= 0 or sG > ratio * d)
    if alpha_small is not None and len(G) <= small_k and s["p"] >= alpha_small:
        ok = False  # small groups need stronger evidence against counting noise
    if trace is not None:
        trace.append(dict(depth=depth, n=n, p=s["p"], k=len(G), S=sG, decoy=d, ok=ok))
    if ok:
        g_groups = ([G[g] for g in analyze(D.sub(G), alpha, ratio, min_pop, depth + 1, None,
                                           exhaustive, stat, scan_fn, criterion, decoy_mode, abs_min,
                                           calibrate, abs_q, abs_test, abs_alpha, abs_p, abs_scale, alpha_small, small_k)]
                    if len(G) >= 2 + min_pop else [G])
        return g_groups + rest_groups
    # not a population of its own: merge G into the closest population of the rest, the one
    # it correlates with as well as that population's members correlate with each other
    if len(rest_groups) == 1:
        return [idx]
    rr = np.where(s["valid"], D.r, np.nan)
    sims = [_affinity(rr, G, g) for g in rest_groups]
    j = int(np.argmax(sims))
    merged = [np.sort(np.concatenate([rest_groups[j], G])) if i == j else g
              for i, g in enumerate(rest_groups)]
    return merged


def detect(r, w, deltas, **kw):
    return analyze(Data(r, w, deltas), **kw)
