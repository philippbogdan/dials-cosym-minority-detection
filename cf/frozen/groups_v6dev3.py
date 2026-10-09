"""Calibrated detection of groups of data sets in a pairwise correlation matrix.

A group is a set of at least two data sets that agree with each other better than one
population of data sets would. It is reported when

1. it is significant against counting noise: the noise of every pairwise correlation is
   the large-sample error of Fisher's z, 1 / sqrt(n_eff - 3), with n_eff the effective number
   of reflections of the pair's sigma-weighted correlation (the quantity cosym keeps as its
   weight; the same error underlies its ``weights=standard_error``), drawn independently for
   every pair (``noise_model=analytic``); or by a bootstrap over reflections (the unique
   reflections are redrawn with Poisson weights, the same draw for all data sets, and every
   correlation recomputed as cosym computes it; ``noise_model=bootstrap``, which overstates
   the noise where pairs share few reflections). The scan below is calibrated against the
   same scan on the noise-only replicates; a pair of data sets, the most numerous and the
   most exposed kind of candidate, needs ``pair_significance``; and
2. it is a population of its own: it stands out from the populations it is separated from
   by more than their own internal structure (data sets of one population differ too: cell,
   hydration, dose; and real sub-populations exist inside one population), that is its
   coherence exceeds ``separation_ratio`` times the coherence of the most coherent subset of
   the same size inside each of those populations (its "decoy"); or its coherence is at
   least ``min_coherence`` x 2 / k for a group of k data sets (correlation units), a
   difference too large for the internal structure of one population, which keeps a group
   visible when the population it is separated from has strong internal structure of its
   own. (The strongest structure inside one population falls about as 1 / k with the size
   of the subset, as does a group's coherence: the rank-1 fit absorbs more of a larger
   group.)

Coherence: residuals of the correlation matrix from a single population (rank-1 fit, as
cosym in one dimension); the coherence of a set is the mean residual over its pairs; a
decoy's is measured from its population's mean.

Data sets with fewer than ``min_usable_pairs`` correlations that have an error (enough
shared reflections) are not judged: they stay with the largest group.

Scan: residuals standardised by their counting-noise standard errors; from every pair of data sets
a group is grown greedily, keeping for every size the best gap (mean within the group minus
the best outsider's mean to the group); each size is standardised by the noise replicates,
the maximum over sizes is the statistic, and its p-value comes from a Gumbel distribution
fitted to the replicates' maxima.

The rest is analysed first (it may hold populations of its own, whose decoys are then
searched separately); a candidate that fails the population test is merged back into the
population it is closest to and the search continues, so a minority is still found when the
strongest structure in the data was a sub-population of the majority.
"""

from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------------------
# Pairwise correlations on a dense layout (as dials.algorithms.symmetry.cosym.target)
# ---------------------------------------------------------------------------------------


def dense_from_target(target):
    """Data sets x unique reflections arrays from a cosym Target (reflections paired by the
    cosym matcher only: epsilon == 1 in the Patterson group)."""
    data = target._data
    pg = target._patterson_group
    hkl = np.array(data.indices(), dtype=np.int64).reshape(-1, 3)
    eps1 = np.array(pg.epsilon(data.indices()), dtype=np.int64) == 1
    ids = np.asarray(target._lattice_ids)
    intensities = np.array(data.data(), dtype=np.float64)
    sigmas = np.array(data.sigmas(), dtype=np.float64)
    n = len(target._lattices)
    hkl, intensities, sigmas, ids = hkl[eps1], intensities[eps1], sigmas[eps1], ids[eps1]
    key = (hkl[:, 0] * 4096 + hkl[:, 1]) * 4096 + hkl[:, 2]
    _, col = np.unique(key, return_inverse=True)
    nref = int(col.max()) + 1 if len(col) else 0
    y = np.zeros((n, nref))
    v = np.zeros((n, nref))
    m = np.zeros((n, nref), dtype=bool)
    y[ids, col] = intensities
    v[ids, col] = sigmas**2
    m[ids, col] = True
    return y, v, m


def weighted_cc_matrix(y, v, m, weights=None, min_pairs=3):
    """Sigma-weighted correlation coefficients for all pairs, as cosym computes them
    (ExtendedDatasetStatistics.weighted_cchalf), with optional per-reflection weights."""
    n = y.shape[0]
    wt = np.ones(y.shape[1]) if weights is None else np.asarray(weights, dtype=float)
    r = np.zeros((n, n))
    for i in range(n - 1):
        mm = m[i][None, :] & m[i + 1 :] & (wt[None, :] > 0)
        cnt = mm.sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(mm, wt[None, :] / (v[i][None, :] + v[i + 1 :]), 0.0)
            nw = w / w.sum(1)[:, None]
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


def bootstrap_replicates(y, v, m, n_bootstrap=100, seed=0):
    """Correlation matrices recomputed with the reflections redrawn (Poisson weights, the
    same for all data sets), minus the unweighted matrix: counting noise only."""
    rng = np.random.default_rng(seed)
    r0 = weighted_cc_matrix(y, v, m)
    out = np.empty((n_bootstrap,) + r0.shape)
    for b in range(n_bootstrap):
        out[b] = weighted_cc_matrix(y, v, m, rng.poisson(1.0, y.shape[1]))
    return out - r0[None]


def effective_sizes(y, v, m):
    """Kish effective number of reflections of every pair's sigma-weighted correlation
    (weights 1 / (v_i + v_j) over the reflections both data sets have), as cosym's wij
    before it rescales them."""
    n = y.shape[0]
    out = np.zeros((n, n))
    for i in range(n - 1):
        mm = m[i][None, :] & m[i + 1 :]
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(mm, 1.0 / (v[i][None, :] + v[i + 1 :]), 0.0)
        s1, s2 = w.sum(1), (w * w).sum(1)
        out[i, i + 1 :] = np.where(s2 > 0, s1**2 / np.where(s2 > 0, s2, 1.0), 0.0)
    return out + out.T


def analytic_replicates(r, n_eff, n_replicates=100, seed=0, min_n_eff=4.0):
    """Correlation matrices with counting noise drawn from the large-sample error of Fisher's
    z, 1 / sqrt(n_eff - 3), independently for every pair, minus the matrix. Pairs with
    n_eff <= min_n_eff get no noise (and are left out by the scan)."""
    n = r.shape[0]
    r = np.array(r, dtype=np.float64)
    np.fill_diagonal(r, 0.0)
    ok = (n_eff > min_n_eff) & (r != 0)
    np.fill_diagonal(ok, False)
    sez = np.where(ok, 1.0 / np.sqrt(np.where(ok, n_eff - 3.0, 1.0)), 0.0)
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(seed)
    out = np.zeros((n_replicates, n, n))
    for b in range(n_replicates):
        eps = np.zeros((n, n))
        eps[iu] = rng.standard_normal(len(iu[0])) * sez[iu]
        eps = eps + eps.T
        out[b] = np.where(ok, np.tanh(z + eps) - r, 0.0)
    return out


def noise_replicates(y, v, m, r, model="analytic", n_replicates=100, seed=0):
    """Counting-noise replicates of the correlation matrix r (replicate - r)."""
    if model == "bootstrap":
        return bootstrap_replicates(y, v, m, n_bootstrap=n_replicates, seed=seed)
    return analytic_replicates(r, effective_sizes(y, v, m), n_replicates, seed)


# ---------------------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------------------


def rank1_fit(r, w, iters=5000, tol=1e-12):
    """Weighted least-squares rank-1 fit (cosym in one dimension), diagonal excluded. Each
    step moves halfway to every rho_i's least-squares value given the others (the plain
    simultaneous update oscillates when the rho_i differ widely)."""
    w = np.array(w, dtype=float)
    np.fill_diagonal(w, 0.0)
    wsum = w.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.sqrt(np.clip(np.where(wsum > 0, (w * r).sum(1) / wsum, 0.0), 1e-6, None))
        for _ in range(iters):
            num = (w * r * rho[None, :]).sum(1)
            den = (w * rho[None, :] ** 2).sum(1)
            new = 0.5 * (np.where(den > 0, num / den, 0.0) + rho)
            if np.max(np.abs(new - rho)) < tol:
                return new
            rho = new
    return rho


def _residuals(r, w):
    rho = rank1_fit(r, w)
    e = r - np.outer(rho, rho)
    np.fill_diagonal(e, 0.0)
    return e, rho


def greedy_growth(x, kmax, counts=None):
    """Greedy growth of every pair of data sets on the symmetric matrix x.

    Returns the order in which members were added (pairs x kmax), and per size k the mean of
    x within the group and the best outsider's mean of x to the group ([pairs, kmax + 1]);
    with `counts` (1 where an element of x is usable, else 0), also the standardised sum of x
    over the group's usable pairs (sum / sqrt(number))."""
    n = x.shape[0]
    seeds = np.array(np.triu_indices(n, 1)).T
    s = len(seeds)
    rows = np.arange(s)
    order = np.zeros((s, kmax), dtype=np.int64)
    order[:, 0], order[:, 1] = seeds[:, 0], seeds[:, 1]
    inset = np.zeros((s, n), dtype=bool)
    inset[rows, seeds[:, 0]] = True
    inset[rows, seeds[:, 1]] = True
    tot = x[seeds[:, 0]] + x[seeds[:, 1]]
    within = np.full((s, kmax + 1), np.nan)
    outside = np.full((s, kmax + 1), np.nan)
    stouffer = np.full((s, kmax + 1), np.nan)
    ssum = x[seeds[:, 0], seeds[:, 1]].copy()
    if counts is not None:
        ctot = counts[seeds[:, 0]] + counts[seeds[:, 1]]
        csum = counts[seeds[:, 0], seeds[:, 1]].astype(float)
    for k in range(2, kmax + 1):
        within[:, k] = ssum / (k * (k - 1) / 2)
        if counts is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                stouffer[:, k] = np.where(csum > 0, ssum / np.sqrt(np.maximum(csum, 1.0)), -np.inf)
        cand = np.where(inset, -np.inf, tot)
        best = cand.argmax(1)
        outside[:, k] = cand[rows, best] / k
        if k == kmax:
            break
        ssum = ssum + tot[rows, best]
        if counts is not None:
            csum = csum + ctot[rows, best]
            ctot = ctot + counts[best]
        inset[rows, best] = True
        tot = tot + x[best]
        order[:, k] = best
    if counts is not None:
        return order, within, outside, stouffer
    return order, within, outside


@dataclass
class Scan:
    p: float
    statistic: float
    members: np.ndarray
    residuals: np.ndarray = field(repr=False)
    rho: np.ndarray = field(repr=False)
    se: np.ndarray = field(repr=False)
    valid: np.ndarray = field(repr=False)


def scan(r, w, deltas, statistic="both"):
    """The most group-like set of data sets and its p-value against counting noise.

    statistic: "both" (default): the gap (mean standardised residual within the group minus
    the best outsider's mean to it) and the standardised sum of the residuals over the group's
    usable pairs, each standardised per size by the noise replicates; whichever stands out
    more counts, against the replicates' maxima over both (the gap misses a group with an
    outsider close to it, the sum alone a group inside a larger coherent block); "gap": the
    gap alone."""
    n = r.shape[0]
    kmax = n // 2
    if n < 4 or kmax < 2:
        return None
    valid = np.asarray(w) > 0
    np.fill_diagonal(valid, False)
    e, rho = _residuals(r, w)
    eb = np.empty_like(deltas)
    for b, d in enumerate(deltas):
        rb = r + d
        np.fill_diagonal(rb, 0.0)
        valid &= ~((rb == 0) & (np.asarray(w) > 0))
        eb[b] = _residuals(rb, w)[0]
    se = eb.std(0)
    valid &= se > 0
    # a correlation without noise of its own (too few reflections for an error) is left
    # out, not given the small spread of the refit
    valid &= np.any(deltas != 0, axis=0)
    if valid.sum() < 6:
        return None
    sse = np.where(valid, se, 1.0)

    counts = valid.astype(float)
    if statistic == "select":
        # choose the candidate of every size by the gap, test it by its standardised block sum;
        # the null applies the same rule to the noise replicates
        def chosen(z):
            order, within, outside, stouffer = greedy_growth(np.where(valid, z, 0.0), kmax, counts)
            gap = within - outside
            best = np.nanargmax(np.where(np.isfinite(gap[:, 2:]), gap[:, 2:], -np.inf), axis=0)
            vals = stouffer[best, np.arange(2, kmax + 1)]
            return np.where(np.isfinite(vals), vals, -np.inf), order, best

        obs, order, best = chosen(e / sse)
        sims = np.array([chosen((x - e) / sse)[0] for x in eb])
        sims = np.where(np.isfinite(sims), sims, np.nan)
        mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
        t_obs = np.where(np.isfinite(obs), (obs - mu) / sd, -np.inf)
        t_sim = np.nan_to_num((sims - mu) / sd, nan=-np.inf).max(1)
        t_sim = t_sim[np.isfinite(t_sim)]
        beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
        loc = t_sim.mean() - 0.5772156649 * beta
        kk = int(np.argmax(t_obs))
        t = float(t_obs[kk])
        p = float(-np.expm1(-np.exp(-(t - loc) / beta)))
        k = kk + 2
        members = np.sort(order[int(best[kk]), :k])
        return Scan(p=p, statistic=t, members=members, residuals=e, rho=rho, se=se, valid=valid)
    kinds = {"both": ("gap", "sum"), "gap": ("gap",), "sum": ("sum",)}[statistic]

    def stats(z):
        order, within, outside, stouffer = greedy_growth(np.where(valid, z, 0.0), kmax, counts)
        return {"gap": within - outside, "sum": stouffer}, order

    vals, order = stats(e / sse)
    sim_vals = [stats((x - e) / sse)[0] for x in eb]
    t_obs, t_sim = [], []
    for kind in kinds:
        obs = np.nanmax(vals[kind][:, 2:], axis=0)
        sims = np.array([np.nanmax(sv[kind][:, 2:], axis=0) for sv in sim_vals])
        sims = np.where(np.isfinite(sims), sims, np.nan)
        mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
        t_obs.append(np.where(np.isfinite(obs), (obs - mu) / sd, -np.inf))
        t_sim.append(np.nan_to_num((sims - mu) / sd, nan=-np.inf).max(1))
    t_obs = np.nan_to_num(np.array(t_obs), nan=-np.inf)  # [kinds, sizes]
    t_sim = np.max(np.array(t_sim), axis=0)
    t_sim = t_sim[np.isfinite(t_sim)]
    beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
    loc = t_sim.mean() - 0.5772156649 * beta
    which, kk = np.unravel_index(int(np.argmax(t_obs)), t_obs.shape)
    t = float(t_obs[which, kk])
    p = float(-np.expm1(-np.exp(-(t - loc) / beta)))
    k = int(kk) + 2
    members = np.sort(order[int(np.nanargmax(vals[kinds[which]][:, k])), :k])
    return Scan(p=p, statistic=t, members=members, residuals=e, rho=rho, se=se, valid=valid)


def _mean_within(e, members):
    k = len(members)
    if k < 2:
        return -np.inf
    v = e[np.ix_(members, members)][np.triu_indices(k, 1)]
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else -np.inf


def _affinity(r, a, b):
    """Mean correlation between sets a and b over the geometric mean of their mean internal
    correlations (data quality cancels)."""
    between = np.nanmean(r[np.ix_(a, b)])
    inside = [np.nanmean(r[np.ix_(x, x)][np.triu_indices(len(x), 1)]) if len(x) > 1 else np.nan
              for x in (a, b)]
    inside = [v for v in inside if np.isfinite(v) and v > 0]
    if not np.isfinite(between) or not inside:
        return -np.inf
    return float(between / np.sqrt(np.prod(inside) if len(inside) == 2 else inside[0] ** 2))


def most_coherent_subset(e, pool, k, exhaustive=200_000):
    """Mean residual within the most coherent size-k subset of `pool` (exhaustive search
    when there are at most `exhaustive` subsets, otherwise greedy growth with swaps)."""
    pool = np.asarray(pool)
    if k < 2 or len(pool) < k:
        return -np.inf
    ez = np.nan_to_num(e[np.ix_(pool, pool)], nan=0.0)
    m = len(pool)
    if math.comb(m, k) <= exhaustive:
        best = -np.inf
        combos = np.array(list(itertools.combinations(range(m), k)))
        for start in range(0, len(combos), 50_000):
            c = combos[start : start + 50_000]
            blk = ez[c[:, :, None], c[:, None, :]]
            vals = (blk.sum((1, 2)) - np.einsum("bii->b", blk)) / (k * (k - 1))
            best = max(best, float(vals.max()))
        return best
    best = -np.inf
    for a in range(m):
        cur = [a]
        while len(cur) < k:
            rest = [x for x in range(m) if x not in cur]
            cur.append(rest[int(np.argmax(ez[np.ix_(rest, cur)].sum(1)))])
        improved = True
        while improved:
            improved = False
            for pos in range(k):
                others = cur[:pos] + cur[pos + 1 :]
                cand = [x for x in range(m) if x not in cur]
                gains = ez[np.ix_(cand, others)].sum(1)
                if gains.max() > ez[cur[pos], others].sum() + 1e-15:
                    cur[pos] = cand[int(np.argmax(gains))]
                    improved = True
        best = max(best, _mean_within(e[np.ix_(pool, pool)], np.array(cur)))
    return best


def decoy(e, populations, k, exhaustive=200_000):
    """Strongest structure inside the given populations (of at least 4 data sets each): the
    most coherent subset of size min(k, |P| // 2), measured from P's mean within-pair
    residual."""
    best = -np.inf
    for p in populations:
        if len(p) < 4:
            continue
        base = _mean_within(e, p)
        val = most_coherent_subset(e, p, min(k, len(p) // 2), exhaustive)
        best = max(best, val - base)
    return best


@dataclass
class GroupDetection:
    """Find the groups of data sets in a pairwise correlation matrix.

    Args:
        rij: pairwise correlation matrix (cosym's rij).
        wij: its weights (cosym's wij); zero where a correlation could not be calculated.
        deltas: bootstrap replicates of rij minus rij (counting noise), shape (B, n, n).
        significance: largest p-value against counting noise for a group to be reported.
        separation_ratio: how many times the coherence of the strongest structure inside the
            populations it is separated from a group's coherence must exceed.
        min_population: smallest number of data sets in which internal structure can be
            judged; a group is only split from a rest at least this large.
        min_coherence: coherence (mean residual within the group, correlation units) at or
            above which a significant pair of data sets is a population of its own whatever
            the structure of the rest; for a group of k data sets the threshold is
            min_coherence x 2 / k. None to judge by separation_ratio alone.
        pair_significance: largest p-value for a group of two data sets (None: as any
            group).
        min_usable_pairs: data sets with fewer correlations that have counting noise are not
            judged (they stay with the largest group).
        statistic: the scan's statistic: "both" (gap or standardised sum over usable pairs,
            whichever stands out more) or "gap".
    """

    rij: np.ndarray
    wij: np.ndarray
    deltas: np.ndarray
    significance: float = 0.01
    separation_ratio: float = 3.5
    min_population: int = 5
    min_coherence: float | None = 0.05
    pair_significance: float | None = 0.001
    min_usable_pairs: int = 3
    statistic: str = "both"
    max_depth: int = 8
    exhaustive: int = 200_000
    log: list = field(default_factory=list)

    def run(self) -> list[np.ndarray]:
        n = self.rij.shape[0]
        keep = self._usable()
        if len(keep) < 2:
            return [np.arange(n)]
        groups = self._analyse(keep, 0)
        rest = np.setdiff1d(np.arange(n), keep)
        if len(rest):
            j = int(np.argmax([len(g) for g in groups]))
            groups[j] = np.sort(np.concatenate([groups[j], rest]))
        return groups

    def _usable(self):
        """Data sets with at least min_usable_pairs correlations that have counting noise,
        dropping the others until every remaining one qualifies."""
        ok = (np.asarray(self.wij) > 0) & np.any(self.deltas != 0, axis=0)
        np.fill_diagonal(ok, False)
        keep = np.ones(len(ok), dtype=bool)
        while True:
            cnt = (ok & keep[None, :] & keep[:, None]).sum(1)
            drop = keep & (cnt < self.min_usable_pairs)
            if not drop.any():
                return np.where(keep)[0]
            keep &= ~drop

    def _analyse(self, idx, depth):
        """Groups of the (sorted) original indices `idx`, as arrays of original indices."""
        n = len(idx)
        if n < 2 + self.min_population or depth > self.max_depth:
            return [idx]
        sub = np.ix_(idx, idx)
        s = scan(self.rij[sub], self.wij[sub], self.deltas[:, idx][:, :, idx], self.statistic)
        if s is None or s.p >= self.significance:
            return [idx]
        members = s.members  # positions within idx
        rest = np.setdiff1d(np.arange(n), members)
        if len(rest) < self.min_population:
            return [idx]
        e = np.where(s.valid, s.residuals, np.nan)
        coherence = _mean_within(e, members)
        rest_groups = self._analyse(idx[rest], depth + 1)
        local = [np.searchsorted(idx, g) for g in rest_groups]
        strongest = decoy(e, local, len(members), self.exhaustive)
        ok = coherence > 0 and (
            strongest <= 0
            or coherence > self.separation_ratio * strongest
            or (
                self.min_coherence is not None
                and coherence >= self.min_coherence * 2.0 / len(members)
            )
        )
        if (
            ok
            and len(members) == 2
            and self.pair_significance is not None
            and s.p >= self.pair_significance
        ):
            ok = False
        self.log.append(dict(depth=depth, n=n, p=s.p, size=len(members), coherence=coherence,
                             decoy=strongest, accepted=bool(ok)))
        if ok:
            if len(members) >= 2 + self.min_population:
                member_groups = self._analyse(idx[members], depth + 1)
            else:
                member_groups = [idx[members]]
            return member_groups + rest_groups
        if len(rest_groups) == 1:
            return [idx]
        # not a population of its own: back into the population it is closest to, i.e. the
        # one it correlates with as well as that population's members correlate with each
        # other (residuals would mislead here: the two ends of one varying population have
        # negative residuals with each other)
        r = np.where(s.valid, self.rij[sub], np.nan)
        closeness = [_affinity(r, members, g) for g in local]
        j = int(np.argmax(closeness))
        rest_groups[j] = np.sort(np.concatenate([rest_groups[j], idx[members]]))
        return rest_groups


def take_out_groups(labels, groups):
    """Cluster labels after taking the detected groups out of the OPTICS clusters: every
    group other than the largest becomes a cluster of its own (the largest, the background
    population, keeps the OPTICS labels); when no group holds more than half of the data
    sets, the detector has split the data into several populations and every group becomes
    a cluster. Clusters left with one data set become outliers (-1); labels are renumbered
    from 0."""
    labels = np.array(labels, dtype=int)
    groups = sorted((np.asarray(g) for g in groups), key=len, reverse=True)
    if len(groups) > 1:
        carved = groups[1:] if 2 * len(groups[0]) > len(labels) else groups
        next_label = labels.max() + 1 if len(labels) else 0
        for g in carved:
            labels[g] = next_label
            next_label += 1
    new = np.full_like(labels, -1)
    k = 0
    for lab in sorted(set(labels.tolist()) - {-1}):
        members = labels == lab
        if members.sum() >= 2:
            new[members] = k
            k += 1
    return new
