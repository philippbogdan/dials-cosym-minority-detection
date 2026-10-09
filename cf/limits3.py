"""Refined limits (PREREGISTRATION.md amendment 4): what the data allow an unsupervised
procedure to find, per mixture trial, with calibrated counting noise.

On the rank-1 residuals e of DIALS's rij (as the procedures use them):
1. Counting noise. Each residual is standardised by its counting-noise standard error: the
   analytic model (cf.bootstrap.analytic_replicates: Fisher's z with DIALS's effective number
   of reflections, refitted per replicate), multiplied by the factor measured against repeat
   measurements of the same crystals in disjoint wedge segments (cf.noisecal, observed /
   predicted, mean over species; 1 where the sweep has no second segment). A group's score is
   the standardised sum of its pairs, Z(G) = sum z_ij / sqrt(#pairs), the most powerful
   statistic against independent noise; the null is the largest Z over all size-k subsets of
   the trial in the noise replicates (Gumbel tail). A minority with p >= alpha lies within the
   counting noise: no procedure keeps false groups from noise below alpha and reports it.
2. Structure inside one population. For every single-population trial of the same collection
   and wedge and every size k: the most coherent size-k subset G' (largest mean residual S),
   its separation ratio (S over the most coherent subset of the rest, measured from the rest's
   own mean, as the procedures' decoy) and its counting-noise p-value as in 1. A minority M is
   within the populations' own structure if some G' of the same size that stands out of the
   counting noise (p < alpha) is at least as coherent (S) and at least as separated (ratio) as
   M is in its trial: any procedure that decides on coherence and separation and reported M
   would report G' in single-population data. Only groups whose rest has at least 5 crystals
   count (the smallest population the procedures judge structure in).
A miss that is neither is attributed to the procedure.

Writes results/limits/<design>_refined.jsonl (one line per mixture trial with k >= 2, and one
per single-population comparator group).
"""

from __future__ import annotations

from cf.result_io import open_result, glob_results

import itertools
import json
import math
import os
import sys

import numpy as np

from cf.analyze import load_design
from cf.bootstrap import load as bload
from cf.cmatrices import load
from cf.limits2 import gumbel_p
from cf.paths import RESULTS, XCAL, base_design, matrix_design, species_of
from cf.proc5 import Data, model
from cf.proc6 import _mean_within, best_subset, decoy
from cf.scan import grow_all

ALPHA = 0.01  # level recorded in the per-trial flags
LIMIT_ALPHA = 0.001  # level of the limit tests (amendment 13)
CAL_BY_SPECIES = os.environ.get("CF_CAL_BY_SPECIES", "1") == "1"  # amendment 16
MIN_REST = 5
TOL = 2.0  # standard errors (amendment 6)
N_COMPARATORS = 20  # most coherent subsets per single-population trial and size


def calibration(design):
    if design.endswith(XCAL):  # one factor learned on the other collection (calibration-transfer check)
        other = {"cryo": "rt", "rt": "cryo"}[base_design(design)]
        f = RESULTS / "noisecal" / f"{other}.json"
        vals = [v["observed_over_analytic"] for v in json.load(open_result(f)).values()
                if v.get("observed_over_analytic") is not None]
        mine = json.load(open_result(RESULTS / "noisecal" / f"{base_design(design)}.json"))
        return {float(k.split("|")[0]): float(np.mean(vals)) for k in mine}
    f = RESULTS / "noisecal" / f"{design}.json"
    if not f.exists():
        return {}
    d = json.load(open_result(f))
    by = {}
    for key, v in d.items():
        w = float(key.split("|")[0])
        if v.get("observed_over_analytic") is not None:
            by.setdefault(w, []).append(v["observed_over_analytic"])
    return {w: float(np.mean(v)) for w, v in by.items()}


def calibration_by_species(design):
    """{wedge: {species: observed / analytic}} from the repeat segments (amendment 16)."""
    if design.endswith(XCAL):  # a procedure does not know the species
        return {}
    f = RESULTS / "noisecal" / f"{design}.json"
    if not f.exists():
        return {}
    out = {}
    for key, v in json.load(open_result(f)).items():
        w, sp = key.split("|")
        if v.get("observed_over_analytic") is not None:
            out.setdefault(float(w), {})[sp] = float(v["observed_over_analytic"])
    return out


def pair_factors(names, species, by_species, mean_factor):
    """Per pair, the geometric mean of the two crystals' species factors (the wedge's mean factor
    for a species without its own; 1 where the wedge has no repeats)."""
    f = np.array([by_species.get(species.get(x), mean_factor) for x in names], dtype=float)
    return np.sqrt(np.outer(f, f))


def stouffer(z, valid, members):
    iu = np.triu_indices(len(members), 1)
    blk = z[np.ix_(members, members)][iu]
    ok = valid[np.ix_(members, members)][iu]
    return float(blk[ok].sum() / np.sqrt(ok.sum())) if ok.any() else -np.inf


def max_stouffer(z, valid, k, exhaustive=5000):
    """Largest Z over size-k subsets: exhaustive when there are at most `exhaustive` subsets,
    otherwise the best of greedy growth from every pair (sum of z within the group, divided by
    the square root of its number of valid pairs)."""
    n = z.shape[0]
    zz = np.where(valid, z, 0.0)
    vv = valid.astype(float)
    if math.comb(n, k) <= exhaustive:
        c = np.array(list(itertools.combinations(range(n), k)))
        blk = zz[c[:, :, None], c[:, None, :]]
        nb = vv[c[:, :, None], c[:, None, :]]
        tot = blk.sum((1, 2)) / 2
        num = nb.sum((1, 2)) / 2
        return float(np.max(np.where(num > 0, tot / np.sqrt(np.maximum(num, 1)), -np.inf)))
    order, within, _ = grow_all(zz, k)
    sums = within[:, k] * (k * (k - 1) / 2)
    best = -np.inf
    for s in np.argsort(sums)[::-1][:20]:  # recount valid pairs for the leading paths
        best = max(best, stouffer(z, valid, order[s, :k]))
    return best


def grow_iv(x, u2, kmax):
    """Greedy growth from every pair on x = z * u (u = 1 / se, zero for pairs left out), adding
    at each step the data set that raises the numerator most; returns, per path and size k,
    the inverse-variance statistic sum(z * u) / sqrt(sum(u^2)) over the pairs of the group
    ([paths, kmax + 1], -inf where undefined)."""
    n = x.shape[0]
    seeds = np.array(np.triu_indices(n, 1)).T
    s = len(seeds)
    ar = np.arange(s)
    inset = np.zeros((s, n), bool)
    inset[ar, seeds[:, 0]] = True
    inset[ar, seeds[:, 1]] = True
    tot = x[seeds[:, 0]] + x[seeds[:, 1]]
    totu = u2[seeds[:, 0]] + u2[seeds[:, 1]]
    num = x[seeds[:, 0], seeds[:, 1]].copy()
    den = u2[seeds[:, 0], seeds[:, 1]].copy()
    out = np.full((s, kmax + 1), -np.inf)
    for k in range(2, kmax + 1):
        out[:, k] = np.where(den > 0, num / np.sqrt(np.maximum(den, 1e-300)), -np.inf)
        if k == kmax:
            break
        cand = np.where(inset, -np.inf, tot)
        best = cand.argmax(1)
        num = num + tot[ar, best]
        den = den + totu[ar, best]
        inset[ar, best] = True
        tot = tot + x[best]
        totu = totu + u2[best]
    return out


def top_subsets(eS, k, N, exhaustive=200_000):
    """The N most coherent size-k subsets (mean residual within; all subsets when there are
    at most `exhaustive`, otherwise the distinct groups of greedy growth from every pair)."""
    n = eS.shape[0]
    ez = np.nan_to_num(eS, nan=0.0)
    if k < 2 or n < k:
        return []
    if math.comb(n, k) <= exhaustive:
        best = []
        combos = np.array(list(itertools.combinations(range(n), k)))
        for start in range(0, len(combos), 50_000):
            c = combos[start : start + 50_000]
            blk = ez[c[:, :, None], c[:, None, :]]
            vals = (blk.sum((1, 2)) - np.einsum("bii->b", blk)) / (k * (k - 1))
            top = np.argsort(vals)[::-1][:N]
            best += [(float(vals[j]), c[j]) for j in top]
        best.sort(key=lambda x: -x[0])
        return best[:N]
    order, within, _ = grow_all(ez, k)
    seen, out = set(), []
    for sidx in np.argsort(within[:, k])[::-1]:
        g = tuple(sorted(order[sidx, :k].tolist()))
        if g in seen:
            continue
        seen.add(g)
        out.append((float(within[sidx, k]), np.array(g)))
        if len(out) == N:
            break
    return out


def subgroup_candidates(T, M, smin, exhaustive=5000):
    """Sub-groups of the minority M that the pre-registered definition of "found" accepts
    (amendment 14: at least max(2, k/2) of its members): per size s from smin to |M| - 1, the
    subset with the largest standardised block sum and the one with the largest
    inverse-variance statistic (all subsets when there are at most `exhaustive`, otherwise
    backward elimination from M, removing the member whose removal raises the statistic most)."""
    M = np.asarray(M)
    stats = (lambda G: stouffer(T.z, T.valid, G), T.iv)
    out = {}
    sizes = range(max(2, smin), len(M))
    if sum(math.comb(len(M), s) for s in sizes) <= exhaustive:
        for s in sizes:
            combos = [M[list(c)] for c in itertools.combinations(range(len(M)), s)]
            for f in stats:
                vals = [f(G) for G in combos]
                G = combos[int(np.argmax(vals))]
                out[tuple(sorted(G.tolist()))] = G
    else:
        for f in stats:
            cur = list(M)
            while len(cur) > max(2, smin):
                vals = [f(np.array(cur[:i] + cur[i + 1:])) for i in range(len(cur))]
                cur.pop(int(np.argmax(vals)))
                G = np.array(cur)
                out[tuple(sorted(G.tolist()))] = G
    return list(out.values())


def _t_any(T, G, iv=False):
    """G's statistic standardised by the noise maxima of its size (as in p_noise_any_size)."""
    k = len(G)
    sims = T.null_all_sizes_iv() if iv else T.null_all_sizes()
    if k < 2 or k - 2 >= sims.shape[1]:
        return -np.inf
    sims = np.where(np.isfinite(sims), sims, np.nan)
    mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
    x = T.iv(G) if iv else stouffer(T.z, T.valid, G)
    return (x - mu[k - 2]) / sd[k - 2]


def blind_best(T):
    """The group a search that does not know the members finds (amendment 15): for each statistic,
    greedy growth from every pair, every size from 2 to n // 2, each size standardised by the noise
    maxima of that size; the group with the largest standardised statistic."""
    n = T.z.shape[0]
    kmax = n // 2
    out = {}
    for iv in (False, True):
        sims = T.null_all_sizes_iv() if iv else T.null_all_sizes()
        sims = np.where(np.isfinite(sims), sims, np.nan)
        mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
        best = (-np.inf, None)
        if iv:
            x = np.where(T.valid, T.z, 0.0) * T.u
            g = grow_iv(x, T.u ** 2, kmax)
            seeds = np.array(np.triu_indices(n, 1)).T
            t = (g[:, 2:] - mu[None]) / sd[None]
            t = np.where(np.isfinite(t), t, -np.inf)
            path, j = np.unravel_index(int(np.argmax(t)), t.shape)
            k = j + 2
            members = [int(seeds[path][0]), int(seeds[path][1])]
            tot = x[members[0]] + x[members[1]]
            while len(members) < k:
                cand = np.where(np.isin(np.arange(n), members), -np.inf, tot)
                b = int(np.argmax(cand))
                members.append(b)
                tot = tot + x[b]
            best = (float(t[path, j]), np.array(sorted(members)))
        else:
            order, within, _ = grow_all(np.where(T.valid, T.z, 0.0), kmax)
            for k in range(2, kmax + 1):
                sums = within[:, k] * (k * (k - 1) / 2)
                for b in np.argsort(sums)[::-1][:5]:
                    G = np.sort(order[b, :k])
                    tk = (stouffer(T.z, T.valid, G) - mu[k - 2]) / sd[k - 2]
                    if tk > best[0]:
                        best = (float(tk), G)
        out["iv" if iv else "stouffer"] = best
    return out


def mixed_candidates(T, bases, A):
    """Groups that "found" accepts beside the minority's own sub-groups (amendment 14): a
    qualifying group with crystals of the majority added, fewer than its minority members. From
    each base group, majority crystals are added one at a time, each time the one that raises
    the group's standardised statistic most (for each statistic); the extension that stands out
    most is kept if it stands out more than its base."""
    out = {}
    for iv in (False, True):
        for G in bases:
            G = [int(x) for x in G]
            nm = len(G)
            cur, rest = list(G), [int(a) for a in A if int(a) not in set(G)]
            best_t, best = _t_any(T, np.array(G), iv), None
            while rest and (len(cur) + 1 - nm) < nm:
                vals = [_t_any(T, np.array(cur + [a]), iv) for a in rest]
                j = int(np.argmax(vals))
                cur.append(rest.pop(j))
                if vals[j] > best_t:
                    best_t, best = vals[j], np.array(cur)
            if best is not None:
                out[tuple(sorted(best.tolist()))] = (best, nm)
    return list(out.values())


class Trial:
    def __init__(self, r, w, deltas, factor):
        D = Data(r, w, deltas)
        self.e, self.rho, se, noise, self.valid = model(D)
        self.se = factor * se
        self.z = np.where(self.valid, self.e / np.where(self.valid, self.se, 1.0), 0.0)
        # standardised noise replicates (the factor cancels: noise and se scale together)
        self.zb = np.where(self.valid[None], noise / np.where(self.valid, se, 1.0)[None], 0.0)
        self.eS = np.where(self.valid, self.e, np.nan)
        self.nb = factor * noise  # calibrated noise replicates of the residuals
        self.factor = factor
        self._null = {}
        # inverse-variance weights of the pairs (amendment 14): u = 1 / se, with the calibrated
        # error (a factor common to all pairs cancels from the statistic)
        self.u = np.where(self.valid, 1.0 / np.where(self.valid, self.se, 1.0), 0.0)

    def se_mean(self, members):
        """Counting-noise standard error of the mean residual over the pairs of `members`."""
        members = np.asarray(members)
        iu = np.triu_indices(len(members), 1)
        ok = self.valid[np.ix_(members, members)][iu]
        if not ok.any():
            return np.nan
        blk = self.nb[:, members][:, :, members][:, iu[0], iu[1]][:, ok]
        return float(blk.mean(1).std())

    def null(self, k):
        if k not in self._null:
            self._null[k] = np.array([max_stouffer(x, self.valid, k) for x in self.zb])
        return self._null[k]

    def p_noise(self, members):
        zm = stouffer(self.z, self.valid, members)
        return zm, gumbel_p(zm, self.null(len(members)))

    def null_all_sizes(self):
        """Per replicate, the largest standardised block sum of every size 2..n//2 (greedy
        growth from every pair, valid pairs only), for the look-elsewhere over sizes."""
        if not hasattr(self, "_null_sizes"):
            n = self.z.shape[0]
            kmax = n // 2
            v = self.valid.astype(float)
            out = []
            for x in self.zb:
                xx = np.where(self.valid, x, 0.0)
                order, within, _ = grow_all(xx, kmax)
                row = []
                for k in range(2, kmax + 1):
                    sums = within[:, k] * (k * (k - 1) / 2)
                    best = np.argsort(sums)[::-1][:5]
                    row.append(max(stouffer(x, self.valid, order[b, :k]) for b in best))
                out.append(row)
            self._null_sizes = np.array(out)  # [B, kmax - 1]
        return self._null_sizes

    def p_noise_any_size(self, members, model_noise=False):
        """p-value of the minority's block sum when the size is not known either: its
        standardised statistic (per-size mean and spread of the noise maxima) against the
        maximum over sizes in the noise replicates (Gumbel tail), as the procedures' scan.
        With model_noise, the noise is the analytic model at face value (what one experiment
        can know without repeat measurements), not rescaled by the repeat calibration."""
        k = len(members)
        sims = self.null_all_sizes()
        if k - 2 >= sims.shape[1] or k < 2:
            return np.nan
        mu, sd = sims.mean(0), sims.std(0) + 1e-12
        z = self.z * self.factor if model_noise else self.z
        t = (stouffer(z, self.valid, members) - mu[k - 2]) / sd[k - 2]
        return gumbel_p(t, ((sims - mu) / sd).max(1))

    def iv(self, members, z=None):
        """Inverse-variance statistic of a group: sum(z_ij u_ij) / sqrt(sum(u_ij^2)) over its
        pairs, the most powerful against a common excess correlation (amendment 14)."""
        z = self.z if z is None else z
        members = np.asarray(members)
        iu = np.triu_indices(len(members), 1)
        zz = z[np.ix_(members, members)][iu]
        uu = self.u[np.ix_(members, members)][iu]
        den = (uu ** 2).sum()
        return float((zz * uu).sum() / np.sqrt(den)) if den > 0 else -np.inf

    def null_all_sizes_iv(self):
        """As null_all_sizes, for the inverse-variance statistic."""
        if not hasattr(self, "_null_sizes_iv"):
            n = self.z.shape[0]
            kmax = n // 2
            u2 = self.u ** 2
            out = []
            for x in self.zb:
                g = grow_iv(np.where(self.valid, x, 0.0) * self.u, u2, kmax)
                out.append(g[:, 2:].max(0))
            self._null_sizes_iv = np.array(out)  # [B, kmax - 1]
        return self._null_sizes_iv

    def p_noise_any_size_iv(self, members, model_noise=False):
        k = len(members)
        sims = self.null_all_sizes_iv()
        if k - 2 >= sims.shape[1] or k < 2:
            return np.nan
        sims = np.where(np.isfinite(sims), sims, np.nan)
        mu, sd = np.nanmean(sims, 0), np.nanstd(sims, 0) + 1e-12
        z = self.z * self.factor if model_noise else self.z
        t = (self.iv(members, z) - mu[k - 2]) / sd[k - 2]
        return gumbel_p(t, np.nanmax((sims - mu) / sd, 1))

    def ratio(self, members, population):
        """S(G) over the most coherent subset of `population` (centred decoy), with the
        counting-noise standard errors of S(G), of the decoy and of the ratio."""
        population = np.asarray(population)
        sG = _mean_within(self.eS, members)
        d = decoy(self.eS, [population], len(members))
        kd = min(len(members), len(population) // 2)
        _, dm = best_subset(self.eS, population, kd)
        se_s = self.se_mean(members)
        se_d = self.se_mean(dm) if dm is not None else np.nan
        r = sG / d if d > 0 else np.inf
        se_r = (abs(r) * np.sqrt((se_s / sG) ** 2 + (se_d / d) ** 2)
                if np.isfinite(r) and sG != 0 and d > 0 else np.nan)
        return sG, r, d, se_s, se_d, se_r


def comparator_sizes(n):
    """Sizes of the single-population groups compared with minorities: every size from 2 to n // 2,
    so that every recorded sub-group of a minority (amendment 14) has comparators of its size
    (amendment 24; before, only the design's minority sizes)."""
    return range(2, n // 2 + 1)


def comparators(T, names, sizes):
    """The N_COMPARATORS most coherent subsets of each size of a single-population trial, with their
    separation from the rest and their standing against the noise."""
    comps = []
    n = len(names)
    for k in sizes:
        for S, G in top_subsets(T.eS, k, N_COMPARATORS):
            rest = np.setdiff1d(np.arange(n), G)
            _, ratio, d, se_s, se_d, se_r = T.ratio(G, rest)
            Z, p = T.p_noise(G)
            comps.append(dict(k=k, S=S, ratio=ratio, decoy=d, Z=Z, p_noise=p, se_S=se_s,
                              se_decoy=se_d, se_ratio=se_r, members=[names[i] for i in G]))
    return comps


class _Context:
    """What every trial of a design needs: design, matrices, calibration."""

    def __init__(self, design):
        self.trials = load_design(design)
        self.mats = load(matrix_design(design))
        self.cal = calibration(design)
        self.cal_sp = calibration_by_species(design)
        self.species = species_of(base_design(design))

    def trial(self, tid, boot):
        t = self.trials[tid]
        m = self.mats[tid]
        names = [str(x) for x in m["names"]]
        factor = self.cal.get(float(t["w"]), 1.0)
        if CAL_BY_SPECIES and float(t["w"]) in self.cal_sp:  # the noise of each pair's species (amendment 16)
            factor = pair_factors(names, self.species, self.cal_sp[float(t["w"])], factor)
        T = Trial(m["rij"].astype(float), m["wij"].astype(float), boot.astype(float), factor)
        return t, m, names, factor, T


def supplement(design, shard=0, nshards=1):
    """Amendment 24: add to the single-population records already written the comparators of the
    sizes comparator_sizes gives and they lack (the same computation as run's)."""
    out = RESULTS / "limits" / f"{design}_refined" / f"shard_{shard:03d}.jsonl"
    if not out.exists():
        return
    recs = [json.loads(l) for l in open_result(out)]
    need = [r for r in recs if r["kind"] == "single"
            and set(comparator_sizes(r["n"])) - {c["k"] for c in r.get("comparators", [])}]
    if not need:
        return
    ctx = _Context(design)
    boots = bload(matrix_design(design), "analytic", ids={r["id"] for r in need})
    for r in need:
        _, _, names, _, T = ctx.trial(r["id"], boots[r["id"]])
        missing = sorted(set(comparator_sizes(r["n"])) - {c["k"] for c in r.get("comparators", [])})
        r["comparators"] = sorted(r.get("comparators", []) + comparators(T, names, missing), key=lambda c: c["k"])
    tmp = out.with_suffix(".tmp")
    with open_result(tmp, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r, default=float) + "\n")
    tmp.replace(out)


def run(design, shard=0, nshards=1):
    ctx = _Context(design)
    trials, mats = ctx.trials, ctx.mats
    out = RESULTS / "limits" / f"{design}_refined" / f"shard_{shard:03d}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {json.loads(l)["id"] for l in open_result(out)} if out.exists() else set()
    ids = sorted(t for t in trials if t in mats and
                 (trials[t]["kind"] == "single" or trials[t]["k"] >= 2))[shard::nshards]
    boots = bload(matrix_design(design), "analytic", ids=set(ids))
    for tid in ids:
        if tid in done or tid not in boots:
            continue
        t, m, names, factor, T = ctx.trial(tid, boots[tid])
        n = len(names)
        rec = dict(id=tid, kind=t["kind"], w=t["w"], noise_factor=float(np.mean(factor)), n=n)
        if t["kind"] == "single":
            rec["comparators"] = comparators(T, names, comparator_sizes(n))
        else:
            pos = {x: i for i, x in enumerate(names)}
            M = np.array([pos[x] for x in t["minority"] if x in pos])
            A = np.array([pos[x] for x in t["majority"] if x in pos])
            rec["k"] = int(t["k"])
            rec["k_used"] = int(len(M))
            if len(M) >= 2 and len(A) >= 4:
                S, ratio, d, se_s, se_d, se_r = T.ratio(M, A)
                # the trial's own majority: its most coherent group of the minority's size
                Sd, Gd = best_subset(T.eS, A, len(M))
                if Gd is not None:
                    _, pd = T.p_noise(Gd)
                    rec.update(majority_group_S=Sd, majority_group_p=pd,
                               majority_group_se_S=T.se_mean(Gd),
                               majority_group=[names[i] for i in Gd])
                Z, p = T.p_noise(M)
                p_any = T.p_noise_any_size(M)
                p_model = T.p_noise_any_size(M, model_noise=True)
                rec.update(Z_iv=T.iv(M), p_noise_any_size_iv=T.p_noise_any_size_iv(M),
                           p_noise_model_iv=T.p_noise_any_size_iv(M, model_noise=True))
                # groups that would count as finding the minority (amendment 14): its own
                # qualifying sub-groups, and those groups with fewer majority crystals added
                pure = subgroup_candidates(T, M, int(math.ceil(max(2, len(M) / 2))))
                mixed = mixed_candidates(T, [M] + pure, A)
                subs = []
                for G, n_min in [(G, len(G)) for G in pure] + mixed:
                    pop = np.array([a for a in A if a not in set(G.tolist())])
                    gS, gr, gd, gse_s, gse_d, gse_r = T.ratio(G, pop)
                    sub = dict(members=[names[i] for i in G], k=int(len(G)), n_minority=int(n_min),
                               S=gS, ratio=gr,
                               decoy=gd, se_S=gse_s, se_decoy=gse_d, se_ratio=gse_r,
                               Z=stouffer(T.z, T.valid, G), Z_iv=T.iv(G),
                               p_noise_any_size=T.p_noise_any_size(G),
                               p_noise_model=T.p_noise_any_size(G, model_noise=True),
                               p_noise_any_size_iv=T.p_noise_any_size_iv(G),
                               p_noise_model_iv=T.p_noise_any_size_iv(G, model_noise=True))
                    Sg, Gg = best_subset(T.eS, pop, len(G))
                    if Gg is not None:
                        sub.update(majority_group_S=Sg, majority_group_p=T.p_noise(Gg)[1],
                                   majority_group_se_S=T.se_mean(Gg))
                    subs.append(sub)
                rec["subgroups"] = subs
                # the group found blind (amendment 15), against the minority's qualifying groups;
                # with the calibrated noise and with the noise one experiment knows (the model's)
                minority, k_need = set(M.tolist()), max(2, len(M) / 2)
                units = [M] + [G for G, _ in [(G, len(G)) for G in pure] + mixed]
                T1 = T if np.all(np.asarray(T.factor) == 1.0) else Trial(
                    m["rij"].astype(float), m["wij"].astype(float), boots[tid].astype(float), 1.0)
                for key, TT in (("blind", T), ("blind_model", T1)):
                    blind = {}
                    for stat, (tb, Gb) in blind_best(TT).items():
                        if Gb is None:
                            continue
                        hit = len(minority & set(Gb.tolist()))
                        t_units = max(_t_any(TT, G, stat == "iv") for G in units)
                        blind[stat] = dict(t=tb, members=[names[i] for i in Gb], k=int(len(Gb)), hit=int(hit),
                                           qualifies=bool(hit >= k_need and hit > len(Gb) / 2),
                                           t_minority_best=float(t_units))
                    rec[key] = blind
                rec.update(S=S, ratio=ratio, decoy=d, Z=Z, p_noise=p, p_noise_any_size=p_any,
                           p_noise_model=p_model,
                           within_model_noise=bool(not p_model < ALPHA),
                           se_S=se_s, se_decoy=se_d, se_ratio=se_r,
                           within_counting_noise=bool(p >= ALPHA),
                           within_counting_noise_any_size=bool(not p_any < ALPHA))
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")


def _not_above(a, b, se_a, se_b, tol):
    """a is not larger than b by more than `tol` standard errors of the difference."""
    if np.isnan(a) or np.isnan(b):
        return False
    if a == -np.inf or b == np.inf:
        return True
    if a == np.inf or b == -np.inf:
        return False
    se = np.sqrt(np.nan_to_num(se_a, nan=0.0) ** 2 + np.nan_to_num(se_b, nan=0.0) ** 2)
    return bool(a - b <= tol * se)


def _pmin(*ps):
    ps = [p for p in ps if p is not None and np.isfinite(p)]
    return min(ps) if ps else None


def _unit(u, comps, wedges, tol, alpha, statistic, tid):
    """Verdicts for one group (the minority or one of its sub-groups): within the counting
    noise (calibrated; and with the model noise), and the single-population group or the trial
    majority's group that dominates it, if any."""
    dom = None
    # structure inside one population seen with at least as much data per crystal is there
    # with less data too (amendment 8)
    pool = [c for w in wedges if w >= float(u["w"]) for c in comps.get((w, u["k"]), [])]
    for c in pool:
        if (_not_above(u["S"], c["S"], u.get("se_S"), c.get("se_S"), tol)
                and _not_above(u["ratio"], c["ratio"], u.get("se_ratio"), c.get("se_ratio"), tol)):
            dom = c["trial"]
            break
    p_s = u.get("p_noise_any_size", u.get("p_noise"))
    p_iv = u.get("p_noise_any_size_iv")
    p_any = {"stouffer": p_s, "iv": p_iv, "either": _pmin(p_s, p_iv)}[statistic]
    test = bool(not (p_any is not None and p_any < alpha))
    pm_s = u.get("p_noise_model", p_s)
    pm_iv = u.get("p_noise_model_iv")
    p_model = {"stouffer": pm_s, "iv": pm_iv, "either": _pmin(pm_s, pm_iv)}[statistic]
    model_noise = bool(not (p_model is not None and p_model < alpha))
    # structure inside the trial's own majority (amendment 11): a group of the same size,
    # standing out of the noise, at least as coherent within noise
    if dom is None and u.get("majority_group_p", 1.0) < alpha:
        if _not_above(u["S"], u["majority_group_S"], u.get("se_S"), u.get("majority_group_se_S"), tol):
            dom = tid + " (its own majority)"
    return dict(noise=test, test=test, model_noise=model_noise, dom=dom)


_CLASSIFIED, _ROWS = {}, {}


def _rows(design):
    """The refined records of a design (read once per process)."""
    if design not in _ROWS:
        import glob

        _ROWS[design] = [json.loads(l) for f in glob_results(str(RESULTS / "limits" / f"{design}_refined" / "shard_*.jsonl"))
                         for l in open_result(f)]
    return _ROWS[design]


def classify(design, tol=TOL, found_ids=frozenset(), alpha=None, statistic="either", subgroups=True,
             membership=True, cross=False, design_sizes_only=False):
    """classify_uncached, computed once per process for the same arguments."""
    key = (design, tol, frozenset(found_ids), alpha, statistic, subgroups, membership, cross, design_sizes_only)
    if key not in _CLASSIFIED:
        _CLASSIFIED[key] = classify_uncached(design, tol, found_ids, alpha, statistic, subgroups, membership,
                                             cross, design_sizes_only)
    return _CLASSIFIED[key]


def classify_uncached(design, tol=TOL, found_ids=frozenset(), alpha=None, statistic="either", subgroups=True,
                      membership=True, cross=False, design_sizes_only=False):
    """Join minorities with the comparators of their collection and wedge (amendments 4, 6, 7,
    10, 11, 13, 14). `alpha` is the level of every counting-noise test (default LIMIT_ALPHA).
    `statistic`: the counting-noise test's statistic, "stouffer" (the standardised block sum,
    amendments 4-13), "iv" (inverse-variance weighted) or "either" (amendment 14: a minority
    is within the noise only if neither stands out). With `subgroups` (amendment 14) the
    minority is a data limit only if it and every recorded sub-group that would count as
    finding it (at least max(2, k/2) members) are: counting noise if all are within the noise,
    structure inside one species otherwise.

    A minority is within the counting noise if its standardised block sum does not stand out
    at ALPHA against the noise replicates' maxima over all subsets and all sizes (an
    unsupervised procedure knows neither), unless a frozen procedure with no more false groups
    than today's DIALS found it
    in that trial (`found_ids`: then it was demonstrably distinguishable); within the
    structure of one population if a qualifying single-population group of the same size is
    at least as coherent and as separated, both within `tol` standard errors (tol = 0: the
    strict comparison of amendment 4)."""
    import glob

    alpha = LIMIT_ALPHA if alpha is None else alpha
    rows = _rows(design)
    # design_sizes_only: comparators of the design's minority sizes only, as before amendment 24
    old_sizes = {t["k"] for t in load_design(base_design(design)).values() if t["kind"] == "mixture"}
    keep_size = (lambda k: k in old_sizes) if design_sizes_only else (lambda k: True)
    comps = {}
    for r in rows:
        if r["kind"] == "single":
            for c in r.get("comparators", []):
                if c["p_noise"] < alpha and r["n"] - c["k"] >= MIN_REST and keep_size(c["k"]):
                    comps.setdefault((float(r["w"]), c["k"]), []).append(dict(c, trial=r["id"]))
    if cross:  # amendment 18: the other collection's single-population groups, at every wedge
        other = {"cryo": "rt", "rt": "cryo"}[base_design(design)]
        suffix = design[len(base_design(design)):]
        for r in _rows(f"{other}{suffix}"):
            if r["kind"] != "single":
                continue
            for c in r.get("comparators", []):
                if c["p_noise"] < alpha and r["n"] - c["k"] >= MIN_REST and keep_size(c["k"]):
                    comps.setdefault((np.inf, c["k"]), []).append(dict(c, trial=r["id"]))
    wedges = sorted({w for w, _ in comps})
    out = {}
    for r in rows:
        if r["kind"] != "mixture":
            continue
        if "S" not in r:
            if r.get("k_used", 0) < 2 <= r.get("k", 0):  # cosym kept fewer than two of its crystals
                out[r["id"]] = dict(r, within_counting_noise=True, within_counting_noise_test=True,
                                    within_single_experiment_noise=False, within_population_structure=False,
                                    dominated_by=None, membership_unresolved=False, data_limit=True,
                                    detectable_unit=None, too_few_crystals_in_data=True)
            continue
        units = [dict(r, k=r["k_used"])] + ([dict(s, w=r["w"]) for s in r.get("subgroups", [])]
                                             if subgroups else [])
        verdicts = [_unit(u, comps, wedges, tol, alpha, statistic, r["id"]) for u in units]
        if r["id"] in found_ids:  # found by a procedure without more false groups: separable
            verdicts = [dict(v, noise=False, model_noise=False, dom=None) for v in verdicts]
        ok = [not v["noise"] and v["dom"] is None for v in verdicts]
        ok_model = [not v["model_noise"] and v["dom"] is None for v in verdicts]
        own = verdicts[0]
        noise = all(v["noise"] for v in verdicts)
        dom = None if any(ok) or noise else next(v["dom"] for v in verdicts if v["dom"] is not None)
        best = int(np.argmax(ok)) if any(ok) else None
        # membership not resolved (amendment 15): searched blind, with either statistic, the most
        # significant group mixes the minority with majority crystals (it would not count as
        # finding it) and stands out more than any group that would
        def blocked(bl):
            return bool(bl and all(not b["qualifies"] and b["hit"] >= 1 and b["t"] > b["t_minority_best"]
                                   for b in bl.values()))
        unresolved = bool(membership and any(ok) and r["id"] not in found_ids and blocked(r.get("blind")))
        # with the noise one experiment knows, every check with that noise (amendment 16)
        unresolved_model = bool(membership and r["id"] not in found_ids and blocked(r.get("blind_model")))
        separable_model = any(ok_model) and not unresolved_model
        out[r["id"]] = dict(r, within_counting_noise=noise, within_counting_noise_test=own["test"],
                            within_single_experiment_noise=bool(any(ok) and not unresolved and not separable_model),
                            within_population_structure=dom is not None, dominated_by=dom,
                            membership_unresolved=unresolved,
                            data_limit=bool(not any(ok) or unresolved),
                            detectable_unit=None if best is None else (
                                "minority" if best == 0 else units[best]["members"]))
    return out


if __name__ == "__main__":
    task, ntasks = int(os.environ.get("CF_SHARD", 0)), int(os.environ.get("CF_NSHARDS", 1))
    if sys.argv[1] == "supplement":  # python -m cf.limits3 supplement DESIGN
        supplement(sys.argv[2], task, ntasks)
    else:
        run(sys.argv[1], task, ntasks)
