"""Post-hoc refinement of the limits (written after the oracle results; see REPORT.md):
what an UNSUPERVISED procedure with a low false-group rate can find at best.

For a mixture trial with minority M (size k), on the rank-1 residuals e of DIALS's rij:
1. Counting noise with look-elsewhere: each residual is standardised by its bootstrap
   standard error (z = e / SE) and a group scored by the standardised sum of its pairs,
   Z(G) = sum z_ij / sqrt(#pairs) (the most powerful statistic for a group against
   independent counting noise). The null distribution is the largest Z over all size-k
   subsets of the trial's crystals in pure counting noise (bootstrap replicates over
   reflections, centred: (e_b - e) / SE); p_noise = upper tail of Z(M) (Gumbel fit to the
   replicate maxima). A minority with p_noise >= alpha cannot be reported by any procedure
   that keeps false groups from counting noise below alpha.
2. Real structure inside one population: the largest S over size-k subsets of every
   single-population trial of the same collection and wedge (all species, parts and
   segments): family_decoy. A minority with S(M) <= family_decoy is no more coherent than a
   group that exists inside one species of the same experiment, so a procedure that reported
   it would report that group too.
Writes results/limits/<design>_unsupervised.jsonl.
"""

from __future__ import annotations

from cf.result_io import open_result

import json
import sys

import numpy as np

from cf.analyze import load_design
from cf.bootstrap import load as bload
from cf.cmatrices import load
from cf.oracle import decoy_max_resid, group_means_nan, rank1_fit
from cf.paths import RESULTS


def residuals(r, w):
    w = np.array(w, float)
    np.fill_diagonal(w, 0)
    rho = rank1_fit(r, w)
    e = r - np.outer(rho, rho)
    e[w <= 0] = np.nan
    np.fill_diagonal(e, np.nan)
    return e


def max_subset(e, k, limit=50_000):
    return decoy_max_resid(e, np.arange(e.shape[0]), k, exhaustive_limit=limit)["max"]


def stouffer(z, members):
    blk = z[np.ix_(members, members)][np.triu_indices(len(members), 1)]
    blk = blk[np.isfinite(blk)]
    return float(blk.sum() / np.sqrt(len(blk))) if len(blk) else -np.inf


def max_stouffer(z, k, limit=5_000):
    """Largest standardised block sum over size-k subsets (exhaustive, or greedy+swaps)."""
    import itertools
    import math

    n = z.shape[0]
    zz = np.nan_to_num(z, nan=0.0)
    cnt = np.isfinite(z).astype(float)
    if math.comb(n, k) <= limit:
        best = -np.inf
        combos = np.array(list(itertools.combinations(range(n), k)))
        for start in range(0, len(combos), 20_000):
            c = combos[start : start + 20_000]
            blk = zz[c[:, :, None], c[:, None, :]]
            nb = cnt[c[:, :, None], c[:, None, :]]
            tot = (blk.sum((1, 2)) - np.einsum("bii->b", blk)) / 2
            num = (nb.sum((1, 2)) - np.einsum("bii->b", nb)) / 2
            vals = np.where(num > 0, tot / np.sqrt(np.maximum(num, 1)), -np.inf)
            best = max(best, float(vals.max()))
        return best
    best = -np.inf
    for a in range(n):
        cur = [a]
        while len(cur) < k:
            rest = [x for x in range(n) if x not in cur]
            cur.append(rest[int(np.argmax(zz[np.ix_(rest, cur)].sum(1)))])
        improved = True
        while improved:
            improved = False
            for pos in range(k):
                others = cur[:pos] + cur[pos + 1 :]
                cand = [x for x in range(n) if x not in cur]
                gains = zz[np.ix_(cand, others)].sum(1)
                if gains.max() > zz[cur[pos], others].sum() + 1e-12:
                    cur[pos] = cand[int(np.argmax(gains))]
                    improved = True
        best = max(best, stouffer(z, np.array(cur)))
    return best


def gumbel_p(x, sample):
    sample = np.asarray(sample)
    beta = np.sqrt(6.0) * sample.std() / np.pi + 1e-15
    loc = sample.mean() - 0.5772156649 * beta
    return float(-np.expm1(-np.exp(-(x - loc) / beta)))


def run(design, alpha=0.01, B=100):
    trials = load_design(design)
    mats = load(design)
    boots = bload(design)
    # family decoys from single-population trials, per (wedge, k)
    singles = [t for t in trials.values() if t["kind"] == "single" and t["id"] in mats]
    ks = sorted({t["k"] for t in trials.values() if t["kind"] == "mixture" and t["k"] >= 2})
    family = {}
    for t in singles:
        m = mats[t["id"]]
        e = residuals(m["rij"].astype(float), m["wij"].astype(float))
        for k in ks:
            if k > len(t["crystals"]) // 2:
                continue
            v = max_subset(e, k)
            key = (t["w"], k)
            if v > family.get(key, (-np.inf, None))[0]:
                family[key] = (v, t["id"])
    out = RESULTS / "limits" / f"{design}_unsupervised.jsonl"
    with open_result(out, "w") as fh:
        for tid, t in trials.items():
            if t["kind"] != "mixture" or t["k"] < 2 or tid not in mats or tid not in boots:
                continue
            m = mats[tid]
            names = [str(x) for x in m["names"]]
            minority = [names.index(x) for x in t["minority"] if x in names]
            if len(minority) < 2:
                continue
            r, w = m["rij"].astype(float), m["wij"].astype(float)
            e = residuals(r, w)
            k = len(minority)
            s_m = float(group_means_nan(e, np.array([minority]))[0])
            eb_all = []
            for d in boots[tid][:B].astype(float):
                rb = r + d
                np.fill_diagonal(rb, 0)
                eb_all.append(residuals(rb, w) - e)
            eb_all = np.array(eb_all)
            se = np.nanstd(eb_all, axis=0)
            ok = np.isfinite(e) & (se > 0)
            z = np.where(ok, e / np.where(ok, se, 1.0), np.nan)
            z_m = stouffer(z, minority)
            reps = [max_stouffer(np.where(ok, x / np.where(ok, se, 1.0), np.nan), k) for x in eb_all]
            p_noise = gumbel_p(z_m, reps)
            fam = family.get((t["w"], k), (np.nan, None))
            rec = dict(id=tid, k=k, S=s_m, Z=z_m, p_noise=p_noise, noise_q99=float(np.quantile(reps, 0.99)),
                       family_decoy=fam[0], family_decoy_trial=fam[1],
                       within_counting_noise=bool(p_noise >= alpha),
                       within_population_structure=bool(np.isfinite(fam[0]) and s_m <= fam[0]))
            rec["unsupervised_limit"] = rec["within_counting_noise"] or rec["within_population_structure"]
            fh.write(json.dumps(rec) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    run(sys.argv[1])
