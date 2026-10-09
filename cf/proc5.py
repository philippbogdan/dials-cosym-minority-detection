"""Procedure candidate P5 (the one proposed for DIALS, prototype on matrices).

Inputs: DIALS's pairwise correlation matrix rij, its weights wij, and B bootstrap replicates
of rij over reflections (counting noise only).

For a set of crystals:
1. Single population: rank-1 fit rho (weights wij); residuals e = r - rho rho^T; counting
   noise of e from the replicates (each refitted): standard errors and noise-only replicates.
2. Scan: from every pair, greedy growth; for each size the best gap
   (mean standardised residual within the group - best outsider's mean to the group);
   standardised per size by the noise replicates; T = max over sizes; p from a Gumbel fit to
   the replicates' maxima. The group G of the maximum is the candidate.
3. Population criterion: the candidate must stand out from the rest by more than the rest's
   own spread: U = delta / sigma_het, with delta the gap in disattenuated units
   (e_ij / rho_i rho_j, i.e. 1 - cos of the angle between crystals) and sigma_het the robust
   spread of the rest's disattenuated residuals in excess of counting noise. If the rest
   itself splits into valid groups, sigma_het is pooled within those groups.
4. A split is kept when p < alpha, U >= c and the rest has at least `min_pop` crystals; both
   sides are then analysed again. Otherwise the set is one population.
"""

from __future__ import annotations

import numpy as np

from cf.oracle import rank1_fit
from cf.scan import grow_all


class Data:
    def __init__(self, r, w, deltas):
        self.r = np.array(r, float)
        self.w = np.array(w, float)
        np.fill_diagonal(self.w, 0.0)
        self.d = np.array(deltas, float)

    def sub(self, idx):
        idx = np.asarray(idx)
        return Data(self.r[np.ix_(idx, idx)], self.w[np.ix_(idx, idx)],
                    self.d[:, idx][:, :, idx])


def _fit(r, w):
    rho = rank1_fit(r, w)
    e = r - np.outer(rho, rho)
    np.fill_diagonal(e, 0.0)
    return e, rho


def model(D: Data):
    n = D.r.shape[0]
    valid = D.w > 0
    np.fill_diagonal(valid, False)
    e, rho = _fit(D.r, D.w)
    eb = np.empty_like(D.d)
    for b in range(len(D.d)):
        rb = D.r + D.d[b]
        np.fill_diagonal(rb, 0.0)
        valid &= ~((rb == 0) & (D.w > 0))
        eb[b] = _fit(rb, D.w)[0]
    se = eb.std(0)
    valid &= se > 0
    # a pair whose correlation has no noise of its own in the replicates (too few
    # reflections for an error) is left out, not given the small spread of the refit
    valid &= np.any(D.d != 0, axis=0)
    return e, rho, se, eb - e[None], valid


def _gaps(z, valid, kmax, stat="gap"):
    order, within, outside = grow_all(np.where(valid, z, 0.0), kmax)
    if stat == "within":
        return within, order
    return within - outside, order


def _gumbel_null(sims):
    mu, sd = sims.mean(0), sims.std(0) + 1e-12
    t_sim = ((sims - mu) / sd).max(1)
    beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
    loc = t_sim.mean() - 0.5772156649 * beta
    return mu, sd, loc, beta


def floor_p(s, G, floor):
    """Look-elsewhere p-value that group G is more coherent than `floor` (correlation units):
    the mean standardised excess (e_ij - floor) / se_ij over its pairs, against the noise
    replicates' largest within-group means per size (scan(..., floor_null=True))."""
    k = len(G)
    if s.get("null_within") is None or k < 2 or k - 2 >= len(s["null_within"][0]):
        return 1.0
    mu, sd, loc, beta = s["null_within"]
    sse = np.where(s["valid"], s["se"], 1.0)
    z = np.where(s["valid"], (s["e"] - floor) / sse, 0.0)
    iu = np.triu_indices(k, 1)
    obs = z[np.ix_(G, G)][iu].sum() / (k * (k - 1) / 2)
    t = (obs - mu[k - 2]) / sd[k - 2]
    return float(-np.expm1(-np.exp(-(t - loc) / beta)))


def scan(D: Data, stat="gap", calibrate=False, floor_null=False):
    n = D.r.shape[0]
    kmax = n // 2
    if n < 4 or kmax < 2:
        return None
    e, rho, se, noise, valid = model(D)
    if valid.sum() < 6:
        return None
    if isinstance(calibrate, str) and calibrate.startswith("robust"):
        # noise level from the residuals' own robust scale (bounds the noise from above),
        # the replicates giving only its pattern over pairs; "robust:f" floors kappa at f
        iu = np.triu_indices(n, 1)
        ratio = np.abs(e[iu] / np.where(valid[iu], se[iu], 1.0))[valid[iu]]
        floor = float(calibrate.split(":")[1]) if ":" in calibrate else 0.0
        kappa = max(floor, float(np.median(ratio)) / 0.6745) if len(ratio) else 1.0
        se = kappa * se
        noise = kappa * noise
    elif calibrate is not True and calibrate:  # a given factor (noise calibrated elsewhere)
        se = float(calibrate) * se
        noise = float(calibrate) * noise
    elif calibrate:  # as dials.algorithms.correlation.groups.scan(calibrate=True)
        iu = np.triu_indices(n, 1)
        ratio = np.abs(e[iu] / np.where(valid[iu], se[iu], 1.0))[valid[iu]]
        kappa = max(1.0, float(np.median(ratio)) / 0.6745) if len(ratio) else 1.0
        se = kappa * se
        noise = kappa * noise
    sse = np.where(valid, se, 1.0)
    g, order = _gaps(np.where(valid, e / sse, 0.0), valid, kmax, stat)
    obs = np.nanmax(g[:, 2:], axis=0)
    sims = np.array([np.nanmax(_gaps(np.where(valid, nb / sse, 0.0), valid, kmax, stat)[0][:, 2:], axis=0)
                     for nb in noise])
    mu, sd = sims.mean(0), sims.std(0) + 1e-12
    t_obs = np.where(np.isfinite(obs), (obs - mu) / sd, -np.inf)
    t_sim = ((sims - mu) / sd).max(1)
    beta = np.sqrt(6.0) * t_sim.std() / np.pi + 1e-12
    loc = t_sim.mean() - 0.5772156649 * beta
    T = float(np.max(t_obs))
    p = float(-np.expm1(-np.exp(-(T - loc) / beta)))
    k = int(np.argmax(t_obs)) + 2
    members = np.sort(order[int(np.nanargmax(g[:, k])), :k])
    null_within = None
    if floor_null:  # null of the within-group mean of standardised residuals, per size
        sims_w = np.array([np.nanmax(_gaps(np.where(valid, nb / sse, 0.0), valid, kmax, "within")[0][:, 2:],
                                     axis=0) for nb in noise])
        null_within = _gumbel_null(sims_w)
    return dict(p=p, T=T, members=members, e=e, rho=rho, se=se, valid=valid, noise=noise,
                null_within=null_within)


def separation(s, G, rest):
    """Disattenuated gap of G from the rest, and the rest's excess spread."""
    rr = np.outer(s["rho"], s["rho"])
    with np.errstate(divide="ignore", invalid="ignore"):
        dz = np.where(s["valid"] & (rr > 0), s["e"] / rr, np.nan)
        dse = np.where(s["valid"] & (rr > 0), s["se"] / rr, np.nan)
    within = np.nanmean(dz[np.ix_(G, G)][np.triu_indices(len(G), 1)])
    outs = np.nanmean(dz[np.ix_(G, rest)], axis=0)
    delta = float(within - np.nanmax(outs))
    return delta, dz, dse


def spread(dz, dse, groups, z_ucb=1.645):
    """Excess spread of disattenuated residuals pooled within `groups` over counting noise,
    as an upper confidence bound: heterogeneity that the data cannot exclude is assumed to
    be present (relative error of the robust variance ~ sqrt(2 / m), m = crystals)."""
    vals, noise = [], []
    m_cryst = sum(len(g) for g in groups if len(g) >= 3)
    for grp in groups:
        if len(grp) < 3:
            continue
        iu = np.triu_indices(len(grp), 1)
        v = dz[np.ix_(grp, grp)][iu]
        nse = dse[np.ix_(grp, grp)][iu]
        ok = np.isfinite(v) & np.isfinite(nse)
        vals.append(v[ok] - np.median(v[ok]) if ok.any() else v[ok])
        noise.append(nse[ok])
    if not vals:
        return 0.0
    v = np.concatenate(vals)
    nse = np.concatenate(noise)
    if len(v) < 3:
        return 0.0
    mad = np.median(np.abs(v)) * 1.4826
    tot2 = mad**2 * (1.0 + z_ucb * np.sqrt(2.0 / max(m_cryst, 2)))
    het2 = tot2 - np.median(nse**2)
    return float(np.sqrt(het2)) if het2 > 0 else 0.0


def analyze(D: Data, alpha=0.01, c=5.0, min_pop=5, depth=0, trace=None):
    """Returns the list of groups (index arrays into D)."""
    n = D.r.shape[0]
    idx = np.arange(n)
    if n < 2 + min_pop or depth > 8:
        return [idx]
    s = scan(D)
    if s is None or s["p"] >= alpha:
        return [idx]
    G = s["members"]
    rest = np.setdiff1d(idx, G)
    if len(rest) < min_pop:
        return [idx]
    delta, dz, dse = separation(s, G, rest)
    sig = spread(dz, dse, [rest])
    ok = delta > 0 and (sig == 0.0 or delta / sig >= c)
    rest_groups = None
    if not ok and delta > 0:
        # the rest may hold populations of its own: judge against their pooled spread
        rest_groups = [rest[g] for g in analyze(D.sub(rest), alpha, c, min_pop, depth + 1)]
        if len(rest_groups) > 1:
            sig = spread(dz, dse, rest_groups)
            ok = sig == 0.0 or delta / sig >= c
    if trace is not None:
        trace.append(dict(depth=depth, n=n, p=s["p"], k=len(G), delta=delta, sigma=sig, ok=ok))
    if not ok:
        return [idx]
    if rest_groups is None:
        rest_groups = [rest[g] for g in analyze(D.sub(rest), alpha, c, min_pop, depth + 1)]
    g_groups = [G[g] for g in analyze(D.sub(G), alpha, c, min_pop, depth + 1)] if len(G) >= 2 + min_pop else [G]
    return g_groups + rest_groups
