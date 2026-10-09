"""Calibration of counting-noise models against independent repeat measurements.

The same crystals measured in disjoint wedge segments (single-population trials, segments
0-3) give independent repeat measurements of every pairwise correlation, each over its own
reflections (segments share 1-4% of their unique reflections at 0.2-0.4 degrees, 10-20% at 1
degree), so the scatter includes both the measurement errors and which reflections a pair happens
to share: the counting noise of one experiment. For every
crystal pair, the scatter of the single-population residual across segments is compared with
the standard error each noise model predicts for one segment:
  reflection  bootstrap over reflections (Poisson weights; cf.bootstrap)
  parametric  measurement errors redrawn from the scaled sigmas around the consensus
              (inverse-variance weighted mean of all crystals), reflections fixed
  analytic    the large-sample error of Fisher's z of each CC, 1 / sqrt(neff - 3), with
              neff DIALS's effective number of reflections of the pair (cf.bootstrap)
Writes results/noisecal/<design>.json.
"""

from __future__ import annotations

import json
import sys

import numpy as np

from cf.analyze import load_design
from cf.fastcc import dense_from_cosym, trial_arrays
from cf.oracle import rank1_fit
from cf.paths import DATA, RESULTS, base_design, representation


def wcc(y, v, m, wt=None):
    import importlib.util

    from cf.paths import ROOT

    path = ROOT / "cf" / "frozen" / "groups_v6_1.py"
    spec = importlib.util.spec_from_file_location("pg", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pg"] = mod
    spec.loader.exec_module(mod)
    return mod.weighted_cc_matrix(y, v, m, wt)


def residual(r, w):
    w = np.array(w, float)
    np.fill_diagonal(w, 0)
    rho = rank1_fit(r, w)
    e = r - np.outer(rho, rho)
    e[w <= 0] = np.nan
    np.fill_diagonal(e, np.nan)
    return e


def parametric_replicates(y, v, m, B=50, seed=0):
    rng = np.random.default_rng(seed)
    w = np.where(m, 1.0 / np.where(m, v, 1.0), 0.0)
    cons = (w * y).sum(0) / np.maximum(w.sum(0), 1e-300)
    sig = np.sqrt(v)
    out = []
    for _ in range(B):
        ys = np.where(m, cons[None, :] + sig * rng.standard_normal(y.shape), 0.0)
        out.append(wcc(ys, v, m))
    return np.array(out)


def reflection_replicates(y, v, m, B=50, seed=0):
    rng = np.random.default_rng(seed)
    return np.array([wcc(y, v, m, rng.poisson(1.0, y.shape[1]).astype(float)) for _ in range(B)])


def main(design):
    trials = load_design(design)
    base, rep = base_design(design), representation(design)
    out = {}
    wedges = sorted({t["w"] for t in trials.values() if t["kind"] == "single" and t.get("seg", 0) > 0})
    for w in wedges:
        for sp in ("bovine", "porcine", "human"):
            es, se_ref, se_par, se_kap, se_ana, names0 = [], [], [], [], [], None
            for seg in range(4):
                tid = f"{base}-single-{sp[0]}-all-w{w:g}" + (f"-s{seg}" if seg else "")
                f = DATA / "trials" / base / f"{tid}.npz"
                if tid not in trials or not f.exists():
                    continue
                d = np.load(f)
                a = trial_arrays(d, rep)
                names = a["names"]
                if names0 is None:
                    names0 = names
                if sorted(names) != sorted(names0):
                    continue
                idx = [names.index(x) for x in names0]
                y, v, m = a["y"], a["v"], a["m"]
                r = a["rij"]
                np.fill_diagonal(r, 0)
                wij = a["wij"]
                e = residual(r, wij)
                if rep == "intensity":
                    ref = np.array([residual(x, wij) for x in reflection_replicates(y, v, m)])
                    par = np.array([residual(x, wij) for x in parametric_replicates(y, v, m)])
                else:  # the bootstrap and parametric models are calibrated for intensities only
                    ref = par = np.full((2,) + r.shape, np.nan)
                from cf.bootstrap import analytic_replicates, neff_matrix

                ana = np.array([residual(r + x, wij) for x in
                                analytic_replicates(r, neff_matrix(y, v, m), B=50).astype(float)])
                sub = np.ix_(idx, idx)
                es.append(e[sub])
                se_ref.append(np.nanstd(ref, axis=0)[sub])
                sp_ = np.nanstd(par, axis=0)
                ok_ = np.isfinite(e) & (sp_ > 0)
                kappa = max(1.0, float(np.median(np.abs(e[ok_] / sp_[ok_])) / 0.6745))
                se_par.append(sp_[sub])
                se_kap.append((kappa * sp_)[sub])
                se_ana.append(np.nanstd(ana, axis=0)[sub])
            if len(es) < 3:
                continue
            es, se_ref, se_par, se_kap, se_ana = map(np.array, (es, se_ref, se_par, se_kap, se_ana))
            iu = np.triu_indices(es.shape[1], 1)
            seg_sd = np.nanstd(es[:, iu[0], iu[1]], axis=0, ddof=1)
            c4 = {3: 0.886, 4: 0.921}[len(es)]  # E[sample SD] = c4 * sigma
            rms_ref = np.sqrt(np.nanmean(se_ref[:, iu[0], iu[1]] ** 2, axis=0))
            rms_par = np.sqrt(np.nanmean(se_par[:, iu[0], iu[1]] ** 2, axis=0))
            rms_kap = np.sqrt(np.nanmean(se_kap[:, iu[0], iu[1]] ** 2, axis=0))
            se_ana = np.where(se_ana > 0, se_ana, np.nan)  # pairs without a usable error
            rms_ana = np.sqrt(np.nanmean(se_ana[:, iu[0], iu[1]] ** 2, axis=0))
            ok = np.isfinite(seg_sd) & (rms_ref > 0) & (rms_par > 0)
            out[f"{w:g}|{sp}"] = dict(
                n_segments=len(es), n_pairs=int(ok.sum()),
                observed_over_reflection=float(np.nanmean(seg_sd[ok]) / c4 / np.nanmean(rms_ref[ok]))
                if ok.any() else None,
                observed_over_parametric=float(np.nanmean(seg_sd[ok]) / c4 / np.nanmean(rms_par[ok]))
                if ok.any() else None,
                observed_over_parametric_kappa=float(np.nanmean(seg_sd[ok]) / c4 / np.nanmean(rms_kap[ok]))
                if ok.any() else None)
            oka = (ok if rep == "intensity" else np.isfinite(seg_sd)) & np.isfinite(rms_ana)
            out[f"{w:g}|{sp}"].update(
                n_pairs_analytic=int(oka.sum()),
                observed_over_analytic=float(np.nanmean(seg_sd[oka]) / c4 / np.nanmean(rms_ana[oka]))
                if oka.any() else None)
            print(w, sp, out[f"{w:g}|{sp}"], flush=True)
    (RESULTS / "noisecal").mkdir(parents=True, exist_ok=True)
    (RESULTS / "noisecal" / f"{design}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
