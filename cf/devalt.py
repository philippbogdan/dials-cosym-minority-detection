"""Development (cryo, in-sample; amendment 14): the v6 detector on pairwise matrices computed from
the same reflections in other ways (amplitudes, weighted amplitudes, normal scores, normalised
amplitudes), with v6's parameters or other coherence thresholds. Run with CF_DATA pointing at the
arrays: python -m cf.devalt VARIANT SHARD NSHARDS [PARAMS]. Writes results/dev/cryo/alt_<variant>[_<params>]/."""

from cf.result_io import open_result, glob_results

import glob, json, os, sys, time
import numpy as np
from cf.analyze import load_design, load_results
from cf.fastcc import trial_arrays
from cf.newproc import load_module, carve, V4, N_BOOTSTRAP, SEED
from cf.paths import DATA, RESULTS

g = load_module("v6")

HKL = {}

def eamp(y, m, nbins=10):
    """|E|: amplitudes of intensities normalised by the trial's mean intensity in shells of
    h^2 + k^2 + l^2 (resolution, cubic cell), equal-count shells."""
    s2 = HKL["s2"]
    edges = np.quantile(s2, np.linspace(0, 1, nbins + 1))
    b = np.clip(np.searchsorted(edges, s2, side="right") - 1, 0, nbins - 1)
    mean = np.array([y[:, b == j][m[:, b == j]].mean() if (m[:, b == j]).any() else 1.0 for j in range(nbins)])
    e2 = y / mean[b][None, :]
    return np.sqrt(np.maximum(e2, 0.0))

def transform(y, v, m, variant):
    if variant == "dials":
        return y, v
    if variant == "unw":
        return y, np.where(m, 1.0, 0.0)
    if variant == "amp":
        return np.sqrt(np.maximum(y, 0.0)), np.where(m, 1.0, 0.0)
    if variant == "ampw":
        F = np.sqrt(np.maximum(y, 0.0))
        sF2 = v / (4.0 * np.maximum(y, np.sqrt(np.maximum(v, 1e-300))))
        return F, np.where(m, sF2, 0.0)
    if variant == "eamp":
        return eamp(y, m), np.where(m, 1.0, 0.0)
    if variant == "nscore":
        from scipy.stats import norm, rankdata
        out = np.zeros_like(y)
        for i in range(y.shape[0]):
            idx = np.where(m[i])[0]
            if len(idx):
                out[i, idx] = norm.ppf(rankdata(y[i, idx]) / (len(idx) + 1))
        return out, np.where(m, 1.0, 0.0)
    raise ValueError(variant)

PARAMS = {"ratio": dict(V4, min_coherence=None), "c0.1": dict(V4, min_coherence=0.1), "c0.08": dict(V4, min_coherence=0.08)}

def run(variant, files, shard, nshards, ptag=None):
    design = "cryo"
    today = load_results(design, "today")
    od = RESULTS / "dev" / design / (f"alt_{variant}" + (f"_{ptag}" if ptag else ""))
    P = PARAMS[ptag] if ptag else V4
    od.mkdir(parents=True, exist_ok=True)
    out = od / f"shard_{shard:03d}.jsonl"
    done = {json.loads(l)["id"] for l in open_result(out)} if out.exists() else set()
    for f in files[shard::nshards]:
        tid = os.path.basename(f)[:-4]
        if tid in done or tid not in today or "error" in today[tid]:
            continue
        d = np.load(f)
        all_names = [str(x) for x in d["names"]]
        a = trial_arrays(d)
        names, y, v, m = a["names"], a["y"], a["v"], a["m"]
        hkl = np.asarray(d["hkl"])[np.asarray(d["eps1"])] if "eps1" in d.files else np.asarray(d["hkl"])
        key = (hkl[:, 0].astype(np.int64) * 4096 + hkl[:, 1]) * 4096 + hkl[:, 2]
        uniq, first = np.unique(key, return_index=True)
        HKL["s2"] = (hkl[first] ** 2).sum(1).astype(float)
        y2, v2 = transform(y, v, m, variant)
        if variant == "dials":
            rij, wij = a["rij"], a["wij"]
        else:
            rij = g.weighted_cc_matrix(y2, v2, m)
            wij = g.effective_sizes(y2, v2, m) * (rij != 0)
        np.fill_diagonal(rij, 0.0)
        t0 = time.process_time()
        na = g.noise_replicates(y2, v2, m, rij, model="analytic", n_replicates=N_BOOTSTRAP, seed=SEED)
        nb = g.noise_replicates(y2, v2, m, rij, model="bootstrap", n_replicates=N_BOOTSTRAP, seed=SEED)
        det = g.GroupDetection(rij, wij, na, **P)
        det_b = g.GroupDetection(rij, wij, nb, **dict(P, pair_significance=None, min_usable_pairs=0))
        combined = g.combine_groups(det.run(), det_b.run(), len(names))
        groups = [[names[i] for i in grp] for grp in combined]
        rec = dict(id=tid, procedure=f"v6_{variant}", clusters=carve(today[tid].get("clusters", []), groups, all_names, "half"),
                   detector_groups=groups, detector_log=det.log + det_b.log, cpu_s=time.process_time() - t0)
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")

if __name__ == "__main__":
    variant, shard, nshards = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    ptag = sys.argv[4] if len(sys.argv) > 4 else None
    files = sorted(glob_results(str(DATA / "trials" / "cryo" / "*.npz")))
    run(variant, files, shard, nshards, ptag)
