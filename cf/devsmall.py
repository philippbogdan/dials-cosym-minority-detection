"""After v6's room-temperature judgement (amendment 20): the groups of two to four crystals v6
reports, measured as the procedure itself can (DIALS's analytic counting noise at face value, as v6's
analytic detector): coherence S (mean rank-1 residual within the group), its standard error over the
noise replicates (each refitted), and the strongest structure of the same size inside the rest (the
detector's decoy). Writes results/dev/<design>/smallgroups/shard_<i>.jsonl.
Usage: python -m cf.devsmall DESIGN SHARD NSHARDS"""

from cf.result_io import open_result

import json
import os
import sys

import numpy as np

from cf.analyze import load_design, load_results
from cf.fastcc import trial_arrays
from cf.newproc import load_module, N_BOOTSTRAP, SEED
from cf.paths import DATA, RESULTS

g = load_module("v6")


def main(design, shard, nshards):
    trials = load_design(design)
    v6 = load_results(design, "v6")
    od = RESULTS / "dev" / design / "smallgroups"
    od.mkdir(parents=True, exist_ok=True)
    out = od / f"shard_{shard:03d}.jsonl"
    done = {json.loads(l)["id"] for l in open_result(out)} if out.exists() else set()
    for tid in sorted(v6)[shard::nshards]:
        r = v6[tid]
        if tid in done or "error" in r:
            continue
        groups = sorted(r.get("detector_groups", []), key=len, reverse=True)
        small = [x for x in groups[1:] if 2 <= len(x) <= 4]
        rec = dict(id=tid, groups=[])
        if small:
            a = trial_arrays(np.load(DATA / "trials" / design / f"{tid}.npz"))
            names = a["names"]
            pos = {x: i for i, x in enumerate(names)}
            rij = a["rij"].copy()
            np.fill_diagonal(rij, 0.0)
            w = np.asarray(a["wij"], dtype=float)
            deltas = g.noise_replicates(a["y"], a["v"], a["m"], rij, model="analytic", n_replicates=N_BOOTSTRAP, seed=SEED)
            e, _ = g._residuals(rij, w)
            valid = (w > 0) & np.any(deltas != 0, axis=0)
            np.fill_diagonal(valid, False)
            eS = np.where(valid, e, np.nan)
            eb = []
            for d in deltas:
                rb = rij + d
                np.fill_diagonal(rb, 0.0)
                eb.append(np.where(valid, g._residuals(rb, w)[0], np.nan))
            for G in small:
                idx = np.array([pos[x] for x in G if x in pos])
                if len(idx) < 2:
                    continue
                rest = np.setdiff1d(np.arange(len(names)), idx)
                S = g._mean_within(eS, idx)
                reps = [g._mean_within(x, idx) for x in eb]
                reps = [v for v in reps if np.isfinite(v)]
                se = float(np.std(reps)) if len(reps) > 2 else float("nan")
                d = g.decoy(eS, [rest], len(idx))
                rec["groups"].append(dict(members=list(G), k=len(idx), S=float(S), se=se,
                                          decoy=float(d), ratio=float(S / d) if d > 0 else float("inf")))
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
