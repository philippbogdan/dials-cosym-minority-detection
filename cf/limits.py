"""Oracle and decoy for every mixture trial (k >= 2), PREREGISTRATION.md amendments 1-2.

For each trial: S(M) of the true minority on the rank-1 residuals of DIALS's rij, its oracle
p-value among the size-k subsets of the majority, and the largest S over those subsets (the
decoy). Writes results/limits/<design>.jsonl."""

from __future__ import annotations

from cf.result_io import open_result

import json
import sys

import numpy as np

from cf.analyze import load_design
from cf.cmatrices import load
from cf.oracle import decoy_max_resid, oracle_p_resid, residuals
from cf.paths import RESULTS


def run(design: str):
    trials = load_design(design)
    mats = load(design)
    out = RESULTS / "limits" / f"{design}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open_result(out, "w") as fh:
        for tid, t in trials.items():
            if t["kind"] != "mixture" or t["k"] < 2 or tid not in mats:
                continue
            m = mats[tid]
            names = [str(x) for x in m["names"]]
            pos = {n: i for i, n in enumerate(names)}
            minority = [pos[x] for x in t["minority"] if x in pos]
            majority = [pos[x] for x in t["majority"] if x in pos]
            rec = dict(id=tid, k=t["k"], n_minority_used=len(minority))
            if len(minority) < 2:
                rec["note"] = "fewer than 2 minority crystals in DIALS's matrix"
            else:
                e, rho = residuals(m["rij"].astype(float), m["wij"].astype(float))
                o = oracle_p_resid(e, minority, majority)
                dmax = decoy_max_resid(e, majority, len(minority))
                rec.update(S=o["S"], p=o["p"], exact=o["exact"], null_q95=o["null_q95"],
                           decoy=dmax["max"], decoy_exact=dmax["exact"],
                           decoy_members=[names[i] for i in dmax["members"]],
                           oracle_rejects=o["p"] < 0.05,
                           beats_decoy=o["S"] > dmax["max"])
                rec["within_noise"] = (not rec["oracle_rejects"]) or (not rec["beats_decoy"])
            fh.write(json.dumps(rec) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    run(sys.argv[1])
