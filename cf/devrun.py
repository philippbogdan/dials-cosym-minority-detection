"""Run a matrix-level procedure configuration on every trial of a design (sharded), from the
stored matrices (results/cmat) and bootstrap replicates (results/boot). Output:
results/dev/<design>/<tag>/shard_<i>.jsonl with the reported groups (crystal names) and cost.

Usage: python -m cf.devrun DESIGN TAG '{"alpha":0.01,"c":6.0}' [--kinds single,mixture]
"""

from __future__ import annotations

from cf.result_io import open_result

import argparse
import json
import os
import time

import numpy as np

from cf.analyze import load_design
from cf.bootstrap import load as bload
from cf.cmatrices import load
from cf.paths import RESULTS
from cf import proc5, proc6, proc7
from cf.proc5 import Data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("design")
    ap.add_argument("tag")
    ap.add_argument("config")
    ap.add_argument("--shard", type=int, default=int(os.environ.get("CF_SHARD", 0)))
    ap.add_argument("--nshards", type=int, default=int(os.environ.get("CF_NSHARDS", 1)))
    a = ap.parse_args()
    cfg = json.loads(a.config)
    proc = {"p5": proc5.analyze, "p6": proc6.analyze, "p7": proc7.analyze,
            "p6u": proc6.analyze_usable}[cfg.pop("proc", "p5")]
    trials = load_design(a.design)
    ids = sorted(t for t in trials if trials[t]["kind"] == "single" or trials[t]["k"] >= 2)
    ids = ids[a.shard :: a.nshards]
    mats = load(a.design)
    boots = bload(a.design, cfg.pop("noise", "reflections"))
    # noise_scale: {wedge: factor} multiplies the replicate noise (e.g. the repeat-measurement
    # calibration of cf.noisecal); an analysis variant that knows w, not a deployable procedure
    noise_scale = {float(k): v for k, v in cfg.pop("noise_scale", {}).items()}
    od = RESULTS / "dev" / a.design / a.tag
    od.mkdir(parents=True, exist_ok=True)
    with open_result(od / "config.json", "w") as fh:
        json.dump(json.loads(a.config), fh)
    out = od / f"shard_{a.shard:03d}.jsonl"
    done = {json.loads(l)["id"] for l in open_result(out)} if out.exists() else set()
    for tid in ids:
        if tid in done or tid not in mats or tid not in boots:
            continue
        m = mats[tid]
        names = [str(x) for x in m["names"]]
        D = Data(m["rij"].astype(float), m["wij"].astype(float), boots[tid].astype(float))
        t0 = time.process_time()
        kw = dict(cfg)
        if noise_scale:
            kw["calibrate"] = noise_scale.get(float(trials[tid]["w"]), 1.0)
        groups = proc(D, **kw)
        cpu = time.process_time() - t0
        rec = dict(id=tid, clusters=[[names[i] for i in g] for g in groups if len(g) >= 2],
                   cpu_s=cpu)
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
