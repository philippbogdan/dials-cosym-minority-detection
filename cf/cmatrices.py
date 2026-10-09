"""DIALS's pairwise correlation matrices (rij, wij) of every trial of a design, from the arrays
saved by cf.run_trials, stored compactly for analysis without the DIALS environment:
results/cmat/<design>/<shard>.npz with, per trial id, rij, wij (float32) and the crystal names
in matrix order (with CF_DEVIATION_CORRELATIONS set, also the superseded statistic of
amendment 1)."""

from __future__ import annotations

from cf.result_io import glob_results

import glob
import json
import os
import sys

import numpy as np

from cf.fastcc import trial_arrays
from cf.oracle import from_npz
from cf.paths import DATA, RESULTS, base_design, representation


def main(design: str, shard: int, nshards: int) -> None:
    files = sorted(glob_results(str(DATA / "trials" / base_design(design) / "*.npz")))[shard::nshards]
    out = {}
    for f in files:
        tid = os.path.basename(f)[:-4]
        d = np.load(f)
        if "eps1" not in d.files:
            continue
        if os.environ.get("CF_DEVIATION_CORRELATIONS"):  # amendment-1 statistic, superseded
            c, npairs = from_npz(d)
            out[f"{tid}|c"] = c.astype(np.float32)
            out[f"{tid}|n"] = npairs.astype(np.int32)
        a = trial_arrays(d, representation(design))
        out[f"{tid}|names"] = np.array(a["names"])  # aligned with rij
        out[f"{tid}|rij"] = np.array(a["rij"], dtype=np.float32)
        out[f"{tid}|wij"] = np.array(a["wij"], dtype=np.float32)
    od = RESULTS / "cmat" / design
    od.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(od / f"shard_{shard:03d}.npz", **out)
    print(len(files), "trials")


def load(design: str) -> dict:
    res = {}
    for f in sorted(glob_results(str(RESULTS / "cmat" / design / "shard_*.npz"))):
        d = np.load(f)
        for key in d.files:
            tid, what = key.split("|")
            res.setdefault(tid, {})[what] = d[key]
    return res


if __name__ == "__main__":
    shard = int(os.environ.get("CF_SHARD", 0))
    nshards = int(os.environ.get("CF_NSHARDS", 1))
    main(sys.argv[1], shard, nshards)
