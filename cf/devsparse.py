"""Development (cryo, in-sample, not adopted): cf/proc7.py's sparse-regime detector where v6's detector is
idle. Usage: CF_DATA=... python -m cf.devsparse SHARD NSHARDS ALPHA -> results/dev/cryo/p7_sparse/."""

from cf.result_io import open_result, glob_results

import glob, json, os, sys, time
import numpy as np
from cf.analyze import load_results, score_clusters
from cf.fastcc import trial_arrays
from cf.newproc import load_module, carve, N_BOOTSTRAP, SEED
from cf.paths import DATA
from cf.proc7 import sparse_groups

g = load_module("v6")
alpha = float(sys.argv[3])
v6 = load_results("cryo", "v6"); today = load_results("cryo", "today")
files = sorted(glob_results(str(DATA / "trials" / "cryo" / "*.npz")))[int(sys.argv[1])::int(sys.argv[2])]
from cf.paths import RESULTS
od = RESULTS / "dev" / "cryo" / "p7_sparse"
od.mkdir(parents=True, exist_ok=True)
out = open_result(od / f"shard_{int(sys.argv[1]):03d}.jsonl", "w")
for f in files:
    tid = os.path.basename(f)[:-4]
    r6 = v6.get(tid)
    if r6 is None or "error" in r6:
        continue
    if len(r6.get("detector_log", [])) > 0:  # v6's detector was active: v6's result stands
        out.write(json.dumps(dict(id=tid, procedure="p7_sparse", clusters=r6["clusters"], sparse=None)) + "\n"); continue
    d = np.load(f); a = trial_arrays(d)
    rij = a["rij"].copy(); np.fill_diagonal(rij, 0)
    na = g.noise_replicates(a["y"], a["v"], a["m"], rij, model="analytic", n_replicates=N_BOOTSTRAP, seed=SEED)
    grp, log = sparse_groups(rij, a["wij"], na, alpha=alpha)
    clusters = r6["clusters"]
    if grp is not None:
        names = a["names"]
        G = [names[i] for i in grp]
        rest = [x for x in names if x not in set(G)]
        clusters = carve(today[tid].get("clusters", []), [rest, G], [str(x) for x in d["names"]], "half")
    out.write(json.dumps(dict(id=tid, procedure="p7_sparse", clusters=clusters, sparse=log)) + "\n")
