"""Development (cryo, in-sample, not adopted): cf/proc7.py sparse_groups_local where v6's detector is idle,
trials of 10 or more, up to three groups. Usage: CF_DATA=... python -m cf.devsparse_local SHARD NSHARDS
ALPHA NMIN ROUNDS -> results/dev/cryo/p7_sparse_local/."""

from cf.result_io import open_result, glob_results

import glob, json, os, sys
import numpy as np
from cf.analyze import load_results
from cf.fastcc import trial_arrays
from cf.newproc import load_module, carve, N_BOOTSTRAP, SEED
from cf.paths import DATA
from cf.proc7 import sparse_groups_local as sparse_groups

g = load_module("v6")
alpha, nmin, rounds = float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
v6 = load_results("cryo", "v6"); today = load_results("cryo", "today")
files = sorted(glob_results(str(DATA / "trials" / "cryo" / "*.npz")))[int(sys.argv[1])::int(sys.argv[2])]
from cf.paths import RESULTS
od = RESULTS / "dev" / "cryo" / "p7_sparse_local"
od.mkdir(parents=True, exist_ok=True)
out = open_result(od / f"shard_{int(sys.argv[1]):03d}.jsonl", "w")
for f in files:
    tid = os.path.basename(f)[:-4]
    r6 = v6.get(tid)
    if r6 is None or "error" in r6:
        continue
    d = np.load(f); a = trial_arrays(d)
    names = a["names"]
    if len(r6.get("detector_log", [])) > 0 or len(names) < nmin:
        out.write(json.dumps(dict(id=tid, procedure="p7_sparse_local", clusters=r6["clusters"], sparse=[])) + "\n"); continue
    rij = a["rij"].copy(); np.fill_diagonal(rij, 0)
    na = g.noise_replicates(a["y"], a["v"], a["m"], rij, model="analytic", n_replicates=N_BOOTSTRAP, seed=SEED)
    idx = np.arange(len(names)); groups, logs = [], []
    for _ in range(rounds):
        sub = np.ix_(idx, idx)
        grp, log = sparse_groups(rij[sub], a["wij"][sub], na[:, idx][:, :, idx], alpha=alpha)
        if log: logs.append(log)
        if grp is None:
            break
        groups.append([names[i] for i in idx[grp]])
        idx = np.setdiff1d(idx, idx[grp])
    clusters = r6["clusters"]
    if groups:
        taken = {x for gg in groups for x in gg}
        rest = [x for x in names if x not in taken]
        clusters = carve(today[tid].get("clusters", []), [rest] + groups, [str(x) for x in d["names"]], "half")
    out.write(json.dumps(dict(id=tid, procedure="p7_sparse_local", clusters=clusters, sparse=logs)) + "\n")
