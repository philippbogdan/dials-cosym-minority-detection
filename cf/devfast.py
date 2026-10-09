"""Cost (after v6): the detector's groups with the reflection bootstrap computed per replicate
(frozen v6) and by weighted sums (the module now), on every trial with arrays; groups must be
identical. Usage: python -m cf.devfast DESIGN SHARD NSHARDS -> results/dev/<design>/fastboot/."""

from cf.result_io import open_result, glob_results

import glob, json, os, sys, time
import numpy as np
from cf.fastcc import trial_arrays
from cf.newproc import load_module, V4, N_BOOTSTRAP, SEED
from cf.paths import DATA, RESULTS, base_design

old, new = load_module("v6"), load_module("live")


def groups(g, rij, wij, na, nb, n):
    det = g.GroupDetection(rij, wij, na, **V4)
    det_b = g.GroupDetection(rij, wij, nb, **dict(V4, pair_significance=None, min_usable_pairs=0))
    return sorted(sorted(int(i) for i in grp) for grp in g.combine_groups(det.run(), det_b.run(), n))


def main(design, shard, nshards):
    files = sorted(glob_results(str(DATA / "trials" / base_design(design) / "*.npz")))[shard::nshards]
    od = RESULTS / "dev" / design / "fastboot"
    od.mkdir(parents=True, exist_ok=True)
    out = od / f"shard_{shard:03d}.jsonl"
    done = {json.loads(l)["id"] for l in open_result(out)} if out.exists() else set()
    for f in files:
        tid = os.path.basename(f)[:-4]
        if tid in done:
            continue
        d = np.load(f)
        if "eps1" not in d.files:
            continue
        a = trial_arrays(d)
        y, v, m, rij, wij = a["y"], a["v"], a["m"], a["rij"].copy(), a["wij"]
        np.fill_diagonal(rij, 0.0)
        na = old.noise_replicates(y, v, m, rij, model="analytic", n_replicates=N_BOOTSTRAP, seed=SEED)
        t0 = time.process_time()
        nb_old = old.noise_replicates(y, v, m, rij, model="bootstrap", n_replicates=N_BOOTSTRAP, seed=SEED)
        t_old = time.process_time() - t0
        t0 = time.process_time()
        nb_new = new.noise_replicates(y, v, m, rij, model="bootstrap", n_replicates=N_BOOTSTRAP, seed=SEED)
        t_new = time.process_time() - t0
        g_old = groups(old, rij, wij, na, nb_old, len(a["names"]))
        g_new = groups(new, rij, wij, na, nb_new, len(a["names"]))
        rec = dict(id=tid, same=g_old == g_new, max_diff=float(np.abs(nb_old - nb_new).max()),
                   cpu_old=t_old, cpu_new=t_new)
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
