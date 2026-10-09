"""Run the trials of a design file with one procedure. Shardable and resumable.

Each trial: select crystals, cut to w degrees, scale with dials.scale defaults, run the
procedure, append one JSON line to results/trials/<design>/<procedure>/shard_<i>.jsonl.
With --save-arrays, the inputs every intensity-based procedure shares (the per-crystal
merged, normalised, resolution-filtered intensities that DIALS computes rij from, plus rij,
wij and the cosym coordinates) are kept in .venv-data/trials/<design>/<id>.npz.
"""

from __future__ import annotations

from cf.result_io import open_result

import argparse
import json
import os
import time
import traceback

import numpy as np

from cf import design, today, trialdata
from cf.paths import DATA, RESULTS

FULL = {"cryo": 10.0, "rt": 60.0, "rtbh": 20.0}
_cc = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".cf_commit")
CODE_COMMIT = open_result(_cc).read().strip() if os.path.exists(_cc) else None


def cosym_inputs(cm):
    """The merged normalised intensities used for rij, in the target's primitive ASU."""
    ca = cm.cosym_analysis
    data = ca.intensities.customized_copy(anomalous_flag=False)
    cb = data.change_of_basis_op_to_primitive_setting()
    data = data.change_basis(cb).map_to_asu()
    hkl = np.array(data.indices(), dtype=np.int32).reshape(-1, 3)
    # DIALS's rij matcher only pairs reflections with epsilon == 1 in the Patterson group
    pg = ca.target._patterson_group
    eps1 = np.array(pg.epsilon(data.indices()), dtype=np.int32) == 1
    return dict(
        hkl=hkl,
        eps1=eps1,
        i=np.array(data.data(), dtype=np.float64),
        sig=np.array(data.sigmas(), dtype=np.float64),
        ds=np.array(ca.dataset_ids, dtype=np.int32),
    )


def run_one(coll_data, t, procedure, save_dir=None, nproc=1):
    names = design.trial_crystals(t)
    w = None if t["w"] >= FULL[t["coll"]] else t["w"]
    t0 = time.time()
    expts, tables = coll_data.trial(names, w, segment=t.get("seg", 0))
    n_refl_in = [len(x) for x in tables]
    sexpts, joint, scale_cpu = trialdata.scale(expts, tables, nproc=nproc)
    if procedure == "today":
        res = today.run(sexpts, joint, names, nproc=nproc)
    elif procedure == "groups":
        # the proposed DIALS change (run with pipeline/run_py_patched.sh), configured by the
        # frozen parameters in dials-patch/frozen.phil
        from dials.util.phil import parse

        params = today.default_params()
        frozen = os.environ.get("FROZEN_PHIL") or os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dials-patch", "frozen.phil")
        from dials.command_line.correlation_matrix import phil_scope

        params = phil_scope.fetch(sources=[parse(open_result(frozen).read())]).extract()
        assert params.significant_clusters.find_small_groups
        res = today.run(sexpts, joint, names, params=params, nproc=nproc)
        # for the equivalence check only: the detector's own groups, recomputed exactly as
        # CorrelationMatrix.find_small_groups computes them (same inputs, same seed)
        from dials.algorithms.correlation.groups import (
            GroupDetection,
            combine_groups,
            dense_from_target,
            noise_replicates,
        )

        cm = res["_cm"]
        gp = params.significant_clusters.groups
        target = cm.cosym_analysis.target
        y, v, m = dense_from_target(target)
        rij = np.array(target.rij_matrix, dtype=np.float64)
        np.fill_diagonal(rij, 0.0)
        wij = np.array(target.wij_matrix, dtype=np.float64)
        models = ("analytic", "bootstrap") if gp.noise_model == "both" else (gp.noise_model,)
        found = []
        for model in models:  # as CorrelationMatrix.find_small_groups
            det = GroupDetection(rij, wij, noise_replicates(y, v, m, rij, model=model,
                                                           n_replicates=gp.n_replicates, seed=gp.seed),
                                 significance=gp.significance, separation_ratio=gp.separation_ratio,
                                 min_population=gp.min_population, min_coherence=gp.min_coherence,
                                 pair_significance=gp.pair_significance if model == "analytic" else None,
                                 min_usable_pairs=gp.min_usable_pairs if model == "analytic" else 0)
            found.append(det.run())
        groups = found[0] if len(found) == 1 else combine_groups(found[0], found[1], len(rij))
        # rows of cosym's matrix: the data sets present in its arrays (it may drop one)
        present = sorted(set(np.asarray(target._lattice_ids).tolist()))
        used = [res["used"][i] for i in present] if len(present) != len(res["used"]) else res["used"]
        res["detector_groups"] = [[used[i] for i in g] for g in groups]
    else:
        raise ValueError(procedure)
    cm = res.pop("_cm")
    arrays = res.pop("_arrays", None)
    if save_dir is not None and arrays is not None:
        ci = cosym_inputs(cm)
        np.savez_compressed(
            save_dir / f"{t['id']}.npz", names=np.array(res["used"]), **arrays, **ci,
            dim_curve=np.array(res["dim_curve"]),
        )
    res.update(
        id=t["id"], procedure=procedure, scale_cpu_s=scale_cpu, code_commit=CODE_COMMIT,
        total_wall_s=time.time() - t0, n_refl_in=n_refl_in,
    )
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("design")
    ap.add_argument("procedure")
    ap.add_argument("--shard-base", type=int, default=0,
                    help="added to CF_SHARD (to split one sharding across places)")
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--nshards", type=int, default=int(os.environ.get("CF_NSHARDS", 1)))
    ap.add_argument("--only", default=None, help="comma-separated trial ids")
    ap.add_argument("--save-arrays", action="store_true")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--retry-errors", action="store_true",
                    help="rerun trials whose recorded result is an error (latest line wins)")
    args = ap.parse_args()
    if args.shard is None:
        args.shard = args.shard_base + int(os.environ.get("CF_SHARD", 0))
    dname = os.path.basename(args.design).replace(".jsonl", "")
    trials = [json.loads(l) for l in open_result(args.design)]
    if args.only:
        keep = set(args.only.split(","))
        trials = [t for t in trials if t["id"] in keep]
    mine = trials[args.shard :: args.nshards]
    outdir = RESULTS / "trials" / dname / (args.tag or args.procedure)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"shard_{args.shard:03d}.jsonl"
    done = set()
    if out.exists():
        for line in open_result(out):
            r = json.loads(line)
            if args.retry_errors and "error" in r:
                done.discard(r["id"])
            else:
                done.add(r["id"])
    save_dir = None
    if args.save_arrays:
        save_dir = DATA / "trials" / dname
        save_dir.mkdir(parents=True, exist_ok=True)
    colls = {}
    for t in mine:
        if t["id"] in done:
            continue
        cd = colls.setdefault(t["coll"], trialdata.Collection(t["coll"]))
        try:
            rec = run_one(cd, t, args.procedure, save_dir)
        except Exception as e:  # recorded, never silently dropped
            rec = dict(id=t["id"], procedure=args.procedure, error=repr(e),
                       traceback=traceback.format_exc()[-2000:])
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
