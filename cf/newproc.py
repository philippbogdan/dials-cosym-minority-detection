"""The frozen procedures (PREREGISTRATION.md amendments 3 and 4) on every trial of a design,
from the arrays saved by the today's-DIALS run.

The detector code is the DIALS change itself, numpy only: `v3` (the proposed procedure,
amendment 6) and its frontier variants run external/dials/src/dials/algorithms/correlation/
groups.py (dials-patch/0001-0004); `v2` and `v1` run the module as frozen in amendments 5 and
3 (cf/frozen/groups_v2.py = dials-patch/0001-0003, cf/frozen/groups_v1.py = dials-patch/0001).
Every configuration sees the same
rij/wij and the same merged intensities DIALS uses (same dense layout), and the same
counting-noise replicates: the bootstrap over reflections with the seed and number of
replicates of `CorrelationMatrix.find_small_groups` (and, for the `v2_analytic` frontier
point, the analytic replicates of cf.bootstrap). Groups other than the largest are taken out
of today's OPTICS clusters of the trial exactly as the patch does. End-to-end runs of the
patched DIALS on a sample of trials check the equivalence (cf.run_trials ... groups).

Output: results/trials/<design>/<config>/shard_<i>.jsonl for every configuration in CONFIGS.
"""

from __future__ import annotations

from cf.result_io import open_result

import importlib.util
import json
import os
import sys
import time

import numpy as np

from cf.analyze import load_design, load_results
from cf.fastcc import trial_arrays
from cf.paths import DATA, RESULTS, ROOT, base_design, representation

N_BOOTSTRAP, SEED = 100, 0
MIN_GROUP_SIZE = {"v7": 5, "v8": 5, "v9": 5}  # groups.min_group_size of the configurations that set it
# groups.small_groups of v8: max_se, z, min_coherence, min_ratio; of v9 (amendment 22): the
# arguments of its module's small_group_admitted
SMALL_GROUPS = {"v8": (0.002, 2.0, 0.04, 2.0),
                "v9": dict(min_signal_to_noise=30.0, z=2.0, min_coherence=0.01, min_ratio=1.5, max_se=None)}
SMALL_GROUPS["v10"] = SMALL_GROUPS["v9"]
# v11 (amendment 26): either v9's condition or v8's (small_groups and small_groups.precise),
# through its module's admitted_by_any
SMALL_GROUPS["v11"] = [SMALL_GROUPS["v9"], dict(min_signal_to_noise=0.0, z=2.0, min_coherence=0.04, min_ratio=2.0,
                                                max_se=0.002)]
# v10 (amendment 23): groups smaller than this are not taken out, decided on all of the detector's
# groups (take_out_groups(min_group_size, keep)) instead of being merged into the largest first
TAKE_OUT_MIN = {"v10": 5, "v11": 5}
V3 = dict(significance=0.01, separation_ratio=3.5, min_population=5, min_coherence=0.05)
V4 = dict(V3, pair_significance=0.001, min_usable_pairs=3)
CONFIGS = {
    # name: (module, detector parameters, noise model[, carve rule])
    "v6": ("v6", V4, "both", "half"),
    # the change as proposed (dials-patch/0001-0009, external/dials): v6's decisions; 0008 fixes
    # DIALS's logging, 0009 computes the bootstrap by weighted sums (its cost is measured here)
    "v6_1": ("live", V4, "both", "half"),
    # v7 (amendment 17, after the room-temperature judgement of v6): the change as proposed
    # (dials-patch/0001-0010) with its defaults: v6's groups, those smaller than 5 (OPTICS's own
    # smallest cluster) returned to the background population
    "v7": ("v7", V4, "both", "half"),
    # v8 (amendment 20, after the room-temperature judgement of v6): v7, and groups smaller than five
    # when their coherence is precisely measured and clearly above the floor (small_groups)
    "v8": ("v8", V4, "both", "half"),
    # v9 (amendment 22): v7, and groups smaller than five whose coherence stands 30 standard errors
    # out of the counting noise (small_groups by signal to noise; dials-patch/0001-0013)
    "v9": ("v9", V4, "both", "half"),
    # v10 (amendment 23): v9, the take-out decided on all groups before the minimum group size
    # (dials-patch/0001-0014)
    "v10": ("v10", V4, "both", "half"),
    # v11 (amendment 26): v10, a small group admitted by v9's condition or by v8's (dials-patch/0001-0015)
    "v11": ("v11", V4, "both", "half"),
    "v5": ("v5", V4, "analytic", "half"),
    "v4": ("v4", V4, "analytic"),
    "v3": ("v3", V3, "bootstrap"),
    "v2": ("v2", V3, "bootstrap"),  # min_coherence 0.05 for every size (amendment 5)
    "v1": ("v1", dict(significance=0.01, separation_ratio=3.5, min_population=5), "bootstrap"),
    "v4_ratio_only": ("v4", dict(V4, min_coherence=None), "analytic"),
    "v4_coherence_only": ("v4", dict(V4, separation_ratio=1e9), "analytic"),
    "v4_bootstrap": ("v4", V4, "bootstrap"),
    "v3_ratio_only": ("v3", dict(V3, min_coherence=None), "bootstrap"),
    "v3_coherence_only": ("v3", dict(V3, separation_ratio=1e9), "bootstrap"),
    "v3_analytic": ("v3", V3, "analytic"),
}
# development variants (cryo only, run on request: python -m cf.newproc cryo v2_f0.04,...)
DEV = {f"v2_f{f:g}": ("v2", dict(V3, min_coherence=f), "bootstrap") for f in (0.035, 0.04, 0.045, 0.05)}
# after v5 (not adopted): the scan statistic "both" (gap or mean within), and a weaker guard
DEV.update({
    "v6dev": ("v6dev", V4, "analytic", "half"),
    "v6dev_gap": ("v6dev", dict(V4, statistic="gap"), "analytic", "half"),
    "v6dev_guard2": ("v6dev", dict(V4, min_usable_pairs=2), "analytic", "half"),
    # the scan's internal-coherence statistic as a standardised sum over usable pairs
    "v6dev2": ("v6dev2", V4, "analytic", "half"),
    "v6dev2_sum": ("v6dev2", dict(V4, statistic="sum"), "analytic", "half"),
    # candidates chosen by the gap, tested by their standardised block sum
    "v6dev3": ("v6dev3", dict(V4, statistic="select"), "analytic", "half"),
})


def load_module(which):
    path = (ROOT / "cf" / "frozen" / f"groups_{which}.py" if which in ("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v6dev", "v6dev2", "v6dev3") else
            ROOT / "cf" / "frozen" / "groups_v6_1.py")
    name = f"groups_{which}"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def patched_groups_module():
    return load_module("v6")


def carve(optics_clusters, groups, names, rule="largest", min_group_size=None, keep=None):
    """As CorrelationMatrix.find_small_groups: groups other than the largest become clusters;
    every other data set keeps its OPTICS label; clusters of one become outliers. With
    rule="half" (v5, dials.algorithms.correlation.groups.take_out_groups) every group becomes
    a cluster when no group holds more than half of the data sets."""
    label = {n: -1 for n in names}
    for i, c in enumerate(optics_clusters):
        for x in c:
            label[x] = i
    nxt = len(optics_clusters)
    ordered = sorted(groups, key=len, reverse=True)
    carved = ordered[1:]
    if rule == "half" and ordered and 2 * len(ordered[0]) <= len(names):
        carved = ordered
    if min_group_size is not None:  # v10: decided on all groups, then the small ones not admitted stay
        carved = [g for g in carved if len(g) >= min_group_size or (keep is not None and keep(g))]
    for g in carved:
        for x in g:
            label[x] = nxt
        nxt += 1
    out = {}
    for x, lab in label.items():
        if lab >= 0:
            out.setdefault(lab, []).append(x)
    return [v for v in out.values() if len(v) >= 2]


def main(design, shard, nshards, configs=None):
    configs = configs or list(CONFIGS)
    CONFIGS.update({c: DEV[c] for c in configs if c in DEV})
    mods = {w: load_module(w) for w in {CONFIGS[c][0] for c in configs} | {"v6"}}
    trials = load_design(design)
    today = load_results(design, "today")
    ids = sorted(trials)[shard::nshards]
    outs, done = {}, {}
    for c in configs:
        od = RESULTS / "trials" / design / c
        od.mkdir(parents=True, exist_ok=True)
        outs[c] = od / f"shard_{shard:03d}.jsonl"
        done[c] = {json.loads(l)["id"] for l in open_result(outs[c])} if outs[c].exists() else set()
    for tid in ids:
        todo = [c for c in configs if tid not in done[c]]
        if not todo:
            continue
        f = DATA / "trials" / base_design(design) / f"{tid}.npz"
        r_today = today.get(tid)
        recs = {c: dict(id=tid, procedure=c) for c in todo}
        if r_today is None or not f.exists():
            for c in todo:
                recs[c]["error"] = "missing today's result or arrays"
        else:
            d = np.load(f)
            all_names = [str(x) for x in d["names"]]  # every crystal of the trial DIALS used
            optics = r_today.get("clusters", []) if "error" not in r_today else []
            t0 = time.process_time()
            a = trial_arrays(d, representation(design))  # rows of DIALS's matrix (cosym may drop a crystal)
            names, rij, wij, y, v, m = a["names"], a["rij"], a["wij"], a["y"], a["v"], a["m"]
            np.fill_diagonal(rij, 0.0)
            dense_cpu = time.process_time() - t0
            noise, noise_cpu = {}, {}
            # the noise models of the DIALS change itself: v6's module for the frozen
            # procedures, the module as proposed for v6_1
            needs = set()
            for c in todo:
                src = CONFIGS[c][0] if CONFIGS[c][0] in ("live", "v7", "v8", "v9", "v10", "v11") else "v6"
                for model in (("analytic", "bootstrap") if CONFIGS[c][2] == "both" else (CONFIGS[c][2],)):
                    needs.add((src, model))
            for src, model in sorted(needs):
                t0 = time.process_time()
                noise[(src, model)] = mods[src].noise_replicates(y, v, m, rij, model=model,
                                                                 n_replicates=N_BOOTSTRAP, seed=SEED)
                noise_cpu[(src, model)] = time.process_time() - t0
            for c in todo:
                which, params, model = CONFIGS[c][:3]
                rule = CONFIGS[c][3] if len(CONFIGS[c]) > 3 else "largest"
                src = which if which in ("live", "v7", "v8", "v9", "v10", "v11") else "v6"
                t0 = time.process_time()
                take_out = {}
                if model == "both":  # v6: the analytic and the bootstrap detector, combined
                    g = mods[which]
                    det = g.GroupDetection(rij, wij, noise[(src, "analytic")], **params)
                    det_b = g.GroupDetection(rij, wij, noise[(src, "bootstrap")],
                                             **dict(params, pair_significance=None, min_usable_pairs=0))
                    combined = g.combine_groups(det.run(), det_b.run(), len(names))
                    if c in MIN_GROUP_SIZE:
                        keep = None
                        if c in SMALL_GROUPS and isinstance(SMALL_GROUPS[c], dict):
                            ana, sg = noise[(src, "analytic")], SMALL_GROUPS[c]

                            def keep(grp, ana=ana, sg=sg, g=g):
                                return g.small_group_admitted(*g.group_evidence(rij, wij, ana, grp), len(grp), **sg)
                        elif c in SMALL_GROUPS:
                            max_se, zz, c0, r0 = SMALL_GROUPS[c]
                            ana = noise[(src, "analytic")]

                            def keep(grp, ana=ana, max_se=max_se, zz=zz, c0=c0, r0=r0):
                                coherence, se, d = g.group_evidence(rij, wij, ana, grp)
                                return bool((max_se is None or se <= max_se) and coherence - zz * se >= c0 * 2 / len(grp)
                                            and (d <= 0 or coherence / d >= r0))
                        # v7's frozen module predates keep (0011)
                        combined = (g.apply_min_group_size(combined, MIN_GROUP_SIZE[c], keep=keep) if keep
                                    else g.apply_min_group_size(combined, MIN_GROUP_SIZE[c]))
                    if c in TAKE_OUT_MIN:  # v10: the small groups its rule admits, by their names
                        ana, sg = noise[(src, "analytic")], SMALL_GROUPS[c]
                        def admit(grp, ana=ana, sg=sg, g=g):  # (not "rule": that is the take-out rule)
                            ev = g.group_evidence(rij, wij, ana, grp)
                            return (g.admitted_by_any(*ev, len(grp), sg) if isinstance(sg, list)
                                    else g.small_group_admitted(*ev, len(grp), **sg))
                        admitted = {tuple(sorted(names[i] for i in grp)) for grp in combined
                                    if len(grp) < TAKE_OUT_MIN[c] and admit(grp)}
                        take_out = dict(min_group_size=TAKE_OUT_MIN[c],
                                        keep=lambda ng, admitted=admitted: tuple(sorted(ng)) in admitted)
                    det.log = det.log + det_b.log
                    groups = [[names[i] for i in grp] for grp in combined]
                else:
                    det = mods[which].GroupDetection(rij, wij, noise[(src, model)], **params)
                    groups = [[names[i] for i in grp] for grp in det.run()]
                cpu = time.process_time() - t0
                recs[c].update(clusters=carve(optics, groups, all_names, rule, **take_out), detector_groups=groups,
                               detector_log=det.log, detector_cpu_s=cpu,
                               noise_cpu_s=dense_cpu + (noise_cpu[(src, "analytic")] + noise_cpu[(src, "bootstrap")]
                                                         if model == "both" else noise_cpu[(src, model)]),
                               optics_cpu_s=r_today.get("cpu_s"))
        for c in todo:
            with open_result(outs[c], "a") as fh:
                fh.write(json.dumps(recs[c], default=float) + "\n")


if __name__ == "__main__":
    cfg = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    main(sys.argv[1], int(os.environ.get("CF_SHARD", 0)), int(os.environ.get("CF_NSHARDS", 1)), cfg)
