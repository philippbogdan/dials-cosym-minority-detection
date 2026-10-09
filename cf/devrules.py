"""The small-group rules (amendments 20, 22, 23), derived exactly from v6's groups: the procedures
that differ from v6 only in which groups of two to four they keep (v7: none; v8, v9, v10: those a
rule admits) are v6's detector groups (results/trials/<design>/v6) with the groups of two to four
measured as the procedure measures them (cf/devsmall.py: coherence S, its standard error se over
DIALS's analytic noise, the rest's strongest structure of the size, `decoy`). Two ways to treat the
groups a rule does not keep: merged into the largest population before the take-out (v7-v9), or left
in the OPTICS clusters with the take-out decided on all groups (v10, v11).

`python -m cf.devrules` searches the rule families of amendments 20 and 22 on cryo and room
temperature (in-sample) and writes results/dev/smallrules.json: per rule, the minorities found and
the single-population trials with a group where today's DIALS has none."""

from __future__ import annotations

from cf.result_io import open_result, glob_results

import glob
import itertools
import json

import numpy as np

from cf.analyze import load_design, load_results, score_clusters
from cf.newproc import carve
from cf.paths import RESULTS

COLLS = ("cryo", "rt")


def features(design):
    out = {}
    for f in glob_results(str(RESULTS / "dev" / design / "smallgroups" / "shard_*.jsonl")):
        for line in open_result(f):
            r = json.loads(line)
            out[r["id"]] = {tuple(sorted(g["members"])): g for g in r["groups"]}
    return out


def _take_out(optics_clusters, carved, names):
    """carve's steps for an explicit set of groups taken out."""
    label = {x: -1 for x in names}
    for i, c in enumerate(optics_clusters):
        for x in c:
            label[x] = i
    nxt = len(optics_clusters)
    for g in carved:
        for x in g:
            label[x] = nxt
        nxt += 1
    out = {}
    for x, lab in label.items():
        if lab >= 0:
            out.setdefault(lab, []).append(x)
    return [v for v in out.values() if len(v) >= 2]


def clusters(design, rule, take_out_first=False, min_group_size=5, cache={}):
    """Final clusters per trial of v6's groups with those of two to four kept by `rule(features)`."""
    if design not in cache:
        cache[design] = (load_design(design), load_results(design, "today"), load_results(design, "v6"),
                         features(design))
    trials, today, v6, feats = cache[design]
    out = {}
    for tid, r in v6.items():
        if "error" in r or tid not in today or "error" in today[tid]:
            continue
        t = trials[tid]
        names = t["crystals"] if "crystals" in t else t["majority"] + t["minority"]
        groups = sorted(r["detector_groups"], key=len, reverse=True)
        optics = today[tid].get("clusters", [])

        def admitted(g):
            if len(g) >= min_group_size:
                return True
            fe = feats.get(tid, {}).get(tuple(sorted(g)))
            return fe is not None and bool(rule(fe))

        if take_out_first:  # v10: v6's take-out on every group, then the small ones not admitted stay
            if not groups:
                out[tid] = carve(optics, [], names, "half")
                continue
            carved = groups if 2 * len(groups[0]) <= len(names) else groups[1:]
            out[tid] = _take_out(optics, [g for g in carved if admitted(g)], names)
        else:  # v7-v9: the small groups not admitted merged into the largest population first
            kept, back = [], list(groups[0]) if groups else []
            for g in groups[1:]:
                if admitted(g):
                    kept.append(g)
                else:
                    back += g
            out[tid] = carve(optics, [back] + kept if kept else [back], names, "half")
    return out


def counts(design, rule, take_out_first=False, cache={}):
    """Minorities found (k >= 2) and single-population trials with a group where today's DIALS has
    none."""
    rows = score_clusters(design, clusters(design, rule, take_out_first))
    if design not in cache:
        today = {k: v.get("clusters", []) for k, v in load_results(design, "today").items() if "error" not in v}
        cache[design] = {r["id"]: r for r in score_clusters(design, today)}
    st = cache[design]
    found = sum(bool(r.get("found")) for r in rows if r["kind"] == "mixture" and r.get("k", 0) >= 2)
    added = sum(1 for r in rows if r["kind"] == "single" and r["false_group"] and not st[r["id"]]["false_group"])
    return found, added


def rule_precise(max_se, z, c0, r0):
    """Amendment 20's family (v8 = 0.002, 2, 0.04, 2)."""
    return lambda g: (np.isfinite(g["se"]) and g["se"] <= max_se and g["S"] - z * g["se"] >= c0 * 2 / g["k"]
                      and (g["decoy"] <= 0 or g["S"] / g["decoy"] >= r0))


def rule_snr(z0, c0, r0):
    """Amendment 22's first family (v9 = 30, 0.01, 1.5)."""
    return lambda g: (np.isfinite(g["se"]) and g["se"] > 0 and g["S"] >= z0 * g["se"]
                      and g["S"] - 2 * g["se"] >= c0 * 2 / g["k"] and (g["decoy"] <= 0 or g["S"] / g["decoy"] >= r0))


def rule_zfloor(z, c0, scale, r0):
    """Amendment 22's second family: S - z se above a floor, no precision condition."""
    return lambda g: (np.isfinite(g["se"]) and g["S"] - z * g["se"] >= (c0 * 2 / g["k"] if scale else c0)
                      and (g["decoy"] <= 0 or g["S"] / g["decoy"] >= r0))


def main():
    res = dict(precise=[], snr=[], zfloor=[], chosen={})
    for max_se, z, c0, r0 in itertools.product((0.001, 0.002, 0.003, 0.005), (0, 1, 2, 3), (0.03, 0.04, 0.05, 0.06),
                                               (1.5, 2, 3, 5)):
        c = {coll: counts(coll, rule_precise(max_se, z, c0, r0)) for coll in COLLS}
        res["precise"].append(dict(max_se=max_se, z=z, c0=c0, r0=r0, **{f"{k}_{coll}": v for coll in COLLS
                                                                       for k, v in zip(("found", "added"), c[coll])}))
    for z0, c0, r0 in itertools.product((10, 15, 20, 25, 30, 40, 50), (0.0, 0.01, 0.02, 0.03, 0.04), (1.25, 1.5, 2, 3)):
        c = {coll: counts(coll, rule_snr(z0, c0, r0)) for coll in COLLS}
        res["snr"].append(dict(z0=z0, c0=c0, r0=r0, **{f"{k}_{coll}": v for coll in COLLS
                                                       for k, v in zip(("found", "added"), c[coll])}))
    for z, c0, scale, r0 in itertools.product((5, 8, 10, 12, 14, 16, 20), (0.01, 0.015, 0.02, 0.025, 0.03),
                                              (True, False), (1.25, 1.5, 2)):
        c = {coll: counts(coll, rule_zfloor(z, c0, scale, r0)) for coll in COLLS}
        res["zfloor"].append(dict(z=z, c0=c0, scale=scale, r0=r0, **{f"{k}_{coll}": v for coll in COLLS
                                                                     for k, v in zip(("found", "added"), c[coll])}))
    never = lambda g: False
    v8r, v9r = rule_precise(0.002, 2, 0.04, 2), rule_snr(30, 0.01, 1.5)
    for name, rule, first in (("v7", never, False), ("v8", v8r, False), ("v9", v9r, False), ("v10", v9r, True),
                              ("v11", lambda g: v8r(g) or v9r(g), True)):
        c = {coll: counts(coll, rule, first) for coll in COLLS}
        res["chosen"][name] = {f"{k}_{coll}": v for coll in COLLS for k, v in zip(("found", "added"), c[coll])}
    out = RESULTS / "dev" / "smallrules.json"
    json.dump(res, open_result(out, "w"), indent=1)
    print("wrote", out, json.dumps(res["chosen"]))


if __name__ == "__main__":
    main()
