"""Development (amendment 14): Union of v6's detector groups and an alternative run's (combine_groups semantics), carved
into today's OPTICS clusters, scored against v6 alone."""

from cf.result_io import open_result, glob_results

import json, sys, glob
from collections import defaultdict
from cf.analyze import load_results, score_clusters
from cf.newproc import carve
import numpy as np

def reported(groups, n):
    groups = sorted(groups, key=len, reverse=True)
    if len(groups) <= 1:
        return []
    return groups[1:] if 2 * len(groups[0]) > n else groups

def combine(first, second, names):
    n = len(names)
    out = [list(g) for g in reported(first, n)]
    taken = {x for g in out for x in g}
    for g in reported(second, n):
        if not set(g) & taken:
            out.append(list(g)); taken |= set(g)
    rest = [x for x in names if x not in taken]
    return ([rest] if rest else []) + out

alt = sys.argv[1]
today = load_results("cryo", "today")
v6 = load_results("cryo", "v6")
recs = {}
for f in glob_results(f"results/dev/cryo/alt_{alt}/shard_*.jsonl"):
    for l in open_result(f):
        r = json.loads(l); recs[r["id"]] = r
design = {json.loads(l)["id"]: json.loads(l) for l in open_result("results/design/cryo.jsonl")}
uni, v6c = {}, {}
for tid, r in recs.items():
    a = v6.get(tid)
    if a is None or "error" in a:
        continue
    names = sorted({x for g in a["detector_groups"] for x in g} | {x for g in r["detector_groups"] for x in g})
    t = design[tid]
    all_names = t["crystals"] if "crystals" in t else t["majority"] + t["minority"]
    groups = combine(a["detector_groups"], r["detector_groups"], names)
    uni[tid] = carve(today[tid].get("clusters", []), groups, all_names, "half")
    v6c[tid] = a["clusters"]
for tag, res in (("v6", v6c), (f"v6+{alt}", uni)):
    rows = score_clusters("cryo", res)
    fd = defaultdict(lambda: [0, 0]); fg = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["kind"] == "single":
            fg[r["w"]][0] += int(r["false_group"]); fg[r["w"]][1] += 1
        elif r.get("k", 0) >= 2:
            fd[(r["pair"], r["w"])][0] += int(r["found"]); fd[(r["pair"], r["w"])][1] += 1
    print(f"{tag:14s} found {sum(x[0] for x in fd.values())}/{sum(x[1] for x in fd.values())} FG {sum(x[0] for x in fg.values())}/{sum(x[1] for x in fg.values())}", dict((k, x[0]) for k, x in sorted(fg.items()) if x[0]))
    print("    ", " ".join(f"{p}{w:g}:{x[0]}" for (p, w), x in sorted(fd.items())))
