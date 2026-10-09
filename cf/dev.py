"""Development harness: run a matrix-level procedure on stored trial matrices and score it
with the pre-registered definitions. Development uses the cryo collection only."""

from __future__ import annotations

import json
import sys
import time
from collections import defaultdict

import numpy as np

from cf.analyze import clopper_pearson, found, load_design
from cf.cmatrices import load


def run(design, procedure, select=None, **kw):
    trials = load_design(design)
    mats = load(design)
    rows = []
    for tid, t in trials.items():
        if select is not None and not select(t):
            continue
        m = mats.get(tid)
        if m is None:
            continue
        names = [str(x) for x in m["names"]]
        t0 = time.process_time()
        groups = procedure(m["rij"].astype(float), m["wij"].astype(float), **kw)
        cpu = time.process_time() - t0
        clusters = [[names[i] for i in g] for g in groups if len(g) >= 2]
        row = dict(id=tid, kind=t["kind"], w=t["w"], cpu=cpu, n_clusters=len(clusters),
                   clusters=clusters)
        if t["kind"] == "mixture":
            row.update(k=t["k"], pair="".join(sorted((t["majority_species"][0], t["minority_species"][0]))))
            if t["k"] >= 2:
                row["found"] = found(t, clusters)
        else:
            row["false_group"] = len(clusters) >= 2
        rows.append(row)
    return rows


def summary(rows):
    cell = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    fg = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["kind"] == "mixture" and r["k"] >= 2:
            c = cell[(r["pair"], r["w"])][r["k"]]
            c[0] += r["found"]
            c[1] += 1
        elif r["kind"] == "single":
            fg[r["w"]][0] += r["false_group"]
            fg[r["w"]][1] += 1
    lines = []
    for key in sorted(cell, key=lambda x: (x[0], -x[1])):
        lines.append(f"{key} " + " ".join(f"k{k}:{v[0]}/{v[1]}" for k, v in sorted(cell[key].items())))
    for w in sorted(fg, reverse=True):
        x, n = fg[w]
        lines.append(f"false groups w={w}: {x}/{n} {np.round(clopper_pearson(x, n), 3)}")
    return "\n".join(lines)
