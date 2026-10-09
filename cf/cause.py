"""Which stage of today's DIALS loses small groups: replay its clustering from each trial's
saved rij/wij (cf.replay, DIALS's own functions) with one stage changed at a time, and record
whether the minority is separable in the coordinates that OPTICS receives.

Variants:
  default         today's DIALS (replayed)
  floor2          OPTICS min_samples search starts at 2 instead of 5
  dims+1 / dims5  one more dimension than DIALS chooses / five dimensions
  floor2_dims+1   both
  heuristic       the paper's fixed heuristic min_samples, no search
  unit            OPTICS on coordinates scaled to unit length
Diagnostics in DIALS's coordinates: silhouette of the true minority/majority labelling
(cosine distance), and whether every minority crystal's nearest neighbour is a minority
crystal.
Writes results/cause/<design>/shard_<i>.jsonl. For the larger room-temperature trials the replay
can be restricted (environment): CF_CAUSE_VARIANTS (e.g. "default,floor2"), CF_CAUSE_WEDGES (e.g.
"60"), CF_CAUSE_KMAX (largest minority).
"""

from __future__ import annotations

from cf.result_io import open_result

import json
import os
import sys

import numpy as np
from sklearn.metrics import silhouette_score

from cf.analyze import found, load_design
from cf.cmatrices import load
from cf.paths import RESULTS
from cf.replay import coordinates, optics_labels
from cf.today import default_params


def clusters_from(labels, names):
    out = {}
    for n, l in zip(names, labels):
        if l >= 0:
            out.setdefault(int(l), []).append(n)
    return list(out.values())


def run(design, shard, nshards):
    trials = load_design(design)
    mats = load(design)
    only_w = os.environ.get("CF_CAUSE_WEDGES")
    kmax = int(os.environ.get("CF_CAUSE_KMAX", 10**6))
    keep_variants = os.environ.get("CF_CAUSE_VARIANTS")
    ids = sorted(t for t in trials if (trials[t]["kind"] == "single"
                 or (trials[t]["kind"] == "mixture" and 2 <= trials[t]["k"] <= kmax))
                 and (not only_w or f"{trials[t]['w']:g}" in only_w.split(",")))
    ids = ids[shard::nshards]
    od = RESULTS / "cause" / design
    od.mkdir(parents=True, exist_ok=True)
    out = od / f"shard_{shard:03d}.jsonl"
    done = set()
    if out.exists():
        done = {json.loads(l)["id"] for l in open_result(out)}
    params = default_params()
    for tid in ids:
        if tid in done or tid not in mats:
            continue
        t = trials[tid]
        m = mats[tid]
        names = [str(x) for x in m["names"]]
        r = m["rij"].astype(float)
        w = m["wij"].astype(float)
        rec = dict(id=tid)
        try:
            np.random.seed(params.seed)
            coords, dim, _ = coordinates(r, w, params)
            rec["dims"] = dim
            minority = t.get("minority", [])
            lab_true = np.array([1 if n in minority else 0 for n in names])
            if 0 < lab_true.sum() < len(names):
                rec["silhouette_true"] = float(silhouette_score(coords, lab_true, metric="cosine"))
                cn = coords / np.linalg.norm(coords, axis=1, keepdims=True)
                sim = cn @ cn.T
                np.fill_diagonal(sim, -np.inf)
                nn = sim.argmax(1)
                mi = np.flatnonzero(lab_true == 1)
                rec["minority_nn_pure"] = bool(np.all(lab_true[nn[mi]] == 1))
            variants = {
                "default": dict(coords=coords, dim=dim, kw={}),
                "floor2": dict(coords=coords, dim=dim, kw=dict(min_samples_floor=2)),
                "heuristic": dict(coords=coords, dim=dim, kw=dict(heuristic_only=True)),
                "unit": dict(coords=coords, dim=dim, kw=dict(unit_vectors=True)),
            }
            for label, d_force in (() if keep_variants else (("dims+1", dim + 1), ("dims5", 5))):
                np.random.seed(params.seed)
                c2, d2, _ = coordinates(r, w, params, dims=d_force)
                variants[label] = dict(coords=c2, dim=d2, kw={})
                if label == "dims+1":
                    variants["floor2_dims+1"] = dict(coords=c2, dim=d2, kw=dict(min_samples_floor=2))
            if keep_variants:
                variants = {k: v for k, v in variants.items() if k in keep_variants.split(",")}
            for label, v in variants.items():
                labels, ms = optics_labels(v["coords"], v["dim"], params, **v["kw"])
                cl = clusters_from(labels, names)
                rec[label] = dict(n_clusters=len(cl), min_samples=int(ms))
                if t["kind"] == "mixture":
                    rec[label]["found"] = found(t, cl)
                else:
                    rec[label]["false_group"] = len(cl) >= 2
        except Exception as e:
            rec["error"] = repr(e)
        with open_result(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    run(sys.argv[1], int(os.environ.get("CF_SHARD", 0)), int(os.environ.get("CF_NSHARDS", 1)))
