"""Score trial results with the pre-registered definitions (PREREGISTRATION.md) and
summarise them per cell: detection probability with exact intervals, detection limits,
false-group rates."""

from __future__ import annotations

from cf.result_io import open_result, glob_results

import glob
import json
import math
from collections import defaultdict

from scipy.stats import beta

from cf.paths import RESULTS, base_design


def clopper_pearson(x: int, n: int, level: float = 0.95):
    if n == 0:
        return (math.nan, math.nan)
    a = (1 - level) / 2
    lo = 0.0 if x == 0 else beta.ppf(a, x, n - x + 1)
    hi = 1.0 if x == n else beta.ppf(1 - a, x + 1, n - x)
    return (float(lo), float(hi))


def found(trial: dict, clusters: list[list[str]]) -> bool:
    """The one pre-registered definition (k >= 2)."""
    m = set(trial["minority"])
    k = len(m)
    need = max(2, k / 2)
    for c in clusters:
        hit = len(m.intersection(c))
        if hit >= need and hit > len(c) / 2:
            return True
    return False


def isolated(trial: dict, clusters: list[list[str]]) -> bool:
    """k = 1: the crystal is in no reported cluster that contains a majority crystal."""
    (x,) = trial["minority"]
    maj = set(trial["majority"])
    return not any(x in c and maj.intersection(c) for c in clusters)


def load_design(name: str) -> dict[str, dict]:
    return {
        t["id"]: t
        for t in (json.loads(l) for l in open_result(RESULTS / "design" / f"{base_design(name)}.jsonl"))
    }


def load_results(design: str, procedure: str) -> dict[str, dict]:
    if procedure == "today":  # today's DIALS (OPTICS on DIALS's own matrix) for every representation
        design = base_design(design)
    out = {}
    for f in sorted(glob_results(str(RESULTS / "trials" / design / procedure / "shard_*.jsonl"))):
        for line in open_result(f):
            r = json.loads(line)
            out[r["id"]] = r
    return out


def score(design: str, procedure: str):
    """Per-trial outcomes. Failed runs count as reporting nothing (not found, no group)."""
    trials = load_design(design)
    res = load_results(design, procedure)
    rows = []
    for tid, t in trials.items():
        r = res.get(tid)
        if r is None:
            continue
        clusters = r.get("clusters", []) if "error" not in r else []
        row = dict(id=tid, kind=t["kind"], coll=t["coll"], w=t["w"], error="error" in r,
                   n_clusters=len(clusters), n_outliers=len(r.get("outliers", [])),
                   cpu_s=r.get("cpu_s"), dims=r.get("dims"))
        if t["kind"] == "mixture":
            row.update(pair="".join(sorted((t["majority_species"][0], t["minority_species"][0]))),
                       direction=t["majority_species"][0] + t["minority_species"][0],
                       k=t["k"], draw=t["draw"], m=len(t["majority"]))
            if t["k"] >= 2:
                row["found"] = found(t, clusters)
            else:
                row["isolated"] = isolated(t, clusters)
        else:
            row.update(species=t["species"], part=t["part"],
                       false_group=len(clusters) >= 2, none=len(clusters) == 0)
        rows.append(row)
    return rows


def detection_table(rows):
    """{(coll, pair, w): {k: (found, n)}} pooled over both directions."""
    cells = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in rows:
        if r["kind"] == "mixture" and r["k"] >= 2:
            c = cells[(r["coll"], r["pair"], r["w"])][r["k"]]
            c[0] += int(r["found"])
            c[1] += 1
    return cells


def limit(per_k: dict, p: float = 0.9):
    """Smallest k >= 2 with detection >= p there and at every larger tested k."""
    ks = sorted(per_k)
    ok = [per_k[k][1] > 0 and per_k[k][0] / per_k[k][1] >= p for k in ks]
    for i, k in enumerate(ks):
        if all(ok[i:]):
            return k
    return None


def false_group_table(rows):
    cells = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["kind"] == "single":
            c = cells[(r["coll"], r["w"])]
            c[0] += int(r["false_group"])
            c[1] += 1
    return cells


def score_clusters(design: str, results: dict[str, list]) -> list[dict]:
    """Score {trial id: reported clusters (lists of crystal names)} like score()."""
    trials = load_design(design)
    rows = []
    for tid, clusters in results.items():
        t = trials[tid]
        row = dict(id=tid, kind=t["kind"], coll=t["coll"], w=t["w"], error=False,
                   n_clusters=len(clusters))
        if t["kind"] == "mixture":
            row.update(pair="".join(sorted((t["majority_species"][0], t["minority_species"][0]))),
                       direction=t["majority_species"][0] + t["minority_species"][0],
                       k=t["k"], draw=t["draw"], m=len(t["majority"]))
            if t["k"] >= 2:
                row["found"] = found(t, clusters)
            else:
                row["isolated"] = isolated(t, clusters)
        else:
            row.update(species=t["species"], part=t["part"],
                       false_group=len(clusters) >= 2, none=len(clusters) == 0)
        rows.append(row)
    return rows


def load_dev(design: str, tag: str) -> dict[str, list]:
    out = {}
    for f in sorted(glob_results(str(RESULTS / "dev" / design / tag / "shard_*.jsonl"))):
        for line in open_result(f):
            r = json.loads(line)
            if "clusters" in r:  # (development checks without clusters are not procedures)
                out[r["id"]] = r["clusters"]
    return out
