"""Every number of FINDINGS.md / COMPLETE.md, computed from the results files.

Writes results/report/numbers.json. Each entry says which procedure, which trials and which
reference (Zenodo species labels, the paper, the limit analyses) it is computed against.
Run: python -m cf.report (then python -m cf.render).
"""

from __future__ import annotations

from cf.result_io import open_result, glob_results, result_exists

import glob
import json
import os
import re
import subprocess
from collections import Counter, defaultdict

import numpy as np

from cf.analyze import (
    clopper_pearson,
    detection_table,
    limit,
    load_design,
    load_dev,
    load_results,
    score,
    score_clusters,
)
from cf.paths import RESULTS, ROOT, manifest

OUT = RESULTS / "report"
COLLS = ("cryo", "rt")
PAIRS = {"bh": "bovine-human", "bp": "bovine-porcine", "hp": "human-porcine"}
# the frozen procedures (PREREGISTRATION.md amendments 3-20)
PROCS = ("today", "v11", "v10", "v9", "v8", "v7", "v6", "v6_1", "v5", "v4", "v3", "v2", "v1", "v4_ratio_only", "v4_coherence_only", "v4_bootstrap",
         "v3_ratio_only", "v3_coherence_only", "v3_analytic", "v6_m3", "v6_m4", "v6_m6", "v6_m8")
PRIMARY = "v11"  # the change as proposed (0001-0015, its defaults); amendment 26
TAKEOUT = "v10"  # v9 with the take-out decided on all groups first (0001-0014); amendment 23
SNR = "v9"  # v7 and small groups far out of their noise, merged before the take-out (0001-0013); amendment 22
FLOOR_ONLY = "v7"  # v6 with OPTICS's floor for every group (0001-0010); amendment 17
PRECISE = "v8"  # v7 and small groups precisely measured (0001-0011); amendment 20
PREREGISTERED = "v6"  # the procedure frozen before room temperature (amendment 12): the held-out judgement
JUDGED = (PREREGISTERED, FLOOR_ONLY, PRECISE, SNR, TAKEOUT, PRIMARY)  # in the order they were frozen
# each judged configuration: amendment, patches, what it does
CONFIGURATIONS = [
    (PREREGISTERED, 12, "0001-0007", "after OPTICS, the calibrated detector's groups down to two data sets taken "
     "out as clusters (frozen before any room-temperature result)"),
    (FLOOR_ONLY, 17, "0001-0010", f"{PREREGISTERED}, with groups smaller than five returned to the background population"),
    (PRECISE, 20, "0001-0011", f"{FLOOR_ONLY}, and a group smaller than five when its coherence is precisely measured "
     "(standard error at most 0.002), two standard errors above 0.04 x 2/k and twice the rest's structure"),
    (SNR, 22, "0001-0013", f"{FLOOR_ONLY}, and a group smaller than five when its coherence stands 30 standard errors "
     "out of the counting noise, two above 0.01 x 2/k and 1.5 times above the rest's structure"),
    (TAKEOUT, 23, "0001-0014", f"{SNR}, with the groups to take out of the OPTICS clusters decided on all of the "
     "detector's groups before the minimum size"),
    (PRIMARY, 26, "0001-0015", f"{TAKEOUT}, and a group smaller than five admitted by {PRECISE}'s condition as well "
     f"as {SNR}'s (small_groups and small_groups.precise)"),
]
ABLATIONS = ("floor2", "heuristic", "unit", "dims+1", "dims5", "floor2_dims+1")


def commit():
    provenance = json.loads((ROOT / "PROVENANCE.json").read_text())
    return provenance["verdict_commit"]


def ci(x, n):
    if n == 0:
        return None
    lo, hi = clopper_pearson(x, n)
    return [x, n, round(lo, 4), round(hi, 4)]


# ---------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------


def data_summary():
    man = manifest()
    out = {}
    for coll in COLLS:
        rows = [r for r in man if r["collection"] == coll]
        procs = [json.load(open_result(f)) for f in glob_results(str(RESULTS / "processing" / coll / "*.json"))]
        out[coll] = dict(
            crystals=len(rows),
            by_species=dict(Counter(r["species"] for r in rows)),
            processed=len(procs),
            index_retry=sum(bool(p.get("index_retry")) for p in procs),
            images_median=float(np.median([p["n_images"] for p in procs])) if procs else None,
            oscillation=procs[0]["oscillation"][1] if procs else None,
            dials_version=procs[0]["dials_version"] if procs else None,
        )
    for d in COLLS:
        tr = load_design(d)
        out[d]["trials"] = len(tr)
        out[d]["mixtures"] = sum(t["kind"] == "mixture" for t in tr.values())
        out[d]["mixtures_k2"] = sum(t["kind"] == "mixture" and t["k"] >= 2 for t in tr.values())
        out[d]["singles"] = sum(t["kind"] == "single" for t in tr.values())
        out[d]["wedges"] = sorted({t["w"] for t in tr.values()}, reverse=True)
    return out


# ---------------------------------------------------------------------------------------
# Fidelity
# ---------------------------------------------------------------------------------------


def multiplex_result(tag):
    d = RESULTS / "fidelity" / tag
    if not (d / "optics_summary.txt").exists():
        return None
    order = [l.strip() for l in open_result(d / "input_order.txt") if l.strip()]
    sp = {r["crystal"]: r["species"] for r in manifest()}
    txt = open_result(d / "optics_summary.txt").read()
    head = re.search(r"OPTICS identified (\d+) clusters and (\d+) outlier", txt)
    if not head:
        return dict(error="no OPTICS summary")
    clusters = []
    for m in re.finditer(r"Cluster (\d+)\n.*?Datasets:([\d,]+)", txt, re.S):
        ids = [int(x) for x in m.group(2).split(",")]
        clusters.append(dict(Counter(sp[order[i]] for i in ids)))
    logs = glob_results(str(d / "*dials.correlation_matrix.log"))
    log = open_result(logs[0]).read() if logs else ""
    dims = re.search(r"Using (\d+) dimensions", log)
    ms = re.search(r"Set Minimum Samples to (\d+)", log)
    res_filter = re.search(r"High resolution limit set to: ([\d.]+)", log)  # clustering's filter
    return dict(n_clusters=int(head.group(1)), n_outliers=int(head.group(2)), clusters=clusters,
                dims=int(dims.group(1)) if dims else None,
                min_samples=int(ms.group(1)) if ms else None,
                resolution_filter=float(res_filter.group(1)) if res_filter else None,
                pure=all(len(c) == 1 for c in clusters))


def purity_table(design):
    trials = load_design(design)
    res = load_results(design, "today")
    sp = {r["crystal"]: r["species"] for r in manifest()}
    out = []
    for tid, t in trials.items():
        r = res.get(tid)
        if r is None:
            continue
        cl = r.get("clusters", []) if "error" not in r else []
        comp = [dict(Counter(sp.get(x, "unclassified") for x in c)) for c in cl]
        out.append(dict(id=tid, part=t.get("part"), w=t["w"], n=len(t["crystals"]), clusters=comp,
                        n_outliers=len(r.get("outliers", [])), dims=r.get("dims"),
                        pure=all(len([k for k in c if k != "unclassified"]) <= 1 for c in comp),
                        error=r.get("error") or r.get("dials_error")))
    return sorted(out, key=lambda x: (x["part"] or "", -x["w"]))


def dials_test_coords():
    f = RESULTS / "fidelity" / "dials_test_coords.json"
    return json.load(open_result(f)) if f.exists() else None


# ---------------------------------------------------------------------------------------
# Detection, false groups, misses
# ---------------------------------------------------------------------------------------


def detection_summary(rows):
    det = detection_table(rows)
    cells = {}
    for (coll, pair, w), per_k in det.items():
        cells[f"{coll}|{pair}|{w:g}"] = dict(
            per_k={str(k): list(per_k[k]) for k in sorted(per_k)},
            k90=limit(per_k, 0.9), k50=limit(per_k, 0.5))
    return cells


def false_groups(rows):
    out = defaultdict(lambda: dict(first=[0, 0], all=[0, 0]))
    for r in rows:
        if r["kind"] != "single":
            continue
        key = f"{r['coll']}|{r['w']:g}"
        out[key]["all"][0] += int(r["false_group"])
        out[key]["all"][1] += 1
        if "-s" not in r["id"].split("-w")[-1]:
            out[key]["first"][0] += int(r["false_group"])
            out[key]["first"][1] += 1
    return dict(out)


def isolated(rows):
    out = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["kind"] == "mixture" and r.get("k") == 1:
            key = f"{r['coll']}|{r['pair']}|{r['w']:g}"
            out[key][0] += int(r["isolated"])
            out[key][1] += 1
    return dict(out)


def singles_complete(design, rows):
    """Whether rows hold every single-population trial of the design: a false-group rate, and so
    being "no more false groups than today's", is established only on all of them."""
    need = {t for t, v in load_design(design).items() if v["kind"] == "single"}
    return need <= {r["id"] for r in rows}


def fg_safe_found(design):
    """Minorities found by a frozen procedure that reports no more false groups than today's
    DIALS on this collection (evidence that they are distinguishable; amendment 6)."""
    procs = [p for p in PROCS if (RESULTS / "trials" / design / p).exists()]
    rows = {p: score(design, p) for p in procs}
    procs = [p for p in procs if singles_complete(design, rows[p])]
    fg = {p: sum(r["false_group"] for r in rows[p] if r["kind"] == "single") for p in procs}
    out = set()
    for p in procs:
        if fg[p] <= fg.get("today", 0):
            out |= {r["id"] for r in rows[p] if r["kind"] == "mixture" and r.get("k", 0) >= 2
                    and r["found"]}
    return out


def fg_safe_found_both(design):
    """Minorities found by a procedure that reports no more false groups than today's DIALS on
    either collection (one configuration for every sample; amendment 18)."""
    procs = [p for p in PROCS if all((RESULTS / "trials" / c / p).exists() for c in COLLS)]
    ok = []
    for p in procs:
        safe = True
        for c in COLLS:
            rows = score(c, p)
            safe &= singles_complete(c, rows)
            fg = sum(r["false_group"] for r in rows if r["kind"] == "single")
            fg0 = sum(r["false_group"] for r in score(c, "today") if r["kind"] == "single")
            safe &= fg <= fg0
        if safe:
            ok.append(p)
    out = set()
    for p in ok:
        out |= {r["id"] for r in score(design, p) if r["kind"] == "mixture" and r.get("k", 0) >= 2 and r["found"]}
    return out


def refined_limits(design):
    """The refined limits (amendments 4-16; main), the strict version of amendment 4, and the
    limits across the collections (amendment 18)."""
    from cf.limits3 import classify

    try:
        found = fg_safe_found(design)
        main = classify(design, found_ids=found)
        strict = classify(design, tol=0.0, found_ids=found, alpha=0.01)
    except Exception:
        return {}, {}, {}
    try:
        cross = classify(design, found_ids=fg_safe_found_both(design), cross=True)
    except Exception:
        cross = {}
    return main, strict, cross


def amendment24_effect(design):
    """Minorities whose classification amendment 24 (comparators of every size) moves from separable
    to a limit of the data, per collection classification and per judged procedure's misses."""
    from cf.limits3 import classify

    try:
        found = fg_safe_found(design)
        found_both = fg_safe_found_both(design)
        views = dict(per_collection=(classify(design, found_ids=found),
                                     classify(design, found_ids=found, design_sizes_only=True)),
                     cross=(classify(design, found_ids=found_both, cross=True),
                            classify(design, found_ids=found_both, cross=True, design_sizes_only=True)))
    except Exception:
        return None
    out = {}
    for name, (new, old) in views.items():
        moved = {t for t in new if t in old and new[t]["data_limit"] and not old[t]["data_limit"]}
        back = {t for t in new if t in old and old[t]["data_limit"] and not new[t]["data_limit"]}
        per = {}
        for p in JUDGED:
            if (RESULTS / "trials" / design / p).exists():
                misses = {r["id"] for r in score(design, p) if r["kind"] == "mixture" and r.get("k", 0) >= 2
                          and not r["found"]}
                per[p] = len(misses & moved)
        out[name] = dict(moved=len(moved), moved_back=len(back), misses=per)
    return out


def preregistered_limits(design):
    f = RESULTS / "limits" / f"{design}.jsonl"
    return {json.loads(l)["id"]: json.loads(l) for l in open_result(f)} if result_exists(f) else {}


CATS = ("found", "miss: counting noise", "miss: structure inside one population",
        "miss: membership not resolved", "miss: noise known only from repeat measurements",
        "miss: procedure", "miss: not classified")


def classify_misses(rows, refined, pre):
    """Per cell: found / data limit (counting noise; structure inside one population) /
    procedure, by the refined rule (amendment 4); and by the pre-registered rule."""
    out = defaultdict(Counter)
    out_pre = defaultdict(Counter)
    for r in rows:
        if r["kind"] != "mixture" or r.get("k", 0) < 2:
            continue
        key = f"{r['coll']}|{r['pair']}|{r['w']:g}"
        q = refined.get(r["id"])
        if r["found"]:
            cat = "found"
        elif q is None:
            cat = "miss: not classified"
        elif q["within_counting_noise"]:
            cat = "miss: counting noise"
        elif q["within_population_structure"]:
            cat = "miss: structure inside one population"
        elif q.get("membership_unresolved"):
            cat = "miss: membership not resolved"
        elif q.get("within_single_experiment_noise"):
            cat = "miss: noise known only from repeat measurements"
        else:
            cat = "miss: procedure"
        out[key][cat] += 1
        p = pre.get(r["id"])
        if r["found"]:
            out_pre[key]["found"] += 1
        elif p is None or "within_noise" not in p:
            out_pre[key]["miss: not classified"] += 1
        elif p["within_noise"]:
            out_pre[key]["miss: data (pre-registered)"] += 1
        else:
            out_pre[key]["miss: procedure (pre-registered)"] += 1
    return {k: dict(v) for k, v in out.items()}, {k: dict(v) for k, v in out_pre.items()}


def data_k90(design, refined, pre, strict=None):
    """Per cell, the k90 the data allow: minorities that are not data limits (refined rule of
    amendment 6; and the strict rule of amendment 4), and by the pre-registered rule (oracle
    and same-trial decoy; sizes at which the oracle cannot reach p < 0.05, fewer than 19
    subsets of the majority to compare with, are left out)."""
    import math

    trials = load_design(design)
    det = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    det_strict = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    det_single = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    det_pre = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for tid, t in trials.items():
        if t["kind"] != "mixture" or t["k"] < 2:
            continue
        pair = "".join(sorted((t["majority_species"][0], t["minority_species"][0])))
        key = f"{design}|{pair}|{t['w']:g}"
        q = refined.get(tid)
        c = det[key][t["k"]]
        c[0] += int(q is not None and not q["data_limit"])
        c[1] += 1
        c = det_single[key][t["k"]]
        c[0] += int(q is not None and not q["data_limit"] and not q.get("within_single_experiment_noise"))
        c[1] += 1
        if strict is not None:
            q = strict.get(tid)
            c = det_strict[key][t["k"]]
            c[0] += int(q is not None and not q["data_limit"])
            c[1] += 1
        if math.comb(len(t["majority"]), t["k"]) >= 19:
            p = pre.get(tid)
            c = det_pre[key][t["k"]]
            c[0] += int(p is not None and "within_noise" in p and not p["within_noise"])
            c[1] += 1
    out = {}
    for key in det:
        out[key] = dict(per_k={str(k): v for k, v in sorted(det[key].items())},
                        k90=limit(det[key], 0.9), k50=limit(det[key], 0.5),
                        k90_single_experiment=limit(det_single[key], 0.9),
                        k90_strict=limit(det_strict[key], 0.9) if strict is not None else None,
                        k90_preregistered=limit(det_pre[key], 0.9),
                        per_k_preregistered={str(k): v for k, v in sorted(det_pre[key].items())})
    return out


# ---------------------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------------------


def cost(design, procedure):
    res = load_results(design, procedure)
    vals, parts = [], defaultdict(list)
    for r in res.values():
        if "error" in r:
            continue
        if procedure == "today":
            c = r.get("cpu_s")
        else:
            c = sum((r.get(x) or 0) for x in ("optics_cpu_s", "noise_cpu_s", "detector_cpu_s"))
            for x in ("optics_cpu_s", "noise_cpu_s", "detector_cpu_s"):
                if r.get(x) is not None:
                    parts[x].append(r[x])
        if c is not None:
            vals.append(c)
    vals = np.array([v for v in vals if np.isfinite(v)])
    return dict(median=float(np.median(vals)) if len(vals) else None, n=int(len(vals)),
                parts={k: float(np.median(v)) for k, v in parts.items()})


# ---------------------------------------------------------------------------------------
# Cause (today's DIALS replayed with one stage changed)
# ---------------------------------------------------------------------------------------


def _cause_rows(design):
    return [json.loads(l) for f in glob_results(str(RESULTS / "cause" / design / "shard_*.jsonl"))
            for l in open_result(f)]


def cause_table(design):
    trials = load_design(design)
    rows = _cause_rows(design)
    if not rows:
        return None
    variants = ("default",) + ABLATIONS
    out = {}
    for v in variants:
        det = defaultdict(lambda: [0, 0])
        fg = defaultdict(lambda: [0, 0])
        for r in rows:
            if "error" in r or v not in r:
                continue
            t = trials[r["id"]]
            if t["kind"] == "mixture":
                d = det[f"{t['w']:g}|{t['k']}"]
                d[0] += int(r[v]["found"])
                d[1] += 1
            else:
                d = fg[f"{t['w']:g}"]
                d[0] += int(r[v]["false_group"])
                d[1] += 1
        out[v] = dict(detection=dict(det), false_groups=dict(fg))
    sep = defaultdict(list)
    for r in rows:
        t = trials[r["id"]]
        if "error" in r or t["kind"] != "mixture" or "silhouette_true" not in r:
            continue
        sep[f"{t['w']:g}|{t['k']}"].append((r["silhouette_true"], r["minority_nn_pure"], r["dims"]))
    out["separability"] = {k: dict(silhouette_median=float(np.median([a[0] for a in v])),
                                   nn_pure=float(np.mean([a[1] for a in v])),
                                   dims_median=float(np.median([a[2] for a in v])), n=len(v))
                           for k, v in sep.items()}
    out["errors"] = sum("error" in r for r in rows)
    out["n"] = len(rows)
    return out


def _cause_as_rows(design, variant):
    trials = load_design(design)
    res = {}
    for r in _cause_rows(design):
        if "error" in r or variant not in r:
            continue
        t = trials[r["id"]]
        row = dict(id=r["id"], kind=t["kind"], coll=t["coll"], w=t["w"])
        if t["kind"] == "mixture":
            row.update(pair="".join(sorted((t["majority_species"][0], t["minority_species"][0]))),
                       k=t["k"], found=bool(r[variant]["found"]))
        else:
            row.update(false_group=bool(r[variant]["false_group"]))
        res[r["id"]] = row
    return list(res.values())


# ---------------------------------------------------------------------------------------
# Frontier
# ---------------------------------------------------------------------------------------


def frontier(design):
    """Every procedure run on every trial of the design, summarised by the fraction of
    minorities (k >= 2) found, the false-group rate and the median CPU seconds per trial; the
    non-dominated set. Development prototypes (cryo) are listed with their detector cost only."""
    procs = {}
    for proc in PROCS:
        if (RESULTS / "trials" / design / proc).exists():
            procs[proc] = (score(design, proc), cost(design, proc)["median"], "frozen")
    for v in ABLATIONS:
        rows = _cause_as_rows(design, v)
        if rows:
            procs[f"today, {v} (replay)"] = (rows, None, "ablation")
    if design == "cryo":
        for d in sorted(glob_results(str(RESULTS / "dev" / design / "*"))):
            tag = os.path.basename(d)
            if os.path.isdir(d) and load_dev(design, tag):
                rows = score_clusters(design, load_dev(design, tag))
                procs[f"prototype {tag} (detector alone)"] = (rows, None, "development")
    trials = load_design(design)
    n_mix = sum(1 for t in trials.values() if t["kind"] == "mixture" and t["k"] >= 2)
    n_single = sum(1 for t in trials.values() if t["kind"] == "single")
    out = {}
    for name, (rows, c, kind) in procs.items():
        mix = [r for r in rows if r["kind"] == "mixture" and r.get("k", 0) >= 2]
        sing = [r for r in rows if r["kind"] == "single"]
        if len(mix) < 0.95 * n_mix or len(sing) < 0.95 * n_single:
            out[name] = dict(incomplete=True, n_mix=len(mix), n_single=len(sing), kind=kind)
            continue
        out[name] = dict(found=sum(r["found"] for r in mix), n_mix=len(mix),
                         false_groups=sum(r["false_group"] for r in sing), n_single=len(sing),
                         cpu_s_median=c, kind=kind)
    complete = {k: v for k, v in out.items() if not v.get("incomplete")}
    for k, v in complete.items():
        dominated = False
        for k2, v2 in complete.items():
            if k2 == k:
                continue
            f1, f2 = v["found"] / v["n_mix"], v2["found"] / v2["n_mix"]
            g1, g2 = v["false_groups"] / v["n_single"], v2["false_groups"] / v2["n_single"]
            c1 = v["cpu_s_median"] if v["cpu_s_median"] is not None else np.inf
            c2 = v2["cpu_s_median"] if v2["cpu_s_median"] is not None else np.inf
            if f2 >= f1 and g2 <= g1 and c2 <= c1 and (f2 > f1 or g2 < g1 or c2 < c1):
                dominated = True
        v["non_dominated"] = not dominated
    return out


def unexplained(design, refined, proc=PRIMARY):
    """The primary procedure's misses that no measurement puts in the data, with what is
    known about each."""
    trials = load_design(design)
    res = load_results(design, proc)
    found_by = defaultdict(list)
    for p in PROCS:
        if p != proc and (RESULTS / "trials" / design / p).exists():
            for r in score(design, p):
                if r["kind"] == "mixture" and r.get("k", 0) >= 2 and r["found"]:
                    found_by[r["id"]].append(p)
    # v6's groups of two to four as the procedure measures them (cf/devsmall.py): for a missed minority
    # that v6 reports as such a group, how far its coherence stands out of its noise
    small = {}
    for f in glob_results(str(RESULTS / "dev" / design / "smallgroups" / "shard_*.jsonl")):
        for line in open_result(f):
            rec = json.loads(line)
            small[rec["id"]] = rec["groups"]
    out = []
    for r in score(design, proc):
        if r["kind"] != "mixture" or r.get("k", 0) < 2 or r["found"]:
            continue
        q = refined.get(r["id"])
        if q is None or q["data_limit"] or q.get("within_single_experiment_noise"):
            continue
        t = trials[r["id"]]
        unit = q.get("detectable_unit")
        u = q if unit in (None, "minority") else next(x for x in q.get("subgroups", []) if x["members"] == unit)
        pick = lambda *ks: min([v for v in (_finite(u.get(k)) for k in ks) if v is not None], default=None)
        out.append(dict(id=r["id"], w=t["w"], k=t["k"], pair=PAIRS[r["pair"]],
                        majority=t["majority_species"], minority=t["minority_species"],
                        unit_size=len(u.get("members", [])) or t["k"],
                        unit_minority=u.get("n_minority", len(u.get("members", [])) or t["k"]),
                        p_noise=pick("p_noise_any_size", "p_noise_any_size_iv"),
                        p_model=pick("p_noise_model", "p_noise_model_iv"),
                        detector_idle=len(res.get(r["id"], {}).get("detector_log", [])) == 0,
                        found_by=found_by.get(r["id"], []),
                        v6_group=next(({"k": g["k"], "S": g["S"], "se": g["se"],
                                        "snr": g["S"] / g["se"] if g["se"] > 0 else None}
                                       for g in sorted(small.get(r["id"], []), key=lambda g: -g["k"])
                                       if set(g["members"]) <= set(t["minority"])), None)))
    return out


def dropped_by_cosym(design):
    """Trials where cosym dropped a crystal after DIALS's own pre-filter (the matrix has fewer
    rows than the crystals DIALS used), and whether today's DIALS still names the dropped
    crystal in a reported cluster (which can only happen if it maps cluster members to the
    wrong crystals)."""
    from cf.cmatrices import load

    mats = load(design)
    if not mats:
        # The report needs matrix row names, not the regenerable correlation arrays.
        with open_result(RESULTS / "matrix_names" / f"{design}.jsonl") as fh:
            mats = {r["id"]: r for r in (json.loads(line) for line in fh)}
    res = load_results(design, "today")
    out = []
    for tid, m in mats.items():
        r = res.get(tid, {})
        used = r.get("used", [])
        names = [str(x) for x in m["names"]]
        if len(used) == len(names):
            continue
        dropped = sorted(set(used) - set(names))
        clusters = r.get("clusters", [])
        reported = {x for c in clusters for x in c} | set(r.get("outliers", []))
        out.append(dict(id=tid, dropped=dropped,
                        position=[used.index(x) for x in dropped], n=len(used),
                        dropped_last=all(used.index(x) == len(used) - 1 for x in dropped),
                        names_shifted=reported != set(names),
                        dropped_named_in_cluster=any(x in c for c in clusters for x in dropped)))
    return dict(n_trials=len(mats), affected=len(out),
                dropped_last=sum(o["dropped_last"] for o in out),
                names_shifted=sum(o["names_shifted"] for o in out),
                dropped_named=sum(o["dropped_named_in_cluster"] for o in out), trials=out)


def paired_false_groups(design, a="today", b=PRIMARY):
    """Single-population trials with a false group under a only, b only, both."""
    ra = {r["id"]: r for r in score(design, a) if r["kind"] == "single"}
    rb = {r["id"]: r for r in score(design, b) if r["kind"] == "single"}
    c = Counter()
    for tid in set(ra) & set(rb):
        c[(ra[tid]["false_group"], rb[tid]["false_group"])] += 1
    return dict(n=sum(c.values()), only_a=c[(True, False)], only_b=c[(False, True)],
                both=c[(True, True)], neither=c[(False, False)], a=a, b=b)


# ---------------------------------------------------------------------------------------
# Equivalence of the end-to-end patched DIALS with the from-arrays runs; tests; noise
# ---------------------------------------------------------------------------------------


def equivalence():
    """The patched DIALS run end to end on a sample of trials (scaling, cosym, clustering and the
    change), against the same module run on the trials' saved arrays: same detector groups, same
    final clusters, same outcome (found / false group). Also the end-to-end run of the frozen v6
    patch (0001-0007, `groups_v6_0007`): how many trials DIALS failed on."""
    from cf.analyze import found as found_rule

    design = load_design("e2e_sample")
    out = {}
    for tag in ("groups_v11", "groups_v10", "groups_v9", "groups_v8", "groups", "groups_v6_0007"):
        e2e = load_results("e2e_sample", tag)
        if not e2e:
            continue
        # failures of the patched DIALS where today's DIALS did not fail on the same trial
        today = {c: load_results(c, "today") for c in {t["coll"] for t in design.values()}}
        new_err = {tid: r["dials_error"] for tid, r in e2e.items() if r.get("dials_error")
                   and not today[design[tid]["coll"]].get(tid, {}).get("dials_error")}
        # the failing trials run alone (results/trials/debug_e2e/<configuration>), where that was done
        alone_tag = {"groups": "v7", "groups_v11": "v11"}.get(tag)
        alone = load_results("debug_e2e", alone_tag) if alone_tag and (RESULTS / "trials" / "debug_e2e" / alone_tag).exists() else {}
        rec_alone = {tid: (alone[tid].get("dials_error") or "no failure") for tid in new_err if tid in alone}
        rec = dict(n_run=len(e2e), dials_errors=len(new_err), alone=rec_alone,
                   dials_errors_as_today=sum(1 for r in e2e.values() if r.get("dials_error")) - len(new_err),
                   error_types=sorted({str(e)[:80] for e in new_err.values()}),
                   n=0, same_clusters=0, same_detector_groups=0, same_outcome=0, differ=[])
        # the runs of 0001-0014, 0001-0013, 0001-0011 and 0001-0010 make v10's, v9's, v8's and v7's clusters; the detector's own groups
        # (before the minimum group size) are v6's; the run of the frozen 0001-0007 is v6 in both
        cl_proc = {"groups_v11": "v11", "groups_v10": "v10", "groups_v9": "v9", "groups_v8": "v8", "groups": "v7"}.get(tag, PREREGISTERED)
        for tid, r in e2e.items():
            t = design[tid]
            fa = load_results(t["coll"], cl_proc).get(tid)
            fd = load_results(t["coll"], PREREGISTERED).get(tid)
            if fa is None or fd is None or "error" in r or "error" in fa:
                continue
            rec["n"] += 1
            norm = lambda cl: sorted(sorted(c) for c in cl)
            same = norm(r.get("clusters", [])) == norm(fa.get("clusters", []))
            rec["same_clusters"] += int(same)
            rec["same_detector_groups"] += int(norm(r.get("detector_groups", [])) == norm(fd.get("detector_groups", [])))
            if t["kind"] == "mixture" and t["k"] >= 2:
                outcome = found_rule(t, r.get("clusters", [])) == found_rule(t, fa.get("clusters", []))
            elif t["kind"] == "single":
                outcome = (len(r.get("clusters", [])) >= 2) == (len(fa.get("clusters", [])) >= 2)
            else:
                outcome = True
            rec["same_outcome"] += int(outcome)
            if not same:
                rec["differ"].append(tid)
        rec["compared_with"] = cl_proc
        out[tag] = rec
    return out or None


def tests():
    out = {}
    for name in ("dials_patch", "cf"):
        f = RESULTS / "tests" / f"{name}.txt"
        if f.exists():
            txt = open_result(f).read()
            m = re.findall(r"(\d+) (passed|failed|skipped|errors?|xfailed|warnings?)", txt)
            out[name] = {k: int(v) for v, k in m}
    return out


def noise_calibration():
    out = {}
    for coll in COLLS:
        for d in (coll, coll + "_amp"):
            f = RESULTS / "noisecal" / f"{d}.json"
            if f.exists():
                out[d] = json.load(open_result(f))
    f = RESULTS / "dev" / "cryo" / "jackknife_cal.json"  # cf/devjack.py
    if f.exists():
        out["cryo_jackknife"] = json.load(open_result(f))
    return out


def limit_robustness(design, proc=PRIMARY):
    """The k90 the data allow (noise one experiment can know) under the rule of amendment 13
    (standardised block sum, whole minority), the rule of amendment 14 (either statistic,
    qualifying sub-groups), the same on the amplitude correlation matrix, and either matrix;
    with the per-minority agreement of the two matrices and the minorities only amplitudes
    separate."""
    from cf.limits3 import classify

    amp = design + "_amp"
    if not (RESULTS / "limits" / f"{amp}_refined").exists():
        return None
    from cf.limits3 import LIMIT_ALPHA

    found = fg_safe_found(design)
    rules = dict(a14=classify(design, found_ids=found),
                 a13=classify(design, found_ids=found, statistic="stouffer", subgroups=False),
                 tol0=classify(design, tol=0.0, found_ids=found),
                 amp=classify(amp, found_ids=found),
                 int_half=classify(design, found_ids=found, alpha=LIMIT_ALPHA / 2),
                 amp_half=classify(amp, found_ids=found, alpha=LIMIT_ALPHA / 2))
    det = lambda q: q is not None and not q["data_limit"] and not q.get("within_single_experiment_noise")
    trials = load_design(design)
    pr = {r["id"]: r for r in score(design, proc)}
    cells = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))
    agree = Counter()
    amp_only = []
    for tid, t in trials.items():
        if t["kind"] != "mixture" or t["k"] < 2:
            continue
        pair = "".join(sorted((t["majority_species"][0], t["minority_species"][0])))
        key = f"{design}|{pair}|{t['w']:g}"
        flags = {name: det(c.get(tid)) for name, c in rules.items()}
        half = (flags.pop("int_half"), flags.pop("amp_half"))
        flags["either"] = half[0] or half[1]  # each matrix at alpha / 2
        for name, fl in flags.items():
            x = cells[key][name][t["k"]]
            x[0] += int(fl)
            x[1] += 1
        if tid in rules["amp"]:
            agree[f"{int(flags['a14'])}{int(flags['amp'])}"] += 1
            if flags["amp"] and not flags["a14"]:
                q = rules["amp"][tid]
                unit = q.get("detectable_unit")
                u = q if unit in (None, "minority") else next(s for s in q["subgroups"] if s["members"] == unit)
                p_u = [x for x in (_finite(u.get("p_noise_any_size")), _finite(u.get("p_noise_any_size_iv"))) if x is not None]
                amp_only.append(dict(id=tid, k=t["k"], w=t["w"], found=bool(pr.get(tid, {}).get("found")),
                                     p=min(p_u) if p_u else None, unit_size=len(u.get("members", [])) or t["k"]))
    k90 = {key: {name: limit(per, 0.9) for name, per in v.items()} for key, v in cells.items()}
    return dict(k90=k90, agree=dict(agree), amp_only=amp_only)


def calibration_transfer(design):
    """Minorities separable only with the noise from their own crystals' repeats (per species):
    their verdict when the counting noise is calibrated instead with one factor learned on the other
    collection (`_xcal`), as a procedure with a frozen calibration could (the review's point 7)."""
    from cf.limits3 import classify

    if not (RESULTS / "limits" / f"{design}_xcal_refined").exists():
        return None
    found = fg_safe_found(design)
    own, x = classify(design, found_ids=found), classify(design + "_xcal", found_ids=found)

    def cat(q):
        if q is None:
            return "not classified"
        if q["within_counting_noise"]:
            return "counting noise"
        if q["within_population_structure"]:
            return "structure"
        if q.get("membership_unresolved"):
            return "membership"
        if q.get("within_single_experiment_noise"):
            return "repeats only"
        return "separable"

    out = Counter(cat(x.get(k)) for k, q in own.items() if cat(q) == "repeats only")
    return dict(out)


def added_false_groups(design, proc=PREREGISTERED):
    """The groups the procedure reports in single-population trials where today's DIALS reports
    none: their sizes, the species, the crystals that recur, and which of those today's DIALS puts
    in another species' cluster when it clusters every public crystal at full data (rt_all)."""
    d = load_design(design)
    st = {r["id"]: r for r in score(design, "today")}
    sp = {r["id"]: r for r in score(design, proc)}
    res = load_results(design, proc)
    sizes, species, crystals = Counter(), Counter(), Counter()
    ps, cohs = [], []
    for tid, r in sp.items():
        if r["kind"] != "single" or not r["false_group"] or st[tid]["false_group"]:
            continue
        species[d[tid]["species"]] += 1
        for l in res[tid].get("detector_log", []):
            if l.get("accepted"):
                ps.append(l["p"])
                cohs.append(l["coherence"])
        groups = sorted(res[tid].get("detector_groups", []), key=len)[:-1]
        for g in groups:
            sizes[len(g)] += 1
            crystals.update(g)
    elsewhere = []
    full = load_results(f"{design}_all", "today").get(f"{design}-all-labelled-w60") if design == "rt" else None
    if full:
        lab = {r["crystal"]: r["species"] for r in manifest()}
        for c in full.get("clusters", []):
            cnt = Counter(lab.get(x) for x in c)
            major = cnt.most_common(1)[0][0]
            elsewhere += [x for x in c if lab.get(x) != major]
    top = [dict(crystal=x, trials=n, species=None, in_other_cluster=x in elsewhere) for x, n in crystals.most_common(8)]
    return dict(sizes=dict(sorted(sizes.items())), species=dict(species), recurring=top,
                n_trials=sum(species.values()), p_min=float(min(ps)) if ps else None,
                p_median=float(np.median(ps)) if ps else None,
                coherence_quartiles=[float(x) for x in np.quantile(cohs, [0.25, 0.75])] if cohs else None)


HELDOUT = "rtbh"  # the third collection (amendment 19)


def heldout(design=HELDOUT, procs=("today",) + JUDGED):
    """The judgement on the untouched collection: per procedure, minorities found, false groups,
    the trials where it reports a group and today's DIALS does not, k90 per cell; and the label
    verification (results/joint/rtbh_verify.json)."""
    if not (RESULTS / "trials" / design / "today").exists():
        return None
    out = dict(procedures={})
    n_trials = len(load_design(design))
    for p in procs:
        if not (RESULTS / "trials" / design / p).exists():
            continue
        rows = score(design, p)
        if len({r["id"] for r in rows}) < n_trials:  # a judgement only on every trial
            continue
        mix = [r for r in rows if r["kind"] == "mixture" and r.get("k", 0) >= 2]
        sing = [r for r in rows if r["kind"] == "single"]
        out["procedures"][p] = dict(
            found=ci(sum(r["found"] for r in mix), len(mix)),
            false_groups=ci(sum(r["false_group"] for r in sing), len(sing)),
            errors=sum(r["error"] for r in rows),
            paired=paired_false_groups(design, b=p) if p != "today" else None,
            detection=detection_summary(rows), cost=cost(design, p))
    f = RESULTS / "joint" / "rtbh_verify.json"
    out["labels"] = json.load(open_result(f)) if f.exists() else None
    out["n_trials"] = len(load_design(design))
    return out


def small_rules():
    """The small-group rule searches (cf/devrules.py, in-sample on both collections): per family,
    the rules searched, those with no single-population trial reporting a group where today's DIALS
    has none, the best of them (most minorities found), and v7-v10 derived exactly from v6's groups."""
    f = RESULTS / "dev" / "smallrules.json"
    if not f.exists():
        return None
    res = json.load(open_result(f))
    out = dict(chosen=res["chosen"])
    for fam in ("precise", "snr", "zfloor"):
        rows = res[fam]
        ok = [r for r in rows if r["added_cryo"] == 0 and r["added_rt"] == 0]
        out[fam] = dict(n=len(rows), fg_free=len(ok),
                        best=max(ok, key=lambda r: r["found_cryo"] + r["found_rt"]) if ok else None)
    return out


def rule_slack(g, family, **p):
    """How far a group of two to four is from a small-group rule's admission: the smallest ratio of
    each condition's value to its threshold (1 or more: admitted). Families as cf/devrules.py."""
    S, se, dec, k = g["S"], g["se"], g["decoy"], g["k"]
    if not (np.isfinite(se) and se > 0):
        return 0.0
    ratio = S / dec if dec > 0 else np.inf
    if family == "precise":  # v8
        parts = [p["max_se"] / se, (S - p["z"] * se) / (p["c0"] * 2 / k), ratio / p["r0"]]
    elif family == "snr":  # v9, v10
        parts = [S / se / p["z0"], ratio / p["r0"]] + ([(S - 2 * se) / (p["c0"] * 2 / k)] if p["c0"] else [])
    else:  # z above a floor
        floor = p["c0"] * 2 / k if p["scale"] else p["c0"]
        parts = [(S - p["z"] * se) / floor, ratio / p["r0"]]
    return float(min(parts))


def small_group_gate(design):
    """v6's groups of two to four (cf/devsmall.py: coherence S, its standard error from DIALS's
    analytic noise, the rest's strongest structure of the size) in single-population trials and in
    mixtures, against the small-group rules of v8 (amendment 20) and v9 (amendment 22, also v10's),
    both through v9's small_group_admitted: groups, the largest S/se, the median se, how many each rule
    keeps; and in single-population trials the within-species group nearest each rule's admission
    (rule_slack) among the trials where today's DIALS reports no group, for these rules and the best
    of each searched family (cf/devrules.py)."""
    from cf.newproc import SMALL_GROUPS, load_module

    files = glob_results(str(RESULTS / "dev" / design / "smallgroups" / "shard_*.jsonl"))
    if not files:
        return None
    admitted = load_module("v9").small_group_admitted
    max_se, z, c0, r0 = SMALL_GROUPS["v8"]
    rules = {"v8": dict(min_signal_to_noise=0.0, z=z, min_coherence=c0, min_ratio=r0, max_se=max_se),
             "v9": SMALL_GROUPS["v9"]}
    d = load_design(design)
    today = {r["id"]: r for r in score(design, "today")}
    out = {}
    for kind in ("single", "mixture"):
        recs = [(r["id"], g) for f in files for l in open_result(f) for r in [json.loads(l)]
                if d.get(r["id"], {}).get("kind") == kind for g in r["groups"]]
        gs = [g for _, g in recs]
        snr = [g["S"] / g["se"] for g in gs if np.isfinite(g["se"]) and g["se"] > 0]
        out[kind] = dict(groups=len(gs), precise=sum(g["se"] <= max_se for g in gs),
                         snr_max=float(max(snr)) if snr else None,
                         se_median=float(np.median([g["se"] for g in gs])) if gs else None,
                         passing={name: sum(bool(admitted(g["S"], g["se"], g["decoy"], g["k"], **kw)) for g in gs)
                                  for name, kw in rules.items()})
        if kind == "single" and gs:
            fams = {"v8": ("precise", dict(max_se=max_se, z=z, c0=c0, r0=r0)),
                    "v9": ("snr", dict(z0=30, c0=0.01, r0=1.5))}
            fams["v11"] = ("union", [fams["v8"], fams["v9"]])  # either condition (v10's is v9's)
            sr = small_rules()
            for fam in ("precise", "snr", "zfloor"):
                if sr and sr[fam]["best"]:
                    fams[f"best_{fam}"] = (fam, {k: v for k, v in sr[fam]["best"].items()
                                                 if k in ("max_se", "z", "z0", "c0", "r0", "scale")})
            nearest = {}
            def slack(g, fam, p):
                return max(rule_slack(g, f, **pp) for f, pp in p) if fam == "union" else rule_slack(g, fam, **p)

            for name, (fam, p) in fams.items():
                counted = [(slack(g, fam, p), tid, g) for tid, g in recs
                           if not today.get(tid, {}).get("false_group")]
                sl, tid, g = max(counted, key=lambda x: x[0])
                nearest[name] = dict(slack=sl, trial=tid, k=g["k"], S=g["S"], se=g["se"], params=str(p),
                                     admitted_where_today_has_a_group=sum(
                                         slack(g2, fam, p) >= 1 for t2, g2 in recs
                                         if today.get(t2, {}).get("false_group")))
            out[kind]["nearest"] = nearest
    return out


def detector_idle(design, proc=PRIMARY):
    """Trials in which the procedure's detector examined no candidate group (no data set with
    enough correlations that have an error, or nothing to scan), per wedge."""
    d = load_design(design)
    out = defaultdict(lambda: [0, 0])
    for tid, r in load_results(design, proc).items():
        if "error" in r or tid not in d:
            continue
        w = f"{d[tid]['w']:g}"
        out[w][0] += int(len(r.get("detector_log", [])) == 0)
        out[w][1] += 1
    return dict(out)


def _finite(x):
    return float(x) if x is not None and np.isfinite(x) else None


# ---------------------------------------------------------------------------------------


def main():
    # These inputs are a fixed snapshot during a report. Reuse repeated reads.
    from functools import lru_cache
    global load_design, load_results, load_dev, score
    load_design = lru_cache(None)(load_design)
    load_results = lru_cache(None)(load_results)
    load_dev = lru_cache(None)(load_dev)
    score = lru_cache(None)(score)
    OUT.mkdir(parents=True, exist_ok=True)
    numbers = dict(commit=commit(), primary=PRIMARY, preregistered=PREREGISTERED, floor_only=FLOOR_ONLY,
                   precise=PRECISE, snr=SNR, takeout=TAKEOUT, judged=list(JUDGED), procedures=list(PROCS),
                   configurations=[dict(name=a, amendment=b, patches=c, does=d) for a, b, c, d in CONFIGURATIONS])
    numbers["data"] = data_summary()
    numbers["fidelity"] = dict(
        multiplex_cryo_full=multiplex_result("cryo_all_default"),
        multiplex_cryo_10_images=multiplex_result("cryo_all_10images"),
        multiplex_cryo_8_images=multiplex_result("cryo_all_8images"),
        dials_test_coords=dials_test_coords(),
        paper=dict(cryo="three species-pure clusters, no outliers, three dimensions, resolution "
                        "filter 1.26 A (section 3.2.1); OPTICS still separates the species down to "
                        "0.8 degrees per crystal, below which data sets become outliers (3.3.3)",
                   rt="on the full room-temperature set (215 data sets, of which 25 per species are "
                      "public): three species-pure clusters of 83, 85 and 28 and 19 outliers (3.2.2)"))
    for d in ("cryo_all", "rt_all"):
        if (RESULTS / "trials" / d / "today").exists():
            numbers["fidelity"][f"today_{d}"] = purity_table(d)
    numbers["noise_calibration"] = noise_calibration()
    for coll in COLLS:
        if not (RESULTS / "design" / f"{coll}.jsonl").exists():
            continue
        refined, strict, cross = refined_limits(coll)
        pre = preregistered_limits(coll)
        numbers[f"{coll}|data_k90"] = data_k90(coll, refined, pre, strict)
        numbers[f"{coll}|data_k90_cross"] = data_k90(coll, cross, pre) if cross else {}
        n_mix = sum(1 for t in load_design(coll).values() if t["kind"] == "mixture" and t["k"] >= 2)
        numbers[f"{coll}|classified"] = [len(refined), n_mix]
        numbers[f"{coll}|limit_robustness"] = limit_robustness(coll)
        numbers[f"{coll}|amendment24"] = amendment24_effect(coll)
        numbers[f"{coll}|calibration_transfer"] = calibration_transfer(coll)
        if (RESULTS / "trials" / coll / PREREGISTERED).exists():
            numbers[f"{coll}|added_false_groups"] = added_false_groups(coll)
        if (RESULTS / "trials" / coll / PRIMARY).exists():
            numbers[f"{coll}|detector_idle"] = detector_idle(coll)
        for proc in PROCS:
            if not (RESULTS / "trials" / coll / proc).exists():
                continue
            rows = score(coll, proc)
            misses, misses_pre = classify_misses(rows, refined, pre)
            misses_strict, _ = classify_misses(rows, strict, pre)
            misses_cross, _ = classify_misses(rows, cross, pre) if cross else ({}, None)
            numbers[f"{coll}|{proc}"] = dict(
                n_trials=len(rows), errors=sum(r["error"] for r in rows),
                detection=detection_summary(rows), false_groups=false_groups(rows),
                isolated=isolated(rows), misses=misses, misses_strict=misses_strict, misses_cross=misses_cross,
                misses_preregistered=misses_pre, cost=cost(coll, proc))
        for p in JUDGED:
            if (RESULTS / "trials" / coll / p).exists():
                numbers[f"{coll}|paired_false_groups_{p}"] = paired_false_groups(coll, b=p)
        if (RESULTS / "trials" / coll / PRIMARY).exists():
            numbers[f"{coll}|paired_false_groups"] = paired_false_groups(coll)
            numbers[f"{coll}|unexplained"] = unexplained(coll, refined)
            numbers[f"{coll}|unexplained_cross"] = unexplained(coll, cross) if cross else None
            ux = numbers[f"{coll}|unexplained_cross"] or []
            gs = [u["v6_group"] for u in ux if u.get("v6_group") and u["v6_group"].get("snr") is not None]
            numbers[f"{coll}|unexplained_summary"] = dict(
                n=len(ux), idle=sum(u["detector_idle"] for u in ux), with_group=len(gs),
                snr=[min(g["snr"] for g in gs), max(g["snr"] for g in gs)] if gs else None,
                se=[min(g["se"] for g in gs), max(g["se"] for g in gs)] if gs else None)
        numbers[f"{coll}|small_group_gate"] = small_group_gate(coll)
        numbers[f"{coll}|dropped_by_cosym"] = dropped_by_cosym(coll)
        numbers[f"{coll}|cause"] = cause_table(coll)
        numbers[f"{coll}|frontier"] = frontier(coll)
    numbers["small_rules"] = small_rules()
    numbers["heldout"] = heldout()
    numbers["equivalence"] = equivalence()
    numbers["tests"] = tests()
    with open_result(OUT / "numbers.json", "w") as fh:
        json.dump(numbers, fh, indent=1, default=str)
    print("wrote", OUT / "numbers.json")


if __name__ == "__main__":
    main()
