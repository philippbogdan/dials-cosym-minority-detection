"""Render FINDINGS.md, REPORT.md and COMPLETE.md from docs-src/*.md.j2 and
results/report/numbers.json (written by cf.report): every number in those documents comes
from the results files of one commit, recorded at the top of each document."""

from __future__ import annotations

from cf.result_io import open_result

import json
from collections import Counter

import jinja2

from cf.analyze import clopper_pearson
from cf.paths import RESULTS, ROOT

PAIR_NAMES = {"bh": "bovine-human", "bp": "bovine-porcine", "hp": "human-porcine"}
WEDGES = {"cryo": (10, 1, 0.5, 0.4, 0.3, 0.2), "rt": (60, 20, 10, 5, 2, 1), "rtbh": (20, 10, 5, 2, 1)}
KS = {"cryo": (2, 3, 4, 5, 6, 8, 10, 12), "rt": (2, 3, 4, 5, 6, 8, 10, 12, 25), "rtbh": (2, 3, 4, 5, 6, 8)}
MAJORITY = {"cryo": 12, "rt": 25, "rtbh": 8}
FULL = {"cryo": 10, "rt": 60, "rtbh": 20}
MISS_CATS = ("found", "miss: counting noise", "miss: structure inside one population",
             "miss: membership not resolved", "miss: noise known only from repeat measurements",
             "miss: procedure", "miss: not classified")
SHORT = {"found": "found", "miss: counting noise": "counting noise",
         "miss: structure inside one population": "structure inside one species",
         "miss: membership not resolved": "membership not resolved",
         "miss: noise known only from repeat measurements": "noise known only from repeats",
         "miss: procedure": "procedure", "miss: not classified": "not classified"}


def ci(x, n):
    if not n:
        return "-"
    lo, hi = clopper_pearson(x, n)
    return f"{x}/{n} [{lo:.2f}, {hi:.2f}]"


def frac(x, n):
    return f"{x}/{n}" if n else "-"


def kfmt(v, coll):
    return str(v) if v is not None else f"> {KS[coll][-1]}"


def kfrac(v, coll):
    """k90 as a number of crystals and as a fraction of the trial, k / (m + k)."""
    if v is None:
        return f"> {KS[coll][-1]}"
    return f"{v} ({v / (MAJORITY[coll] + v):.0%})"


def k90_table(n, coll, proc, key="k90"):
    d = n.get(f"{coll}|{proc}", {}).get("detection", {})
    rows = ["| species pair | " + " | ".join(f"{w:g} deg" for w in WEDGES[coll]) + " |",
            "|---|" + "---|" * len(WEDGES[coll])]
    for pair, name in PAIR_NAMES.items():
        cells = []
        for w in WEDGES[coll]:
            c = d.get(f"{coll}|{pair}|{w:g}")
            cells.append("-" if c is None else kfrac(c[key], coll))
        rows.append(f"| {name} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def k90_compare(n, coll, procs):
    """Per cell: the data's k90 (repeat-calibrated noise; with the noise one experiment can know;
    and the same across the collections, amendment 18) against each procedure's k90."""
    dk = n.get(f"{coll}|data_k90", {})
    dx = n.get(f"{coll}|data_k90_cross", {}) or {}
    head = ("| species pair | wedge | data (noise from repeats) | data (noise one experiment can know) | "
            "data, across the collections | " + " | ".join(procs) + " |")
    rows = [head, "|---|---|---|---|---|" + "---|" * len(procs)]
    for pair, name in PAIR_NAMES.items():
        for w in WEDGES[coll]:
            key = f"{coll}|{pair}|{w:g}"
            d = dk.get(key)
            if d is None:
                continue
            cells = [kfmt(n.get(f"{coll}|{p}", {}).get("detection", {}).get(key, {}).get("k90"), coll)
                     for p in procs]
            x = dx.get(key, {}).get("k90_single_experiment") if dx.get(key) else None
            rows.append(f"| {name} | {w:g} deg | {kfmt(d['k90'], coll)} | "
                        f"{kfmt(d.get('k90_single_experiment'), coll)} | {kfmt(x, coll) if dx.get(key) else '-'} | "
                        + " | ".join(cells) + " |")
    return "\n".join(rows)


def k90_compare3(n, coll, procs):
    """Per cell, the k90 of each procedure (crystals; fraction of the trial)."""
    rows = ["| species pair | wedge | " + " | ".join(procs) + " |", "|---|---|" + "---|" * len(procs)]
    for pair, name in PAIR_NAMES.items():
        for w in WEDGES[coll]:
            key = f"{coll}|{pair}|{w:g}"
            if key not in n.get(f"{coll}|today", {}).get("detection", {}):
                continue
            cells = [kfrac(n.get(f"{coll}|{p}", {}).get("detection", {}).get(key, {}).get("k90"), coll) for p in procs]
            rows.append(f"| {name} | {w:g} deg | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def gap_cells(n, coll, proc, cross=False):
    """Cells where the procedure's k90 is above the data's (noise one experiment can know; with
    cross, the data's across the collections, amendment 18)."""
    dk = n.get(f"{coll}|data_k90_cross" if cross else f"{coll}|data_k90", {}) or {}
    det = n.get(f"{coll}|{proc}", {}).get("detection", {})
    out = []
    for key, d in dk.items():
        dv = d.get("k90_single_experiment")
        pv = det.get(key, {}).get("k90")
        big = 10**6
        if (pv if pv is not None else big) > (dv if dv is not None else big):
            _, pair, w = key.split("|")
            out.append(f"{PAIR_NAMES[pair]} at {w} deg ({kfmt(pv, coll)} against {kfmt(dv, coll)})")
    return out


def pooled_detection_table(n, coll, proc):
    """Found / trials per wedge and minority size, pooled over the three species pairs, with
    exact 95% intervals."""
    d = n.get(f"{coll}|{proc}", {}).get("detection", {})
    rows = ["| wedge | " + " | ".join(f"k={k}" for k in KS[coll]) + " |",
            "|---|" + "---|" * len(KS[coll])]
    for w in WEDGES[coll]:
        cells = []
        for k in KS[coll]:
            x = m = 0
            for pair in PAIR_NAMES:
                c = d.get(f"{coll}|{pair}|{w:g}", {}).get("per_k", {}).get(str(k))
                if c:
                    x += c[0]
                    m += c[1]
            cells.append(ci(x, m))
        rows.append(f"| {w:g} deg | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def false_group_table(n, coll, procs):
    rows = ["| wedge | " + " | ".join(f"{p}" for p in procs) + " |", "|---|" + "---|" * len(procs)]
    for w in WEDGES[coll]:
        cells = []
        for p in procs:
            c = n.get(f"{coll}|{p}", {}).get("false_groups", {}).get(f"{coll}|{w:g}")
            cells.append(ci(*c["all"]) if c else "-")
        rows.append(f"| {w:g} deg | " + " | ".join(cells) + " |")
    tot = []
    for p in procs:
        x, m = totals(n, coll, p)["false_groups"]
        tot.append(ci(x, m))
    rows.append("| all | " + " | ".join(tot) + " |")
    return "\n".join(rows)


def totals(n, coll, proc):
    d = n.get(f"{coll}|{proc}", {})
    fg = d.get("false_groups", {})
    x = sum(v["all"][0] for v in fg.values())
    m = sum(v["all"][1] for v in fg.values())
    det = d.get("detection", {})
    fx = sum(c[0] for cell in det.values() for c in cell["per_k"].values())
    fn = sum(c[1] for cell in det.values() for c in cell["per_k"].values())
    return dict(false_groups=(x, m), found=(fx, fn))


def miss_totals(n, coll, proc, key="misses"):
    tot = Counter()
    for v in n.get(f"{coll}|{proc}", {}).get(key, {}).values():
        tot.update(v)
    return {SHORT.get(k, k): tot.get(k, 0) for k in MISS_CATS if tot.get(k, 0) or k == "miss: procedure"}


def misses_by_procedure(n, coll, procs, key="misses"):
    cats = [c for c in MISS_CATS if c != "miss: not classified"]
    rows = ["| procedure | " + " | ".join(SHORT[c] for c in cats) + " |", "|---|" + "---|" * len(cats)]
    for p in procs:
        tot = Counter()
        for v in n.get(f"{coll}|{p}", {}).get(key, {}).values():
            tot.update(v)
        rows.append(f"| {p} | " + " | ".join(str(tot.get(c, 0)) for c in cats) + " |")
    return "\n".join(rows)


def misses_table(n, coll, proc):
    d = n.get(f"{coll}|{proc}", {}).get("misses", {})
    cats = [c for c in MISS_CATS if c != "miss: not classified"]
    rows = ["| species pair | wedge | " + " | ".join(SHORT[c] for c in cats) + " |",
            "|---|---|" + "---|" * len(cats)]
    for pair, name in PAIR_NAMES.items():
        for w in WEDGES[coll]:
            c = d.get(f"{coll}|{pair}|{w:g}")
            if c:
                rows.append(f"| {name} | {w:g} deg | " + " | ".join(str(c.get(k, 0)) for k in cats) + " |")
    return "\n".join(rows)


def cause_table(n, coll, variants=("default", "floor2", "heuristic", "unit", "dims+1", "dims5", "floor2_dims+1")):
    """Today's DIALS replayed with one stage changed: minorities of 2-4 found at full data
    and at the next wedge, and false groups over all single-population trials; variants and
    columns that were not replayed are left out."""
    c = n.get(f"{coll}|cause")
    if not c:
        return "(not measured)"
    w1, w2 = WEDGES[coll][0], WEDGES[coll][1]
    cols = [(f"k=2-4 found, {w1:g} deg", w1, (2, 3, 4)), (f"k=5-12 found, {w1:g} deg", w1, (5, 6, 8, 10, 12)),
            (f"k=2-4 found, {w2:g} deg", w2, (2, 3, 4))]
    body = []
    for v in variants:
        if v not in c:
            continue
        det = c[v]["detection"]

        def s(w, ks):
            x = sum(det.get(f"{w:g}|{k}", [0, 0])[0] for k in ks)
            m = sum(det.get(f"{w:g}|{k}", [0, 0])[1] for k in ks)
            return frac(x, m) if m else None

        fg = c[v]["false_groups"]
        fx = sum(a[0] for a in fg.values())
        fm = sum(a[1] for a in fg.values())
        body.append((v, [s(w, ks) for _, w, ks in cols], ci(fx, fm) if fm else None))
    keep = [i for i in range(len(cols)) if any(r[1][i] is not None for r in body)]
    body = [r for r in body if any(r[1][i] is not None for i in keep) or r[2] is not None]
    rows = ["| variant | " + " | ".join(cols[i][0] for i in keep) + " | false groups (single-population trials replayed) |",
            "|---|" + "---|" * (len(keep) + 1)]
    for v, cells, fg in body:
        rows.append(f"| {v} | " + " | ".join(cells[i] or "-" for i in keep) + f" | {fg or '-'} |")
    return "\n".join(rows)


def separability_table(n, coll):
    c = (n.get(f"{coll}|cause") or {}).get("separability", {})
    rows = ["| wedge | minority size | silhouette of the true split (median) | minority's nearest neighbours all minority | dimensions (median) | trials |",
            "|---|---|---|---|---|---|"]
    for w in WEDGES[coll][:3]:
        for k in (2, 3, 4):
            d = c.get(f"{w:g}|{k}")
            if d:
                rows.append(f"| {w:g} deg | {k} | {d['silhouette_median']:.2f} | {d['nn_pure']:.0%} | "
                            f"{d['dims_median']:g} | {d['n']} |")
    return "\n".join(rows)


def frontier_table(n, coll, kinds=("frozen", "ablation", "development")):
    d = n.get(f"{coll}|frontier", {})
    rows = ["| procedure | kind | minorities found | false groups | median CPU s per trial | non-dominated |",
            "|---|---|---|---|---|---|"]
    for name, v in d.items():
        if v.get("incomplete") or v.get("kind") not in kinds:
            continue
        c = v["cpu_s_median"]
        rows.append(f"| {name} | {v.get('kind')} | {frac(v['found'], v['n_mix'])} | "
                    f"{frac(v['false_groups'], v['n_single'])} | {'%.2f' % c if c is not None else '-'} | "
                    f"{'yes' if v.get('non_dominated') else ''} |")
    return "\n".join(rows)


def noise_table(n, coll):
    d = n.get("noise_calibration", {}).get(coll, {})
    da = n.get("noise_calibration", {}).get(coll + "_amp", {})
    f2 = lambda x: "%.2f" % x if x is not None else "-"
    rows = ["| wedge | species | bootstrap over reflections | measurement errors only | analytic (Fisher z) | "
            "amplitude correlations: analytic (Fisher z) |",
            "|---|---|---|---|---|---|"]
    for key, v in sorted(d.items(), key=lambda kv: (-float(kv[0].split("|")[0]), kv[0].split("|")[1])):
        w, sp = key.split("|")
        rows.append(f"| {w} deg | {sp} | {f2(v.get('observed_over_reflection'))} | {f2(v.get('observed_over_parametric'))} | "
                    f"{f2(v.get('observed_over_analytic'))} | {f2(da.get(key, {}).get('observed_over_analytic'))} |")
    return "\n".join(rows)


def misses_by_wedge(n, coll, proc, cat="miss: procedure", key="misses"):
    """The minorities of one miss category, by wedge (largest first): [(wedge, count)]."""
    c = Counter()
    for cell, v in n.get(f"{coll}|{proc}", {}).get(key, {}).items():
        c[float(cell.split("|")[2])] += v.get(cat, 0)
    return [(f"{w:g}", c[w]) for w in sorted(c, reverse=True) if c[w]]


def present(n, coll):
    """The procedures measured on a collection."""
    return [p for p in n.get("procedures", []) if f"{coll}|{p}" in n]


def heldout_table(n):
    """The untouched collection (amendment 19): per procedure, found, false groups, paired."""
    h = n.get("heldout")
    if not h:
        return "(not yet measured)"
    rows = ["| procedure | minorities found (k >= 2) | single-population trials with a false group | "
            "of them where today's DIALS reports none | where only today's does |", "|---|---|---|---|---|"]
    for p, v in h["procedures"].items():
        pa = v.get("paired") or {}
        rows.append(f"| {p} | {ci(*v['found'][:2]) if v['found'] else '-'} | {ci(*v['false_groups'][:2]) if v['false_groups'] else '-'} | "
                    f"{pa.get('only_b', '-') if pa else '-'} | {pa.get('only_a', '-') if pa else '-'} |")
    return "\n".join(rows)


def heldout_k90(n):
    h = n.get("heldout")
    if not h:
        return ""
    procs = list(h["procedures"])
    rows = ["| species pair | wedge | " + " | ".join(procs) + " |", "|---|---|" + "---|" * len(procs)]
    keys = sorted({k for v in h["procedures"].values() for k in v["detection"]},
                  key=lambda k: -float(k.split("|")[2]))
    for key in keys:
        _, pair, w = key.split("|")
        cells = [kfrac(h["procedures"][p]["detection"].get(key, {}).get("k90"), "rtbh") for p in procs]
        rows.append(f"| {PAIR_NAMES[pair]} | {w} deg | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def robustness_table(n, coll, proc=None):
    """The k90 the data allow under each limit rule and matrix, with the procedure's, for the
    cells where any of them differ."""
    r = n.get(f"{coll}|limit_robustness")
    if not r:
        return "(not computed)"
    proc = proc or n["primary"]
    det = n.get(f"{coll}|{proc}", {}).get("detection", {})
    rows = ["| species pair | wedge | amendment 13 (whole minority, block sum) | amendment 14 (the rule) | "
            "the rule, no tolerance in the structure test | amplitude correlations | either matrix (each at "
            f"half the level) | {proc} |", "|---|---|---|---|---|---|---|---|"]
    for key in sorted(r["k90"], key=lambda k: (k.split("|")[1], -float(k.split("|")[2]))):
        v = r["k90"][key]
        pv = det.get(key, {}).get("k90")
        if len({v.get("a13"), v.get("a14"), v.get("tol0"), v.get("amp"), v.get("either"), pv}) == 1:
            continue
        _, pair, w = key.split("|")
        rows.append(f"| {PAIR_NAMES[pair]} | {w} deg | {kfmt(v.get('a13'), coll)} | {kfmt(v.get('a14'), coll)} | "
                    f"{kfmt(v.get('tol0'), coll)} | {kfmt(v.get('amp'), coll)} | {kfmt(v.get('either'), coll)} | "
                    f"{kfmt(pv, coll)} |")
    if len(rows) == 2:
        return "(the same in every cell)"
    return "\n".join(rows)


def purity_rows(rows):
    out = ["| trial | crystals | wedge | clusters (species counts) | outliers | dimensions |", "|---|---|---|---|---|---|"]
    for r in rows:
        if r.get("error"):  # DIALS failed on the trial: it reports nothing
            err = str(r["error"]).split("(")[0]
            out.append(f"| {r['id']} | {r['n']} | {r['w']:g} deg | DIALS failed ({err}) | - | - |")
            continue
        cl = "; ".join(", ".join(f"{k} {v}" for k, v in c.items()) for c in r["clusters"]) or "none"
        out.append(f"| {r['id']} | {r['n']} | {r['w']:g} deg | {cl} | {r['n_outliers']} | {r.get('dims')} |")
    return "\n".join(out)


def coords_summary(n):
    d = n.get("fidelity", {}).get("dials_test_coords") or {}
    rows = ["| DIALS reference coordinates | data sets | expected (paper) | initial min_samples | runs identical to expected | V-measure | ARI |",
            "|---|---|---|---|---|---|---|"]
    for name, rec in d.items():
        exp = ", ".join(f"{('outliers' if k == '-1' else 'cluster ' + k)} {v}" for k, v in rec["expected"].items())
        for label, run in rec.items():
            if not isinstance(run, dict) or "runs" not in run:
                continue
            ident = sum(r["identical"] for r in run["runs"])
            vm = min(r["v_measure"] for r in run["runs"]), max(r["v_measure"] for r in run["runs"])
            ari = min(r["ari"] for r in run["runs"]), max(r["ari"] for r in run["runs"])
            rows.append(f"| {name} ({label}) | {rec['n']} | {exp} | {run['initial_min_samples']} | "
                        f"{ident}/{len(run['runs'])} | {vm[0]:.3f}-{vm[1]:.3f} | {ari[0]:.3f}-{ari[1]:.3f} |")
    return "\n".join(rows)


def clusters_text(clusters):
    """[{'bovine': 12}, {'porcine': 11, 'human': 1}] -> '12 bovine; 11 porcine + 1 human'."""
    return "; ".join(" + ".join(f"{v} {k}" for k, v in c.items()) for c in clusters) or "none"


def coords_changes(n, key="data_4 (room temperature, 215)"):
    """Data sets whose label differs from the paper's on DIALS's reference coordinates (first
    run with min_samples as in DIALS's test), as 'expected -> obtained: count'."""
    rec = (n.get("fidelity", {}).get("dials_test_coords") or {}).get(key, {})
    run = rec.get("min_samples as in DIALS's test", {}).get("runs", [{}])[0]
    cont = run.get("contingency", {})
    moved = {k: v for k, v in cont.items() if k.split("->")[0] != k.split("->")[1]}
    same = sum(v for k, v in cont.items() if k.split("->")[0] == k.split("->")[1])
    return dict(moved=moved, same=same, total=sum(cont.values()))


def main():
    numbers = json.load(open_result(RESULTS / "report" / "numbers.json"))
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(ROOT / "docs-src")),
                             undefined=jinja2.StrictUndefined, keep_trailing_newline=True)
    env.filters["clusters"] = clusters_text
    ctx = dict(n=numbers, ci=ci, frac=frac, kfmt=kfmt, kfrac=kfrac, k90_table=k90_table,
               k90_compare=k90_compare, gap_cells=gap_cells, totals=totals, miss_totals=miss_totals,
               misses_by_procedure=misses_by_procedure, pooled_detection_table=pooled_detection_table,
               false_group_table=false_group_table, misses_table=misses_table, cause_table=cause_table,
               separability_table=separability_table, frontier_table=frontier_table,
               noise_table=noise_table, purity_rows=purity_rows, coords_summary=coords_summary,
               coords_changes=coords_changes, robustness_table=robustness_table, k90_compare3=k90_compare3,
               heldout_table=heldout_table, heldout_k90=heldout_k90, present=present,
               misses_by_wedge=misses_by_wedge,
               PAIR_NAMES=PAIR_NAMES, WEDGES=WEDGES, FULL=FULL, MAJORITY=MAJORITY)
    for name in ("FINDINGS.md", "REPORT.md", "COMPLETE.md"):
        src = ROOT / "docs-src" / f"{name}.j2"
        if src.exists():
            (ROOT / name).write_text(env.get_template(f"{name}.j2").render(**ctx))
            print("rendered", name)


if __name__ == "__main__":
    main()
