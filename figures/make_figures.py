#!/usr/bin/env python3
"""Draw the README figures from the recorded results, without DIALS or report writes.

Run: python figures/make_figures.py
Dependencies: numpy, scipy, matplotlib, Pillow. No fitting or diffraction reruns.
Trial scoring and retry precedence reproduce cf/analyze.py. CSVs retain every
plotted count and coordinate diagnostic; manifest.json records input checksums.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from matplotlib.lines import Line2D
from PIL import Image
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
INPUTS = {}
AUDIT = {}
BLUE = "#0072B2"
ORANGE = "#A64B00"
GREY = "#555555"
METHODS = ("today", "v6", "v11")
COLOURS = (GREY, ORANGE, BLUE)
MARKERS = ("o", "^", "s")
LABELS = ("DIALS 3.30.1", "Earlier v6", "Candidate v11")


def read_rows(path):
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def records(directory):
    # Match cf.result_io: prefer an uncompressed shard if both forms exist.
    paths = {p.name.removesuffix(".gz"): p for p in directory.glob("shard_*.jsonl.gz")}
    paths.update({p.name: p for p in directory.glob("shard_*.jsonl")})
    return [row for name in sorted(paths) for row in read_rows(paths[name])]


def design(collection):
    return {row["id"]: row for row in read_rows(ROOT / "results/design" / f"{collection}.jsonl")}


def found(trial, clusters):
    minority = set(trial["minority"])
    return any(
        len(minority.intersection(group)) >= max(2, len(minority) / 2)
        and len(minority.intersection(group)) > len(group) / 2
        for group in clusters
    )


def interval(hits, total):
    # Nominal exact binomial intervals, as in the report. Reused crystals mean
    # these are descriptive, not independent-collection uncertainty bounds.
    return (
        0 if hits == 0 else float(beta.ppf(0.025, hits, total - hits + 1)),
        1 if hits == total else float(beta.ppf(0.975, hits + 1, total - hits)),
    )


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export(fig, stem):
    fig.savefig(OUT / f"{stem}.svg", facecolor="white", metadata={"Date": None})
    fig.savefig(OUT / f"{stem}.png", dpi=180, facecolor="white")
    plt.close(fig)


def cause_figure():
    trials = design("cryo")
    all_rows = records(ROOT / "results/cause/cryo")
    assert len({r["id"] for r in all_rows}) == len(all_rows)
    rows = [r for r in all_rows if "error" not in r and
            trials[r["id"]]["kind"] == "mixture" and
            trials[r["id"]]["w"] == 10 and 2 <= trials[r["id"]]["k"] <= 4]
    assert len(rows) == 78
    variants = ("default", "dims+1", "dims5", "unit", "heuristic", "floor2", "floor2_dims+1")
    names = ("Default", "One extra dimension", "Five dimensions", "Unit-length coordinates",
             "Fixed min_samples heuristic", "min_samples floor = 2", "Floor = 2 + extra dimension")
    counts = []
    for variant in variants:
        subset = [r for r in rows if variant in r]
        null = [r for r in all_rows if "error" not in r and variant in r
                and trials[r["id"]]["kind"] == "single"]
        counts.append(dict(variant=variant, found=sum(r[variant]["found"] for r in subset),
                           trials=len(subset), false_groups=sum(r[variant]["false_group"] for r in null),
                           single_population_trials=len(null)))
    assert [r["found"] for r in counts] == [0, 0, 0, 0, 0, 73, 73]
    assert counts[5]["false_groups"] == 104 and counts[5]["single_population_trials"] == 184
    write_csv("cause_counts.csv", counts)
    diagnostics = [dict(id=r["id"], k=trials[r["id"]]["k"],
                        silhouette_true=r["silhouette_true"],
                        minority_nn_pure=r["minority_nn_pure"], dims=r["dims"])
                   for r in sorted(rows, key=lambda r: r["id"])]
    write_csv("coordinate_diagnostics.csv", diagnostics)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1, 1.45]})
    fig.subplots_adjust(left=0.075, right=0.97, bottom=0.19, top=0.79, wspace=0.82)
    fig.suptitle("Small groups are separated before OPTICS, but its size floor blocks them", x=0.04,
                 ha="left", fontsize=14, weight="bold", y=0.98)
    fig.text(0.04, 0.875, "Cryogenic insulin, 10 degrees per crystal. All 78 trials with minorities of 2 to 4 crystals.",
             fontsize=10, color=GREY)
    ax = axes[0]
    rng = np.random.default_rng(20261009)
    for position, k in enumerate((2, 3, 4)):
        group = [r for r in diagnostics if r["k"] == k]
        values = np.array([r["silhouette_true"] for r in group])
        ax.scatter(values, position + rng.uniform(-0.13, 0.13, len(group)), s=22,
                   color=BLUE, alpha=0.75, edgecolors="none")
        ax.plot(np.median(values), position, "|", color="black", markersize=23, markeredgewidth=2)
        pure = sum(r["minority_nn_pure"] for r in group)
        ax.text(-0.94, position + 0.28, f"{pure}/{len(group)} nearest-neighbour checks pass", fontsize=8.5)
    ax.set(yticks=range(3), yticklabels=("2 crystals", "3 crystals", "4 crystals"),
           xlim=(-1, 1.06), ylim=(2.5, -0.45), xticks=(-1, 0, 1),
           xlabel="Silhouette of the true split (cosine distance)")
    ax.set_title("a  Separation in the fitted coordinates", loc="left", fontsize=10, pad=12)
    ax.axvline(0, color="#CCCCCC", linewidth=0.8, zorder=0)
    ax.text(0, -0.21, "Each dot is one trial. Black marks are medians.", transform=ax.transAxes,
            fontsize=8.5, color=GREY)

    ax = axes[1]
    for i, row in enumerate(counts):
        colour = BLUE if row["variant"].startswith("floor2") else GREY
        ax.barh(i, 100 * row["found"] / row["trials"], height=0.48, color=colour)
        ax.plot(0, i, "|", color=colour, markersize=10)
        ax.text(100 * row["found"] / row["trials"] + 2, i,
                f"{row['found']}/{row['trials']}", va="center", fontsize=9)
    ax.set(yticks=range(len(names)), yticklabels=names, ylim=(len(names)-0.4, -0.6),
           xlim=(0, 112), xticks=(0, 50, 100), xlabel="Minorities found (%)")
    ax.set_title("b  Change one stage of the saved-matrix replay", loc="left", fontsize=10, pad=12)
    export(fig, "cause")


def comparison_figure():
    summaries = []
    expected = {
        "cryo": ((106, 0), (416, 0), (252, 0)),
        "rt": ((416, 11), (971, 33), (811, 11)),
        "rtbh": ((26, 0), (85, 0), (58, 0)),
    }
    for collection in expected:
        trials = design(collection)
        for index, method in enumerate(METHODS):
            raw = records(ROOT / f"results/trials/{collection}/{method}")
            # Earlier failed attempts precede retries. The final record per ID is
            # authoritative, exactly as in cf.analyze.load_results.
            final = {r["id"]: r for r in raw}
            assert set(final) == set(trials)
            AUDIT[f"{collection}/{method}"] = dict(raw_rows=len(raw), final_trials=len(final),
                                                   errors=sum("error" in r for r in final.values()))
            counts = {"detection": [0, 0], "false_groups": [0, 0]}
            for tid, result in final.items():
                t = trials[tid]
                clusters = result.get("clusters", []) if "error" not in result else []
                if t["kind"] == "mixture" and t["k"] >= 2:
                    kind, hit = "detection", found(t, clusters)
                elif t["kind"] == "single":
                    kind, hit = "false_groups", len(clusters) >= 2
                else:
                    continue
                counts[kind][0] += int(hit)
                counts[kind][1] += 1
            assert tuple(counts[k][0] for k in ("detection", "false_groups")) == expected[collection][index]
            for kind, (hits, total) in counts.items():
                lo, hi = interval(hits, total)
                summaries.append(dict(collection=collection, method=method, metric=kind,
                                      hits=hits, trials=total, fraction=hits/total,
                                      nominal_95_low=lo, nominal_95_high=hi))
    write_csv("comparison_counts.csv", summaries)
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.1), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.12, top=0.77, hspace=0.32, wspace=0.22)
    fig.suptitle("More minorities found, with a small held-out test of false groups",
                 x=0.04, ha="left", fontsize=14, weight="bold", y=0.985)
    handles = [Line2D([], [], color=c, marker=m, linestyle="none", label=l)
               for c, m, l in zip(COLOURS, MARKERS, LABELS)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.035, 0.94), frameon=False,
               ncol=3, columnspacing=2.5)
    for col, (collection, title) in enumerate((
        ("cryo", "Cryo\nUsed to choose v11"),
        ("rt", "Room temperature\nUsed to choose v11"),
        ("rtbh", "Third collection\nHeld out for v11"),
    )):
        axes[0, col].set_title(title, fontsize=11, pad=13)
        for row, metric in enumerate(("detection", "false_groups")):
            ax = axes[row, col]
            for i, (method, colour, marker) in enumerate(zip(METHODS, COLOURS, MARKERS)):
                r = next(r for r in summaries if r["collection"] == collection
                         and r["method"] == method and r["metric"] == metric)
                y = 100*r["fraction"]
                ax.errorbar(i, y, yerr=[[y-100*r["nominal_95_low"]],
                                      [100*r["nominal_95_high"]-y]],
                            fmt=marker, color=colour, markersize=7, capsize=3, linewidth=1.3)
                ax.text(i, 100*r["nominal_95_high"] + (3 if row == 0 else 1.1),
                        f"{r['hits']}/{r['trials']}", ha="center", fontsize=9)
            ax.set(xlim=(-0.55, 2.55), xticks=range(3), xticklabels=("DIALS", "v6", "v11"))
            ax.set_ylim((-3, 102) if row == 0 else (-1, 30))
            ax.set_yticks((0, 25, 50, 75, 100) if row == 0 else (0, 10, 20, 30))
            ax.grid(axis="y", color="#E8E8E8", linewidth=0.7, zorder=0)
    axes[0, 0].set_ylabel("Minorities found (%)")
    axes[1, 0].set_ylabel("Trials with a false group (%)")
    fig.text(0.075, 0.025, "Whiskers: nominal 95% exact binomial intervals. Trials reuse crystals. See caption.",
             fontsize=9, color=GREY)
    export(fig, "comparison")


def main():
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelsize": 9, "axes.linewidth": 0.7,
                         "svg.fonttype": "none", "svg.hashsalt": "cosym-readme-20261009",
                         "savefig.facecolor": "white"}):
        cause_figure()
        comparison_figure()
    outputs = {}
    for p in sorted(OUT.glob("*.png")):
        with Image.open(p) as img:
            outputs[p.name] = dict(pixels=img.size, mode=img.mode, dpi=img.info.get("dpi"),
                                   opaque=img.mode != "RGBA" or img.getchannel("A").getextrema() == (255, 255))
    manifest = dict(inputs=INPUTS, scoring_audit=AUDIT, outputs=outputs,
                    versions=dict(python=platform.python_version(), numpy=np.__version__,
                                  scipy=scipy.__version__, matplotlib=matplotlib.__version__),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    transformations=[
                        "Final trial record per ID in sorted shard order; failed trials report no clusters.",
                        "Detection: k >= 2. A cluster holds at least max(2,k/2) minority crystals, which outnumber its other members.",
                        "False group: two or more clusters in a single-species trial.",
                        "Cause: successful cryo replays at 10 degrees, k=2,3,4; all 78 retained.",
                        "Coordinate dots use deterministic vertical jitter only; seed 20261009.",
                        "Comparison pools all tested wedges and minority sizes; no weighting or smoothing.",
                        "Nominal Clopper-Pearson 95% intervals do not correct dependence from reused crystals.",
                    ])
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"outputs": outputs, "scoring_audit": AUDIT}, indent=2))


if __name__ == "__main__":
    main()
