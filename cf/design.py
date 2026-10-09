"""Trial design (PREREGISTRATION.md): mixture trials with disjoint minority draws, and
single-population trials. Deterministic: the same seed always gives the same trials."""

from __future__ import annotations

from cf.result_io import open_result

import hashlib
import json
import sys

import numpy as np

from cf.paths import RESULTS, manifest

SPECIES = ("bovine", "porcine", "human")
SEED = 20261002

# Fixed in PREREGISTRATION.md amendment 1 (after the data-quality pilot, before any trial).
KS = {"cryo": (1, 2, 3, 4, 5, 6, 8, 10, 12), "rt": (1, 2, 3, 4, 5, 6, 8, 10, 12, 25),
      "rtbh": (1, 2, 3, 4, 5, 6, 8)}  # amendment 19
WEDGES = {"cryo": (10, 1, 0.5, 0.4, 0.3, 0.2), "rt": (60, 20, 10, 5, 2, 1), "rtbh": (20, 10, 5, 2, 1)}
FULL = {"cryo": 10.0, "rt": 60.0, "rtbh": 20.0}
PARTS = {"rtbh": ("all",)}  # eight crystals per species: no halves
SEGMENTS = (0, 1, 2, 3)  # amendment 1: later disjoint wedge segments for single trials


def _rng(*key) -> np.random.Generator:
    h = hashlib.sha256(json.dumps([SEED, *key]).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def crystals_by_species(coll: str) -> dict[str, list[str]]:
    out = {s: [] for s in SPECIES}
    for r in manifest(coll):
        if r["species"] in out:
            out[r["species"]].append(r["crystal"])
    return {s: sorted(v) for s, v in out.items()}


def mixture_trials(coll: str, wedges=None, ks=None):
    wedges = WEDGES[coll] if wedges is None else wedges
    ks = KS[coll] if ks is None else ks
    by = crystals_by_species(coll)
    for a in SPECIES:
        for b in SPECIES:
            if a == b or not by[a] or not by[b]:
                continue
            perm = list(_rng(coll, a, b).permutation(by[b]))
            for k in ks:
                for d in range(len(perm) // k):
                    minority = sorted(perm[d * k : (d + 1) * k])
                    for w in wedges:
                        yield dict(
                            kind="mixture", coll=coll, majority_species=a,
                            minority_species=b, k=k, draw=d, w=w,
                            majority=by[a], minority=minority,
                            id=f"{coll}-mix-{a[0]}{b[0]}-k{k:02d}-d{d:02d}-w{w:g}",
                        )


def single_trials(coll: str, wedges=None):
    wedges = WEDGES[coll] if wedges is None else wedges
    by = crystals_by_species(coll)
    for s in SPECIES:
        if not by[s]:
            continue
        perm = list(_rng(coll, s, "halves").permutation(by[s]))
        h = len(perm) // 2
        sets = {"all": by[s], "half0": sorted(perm[:h]), "half1": sorted(perm[h:])}
        sets = {k: v for k, v in sets.items() if k in PARTS.get(coll, sets)}
        for part, members in sets.items():
            for w in wedges:
                for seg in SEGMENTS:
                    if (seg + 1) * w > FULL[coll]:
                        continue
                    yield dict(
                        kind="single", coll=coll, species=s, part=part, w=w, seg=seg,
                        crystals=members,
                        id=f"{coll}-single-{s[0]}-{part}-w{w:g}" + (f"-s{seg}" if seg else ""),
                    )


def all_trials(coll: str, wedges=None):
    """Three-population sets (the paper's experiments): every labelled crystal of the
    collection, and for rt also with the 19 unclassified crystals."""
    wedges = WEDGES[coll] if wedges is None else wedges
    by = crystals_by_species(coll)
    labelled = by["bovine"] + by["porcine"] + by["human"]
    sets = {"labelled": labelled}
    if coll == "rt":
        sets["with_unclassified"] = labelled + sorted(
            r["crystal"] for r in manifest(coll) if r["species"] == "unclassified")
    for part, members in sets.items():
        for w in wedges:
            yield dict(kind="all", coll=coll, part=part, w=w, crystals=members,
                       id=f"{coll}-all-{part}-w{w:g}")


def trial_crystals(t: dict) -> list[str]:
    return t["majority"] + t["minority"] if t["kind"] == "mixture" else t["crystals"]


def write(coll: str) -> None:
    trials = list(mixture_trials(coll)) + list(single_trials(coll))
    out = RESULTS / "design" / f"{coll}.jsonl"
    with open_result(RESULTS / "design" / f"{coll}_all.jsonl", "w") as fh:
        for t in all_trials(coll):
            fh.write(json.dumps(t) + "\n")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open_result(out, "w") as fh:
        for t in trials:
            fh.write(json.dumps(t) + "\n")
    print(coll, len(trials), "trials ->", out)


if __name__ == "__main__":
    write(sys.argv[1])
