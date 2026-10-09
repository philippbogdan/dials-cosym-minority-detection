"""Verification of the rtbh labels (PREREGISTRATION.md amendment 19), before any trial: the 16 rtbh
crystals and the labelled room-temperature bovine and human crystals, indexed together by dials.cosym
(pipeline/verify_rtbh.sh), each room-temperature crystal cut to its first 20 degrees, scaled together
with dials.scale defaults; DIALS's pairwise correlation matrix (dials.correlation_matrix, as today's
DIALS computes it); for every rtbh crystal its mean correlation with each labelled species. Writes
results/joint/rtbh_verify.json; exits 1 if any crystal correlates more with the other species."""

from __future__ import annotations

import copy
import json
import sys

import numpy as np
from dials.array_family import flex
from dials.util.slice import slice_experiments, slice_reflections
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory

from cf import today, trialdata
from cf.paths import DATA, RESULTS, manifest


def main():
    j = DATA / "joint" / "rtbh_verify"
    expts = ExperimentListFactory.from_json_file(str(j / "symmetrized.expt"), check_format=False)
    refl = flex.reflection_table.from_file(str(j / "symmetrized.refl"))
    rows = {(r["collection"], r["crystal"]): r for r in manifest()}
    ident = {}
    for (coll, name), r in rows.items():
        if coll in ("rt", "rtbh"):
            p = DATA / "proc" / coll / name / "integrated.expt"
            if p.exists():
                ident[ExperimentListFactory.from_json_file(str(p), check_format=False)[0].identifier] = (coll, name)
    el, tables, names = ExperimentList(), [], []
    for i, e in enumerate(expts):
        coll, name = ident[e.identifier]
        sub = refl.select_on_experiment_identifiers([e.identifier])
        sub.reset_ids()
        one = ExperimentList()
        one.append(copy.deepcopy(e))
        if coll == "rt":  # its first 20 degrees, as long as an rtbh sweep
            first, _ = one[0].scan.get_image_range()
            osc = one[0].scan.get_oscillation()[1]
            rng = (first, first + int(round(20 / osc)) - 1)
            one = slice_experiments(one, [rng])
            sub = slice_reflections(sub, [rng])
        sub["id"] = flex.int(len(sub), i)
        for k in list(sub.experiment_identifiers().keys()):
            del sub.experiment_identifiers()[k]
        sub.experiment_identifiers()[i] = one[0].identifier
        el.append(one[0])
        tables.append(sub)
        names.append(f"{coll}:{name}")
    sexpts, joint, _ = trialdata.scale(el, tables, nproc=8)
    res = today.run(sexpts, joint, names, nproc=8)
    rij = np.array(res["_arrays"]["rij"])
    order = [str(x) for x in res["_arrays"]["order"]]
    if rij.shape[0] != len(order):  # cosym dropped a data set: rows are the lattices still present
        rows_present = np.unique(np.asarray(res["_cm"].cosym_analysis.target._lattice_ids))
        order = [order[i] for i in rows_present]
    species = {f"{c}:{n}": rows[(c, n)]["species"] for c, n in rows}
    ref = {s: [i for i, x in enumerate(order) if x.startswith("rt:") and species[x] == s] for s in ("bovine", "human")}
    out, bad = {}, []
    for i, x in enumerate(order):
        if not x.startswith("rtbh:"):
            continue
        m = {s: float(np.mean([rij[i, j] for j in idx if rij[i, j] != 0])) for s, idx in ref.items()}
        label = species[x]
        other = "human" if label == "bovine" else "bovine"
        ok = m[label] > m[other]
        out[x.split(":")[1]] = dict(label=label, mean_cc=m, agrees=ok)
        if not ok:
            bad.append(x)
    # the matrix itself, and mean correlations within and between the two deposit batches (image
    # numbers 1405xx-1408xx and 1434xx-1438xx) and each labelled species: whether the batches are two
    # populations at all, whichever their labels
    groups = {"batch_1": [i for i, x in enumerate(order) if x.startswith("rtbh:") and x.split(":")[1] < "image_1430"],
              "batch2": [i for i, x in enumerate(order) if x.startswith("rtbh:") and x.split(":")[1] > "image_1430"],
              "rt_bovine": ref["bovine"], "rt_human": ref["human"]}

    def block(a, b):
        v = [rij[i, j] for i in groups[a] for j in groups[b] if i != j and rij[i, j] != 0]
        return float(np.mean(v)) if v else None

    blocks = {f"{a}|{b}": block(a, b) for a in groups for b in groups if a <= b}
    (RESULTS / "joint").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(RESULTS / "joint" / "rtbh_verify_matrix.npz", rij=rij, order=np.array(order))
    (RESULTS / "joint" / "rtbh_verify.json").write_text(json.dumps(dict(crystals=out, disagree=bad, blocks=blocks), indent=1))
    print(json.dumps(blocks, indent=1))
    print(json.dumps(out, indent=1))
    print("disagree:", bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
