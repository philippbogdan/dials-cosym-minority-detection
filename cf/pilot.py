"""Data-quality pilot (no procedure is evaluated): for each wedge width, scale all labelled
crystals of a collection together and compute DIALS's pairwise sigma-weighted CC matrix and
effective sample sizes. Summarises within- and between-species CC and neff per wedge, which
fixes the wedge grid of the trial design (PREREGISTRATION.md amendment 1)."""

from __future__ import annotations

import json
import sys

import numpy as np

from cf import design, trialdata
from cf.paths import RESULTS
from cf.run_trials import FULL, cosym_inputs


def rij_for(coll, names, w):
    from dials.algorithms.correlation.analysis import CorrelationMatrix
    from dials.util.multi_dataset_handling import (
        assign_unique_identifiers,
        parse_multiple_datasets,
    )

    from cf.today import default_params

    cd = trialdata.Collection(coll)
    expts, tables = cd.trial(names, None if w >= FULL[coll] else w)
    sexpts, joint, _ = trialdata.scale(expts, tables)
    ident = {e.identifier: n for e, n in zip(expts, names)}
    params = default_params()
    params.nproc = 1
    refl = parse_multiple_datasets([joint])
    sexpts, refl = assign_unique_identifiers(sexpts, refl)
    cm = CorrelationMatrix(sexpts, refl, params)
    cm.cosym_analysis._intialise_target()
    order = [ident[cm.ids_to_identifiers_map[i]] for i in cm.labels]
    t = cm.cosym_analysis.target
    # raw neff: undo the normalisation of the weights (scale so sum == count nonzero)
    ci = cosym_inputs(cm)
    n_unique = [int((ci["ds"] == d).sum()) for d in range(len(order))]
    return order, np.array(t.rij_matrix), np.array(t.wij_matrix), n_unique, cm


def main(coll):
    by = design.crystals_by_species(coll)
    names = by["bovine"] + by["porcine"] + by["human"]
    sp = {n: s for s, v in by.items() for n in v}
    grid = {"cryo": [10, 5, 3, 2, 1.5, 1, 0.7, 0.5, 0.3, 0.2],
            "rt": [60, 30, 20, 10, 5, 3, 2, 1, 0.5, 0.3]}[coll]
    out = []
    for w in grid:
        try:
            order, rij, wij, n_unique, cm = rij_for(coll, names, w)
        except Exception as e:
            out.append(dict(w=w, error=repr(e)))
            print(w, "error", e, flush=True)
            continue
        lab = np.array([sp[n] for n in order])
        iu = np.triu_indices(len(order), 1)
        rec = dict(w=w, n=len(order), d_min=getattr(cm.cosym_analysis, "d_min", None),
                   n_unique_median=float(np.median(n_unique)))
        for a in design.SPECIES:
            for b in design.SPECIES:
                if a > b:
                    continue
                m = (lab[iu[0]] == a) & (lab[iu[1]] == b) | (lab[iu[0]] == b) & (lab[iu[1]] == a)
                rec[f"cc_{a[0]}{b[0]}_mean"] = float(np.mean(rij[iu][m]))
                rec[f"cc_{a[0]}{b[0]}_sd"] = float(np.std(rij[iu][m]))
        rec["wij_mean"] = float(np.mean(wij[iu]))
        out.append(rec)
        print(json.dumps(rec), flush=True)
        np.savez_compressed(RESULTS / "pilot" / f"{coll}_w{w:g}.npz", rij=rij, wij=wij,
                            order=np.array(order))
    (RESULTS / "pilot" / f"{coll}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    (RESULTS / "pilot").mkdir(parents=True, exist_ok=True)
    main(sys.argv[1])
