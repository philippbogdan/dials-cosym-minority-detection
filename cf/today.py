"""Today's DIALS clustering: the OPTICS significant clusters of dials.correlation_matrix
(DIALS 3.30.1, every default), called the way the command-line program calls it."""

from __future__ import annotations

import copy
import random
import time
import traceback

import numpy as np
from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.array_family import flex
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)


def default_params():
    from dials.command_line.correlation_matrix import phil_scope

    return phil_scope.extract()


def run(expts, joint, names, params=None, nproc: int = 1, keep_arrays: bool = True):
    """Run dials.correlation_matrix's analysis on scaled data of a trial.

    `names[i]` is the crystal of experiment i (experiment identifiers are used to map back).
    Returns a dict with per-crystal labels (cluster index, or -1 for outlier/noise)."""
    if params is None:
        params = default_params()
    params = copy.deepcopy(params)
    params.nproc = nproc
    ident_to_name = {e.identifier: n for e, n in zip(expts, names)}
    if params.seed is not None:  # as dials.command_line.correlation_matrix.run does
        flex.set_random_seed(params.seed)
        np.random.seed(params.seed)
        random.seed(params.seed)
    t0 = time.process_time()
    w0 = time.time()
    reflections = parse_multiple_datasets([copy.deepcopy(joint)])
    experiments, reflections = assign_unique_identifiers(expts, reflections)
    cm = CorrelationMatrix(experiments=experiments, reflections=reflections, params=params)
    try:
        cm.calculate_matrices()
    except Exception as e:  # today's DIALS failed: it reports no clusters
        id_to_name = {i: ident_to_name[ident] for i, ident in cm.ids_to_identifiers_map.items()}
        out = dict(labels={n: -1 for n in names}, clusters=[], outliers=[], dims=None,
                   dim_curve=[], cpu_s=time.process_time() - t0, wall_s=time.time() - w0,
                   used=[id_to_name[i] for i in cm.labels], dials_error=repr(e),
                   dials_traceback=traceback.format_exc()[-3000:])
        if keep_arrays and hasattr(cm.cosym_analysis, "target"):
            t = cm.cosym_analysis.target
            out["_arrays"] = dict(rij=np.array(t.rij_matrix), wij=np.array(t.wij_matrix),
                                  coords=np.zeros((len(names), 0)), order=np.array(out["used"]))
        out["_cm"] = cm
        return out
    cpu = time.process_time() - t0
    wall = time.time() - w0
    id_to_name = {i: ident_to_name[ident] for i, ident in cm.ids_to_identifiers_map.items()}
    labels = {n: -1 for n in names}
    clusters = []
    for c in cm.significant_clusters:
        members = [id_to_name[i] for i in c.labels]
        clusters.append(members)
        for m in members:
            labels[m] = int(c.cluster_id)
    used = [id_to_name[i] for i in cm.labels]
    out = dict(
        labels=labels,
        clusters=clusters,
        outliers=[id_to_name[i] for i in cm.outliers],
        used=used,
        dims=int(cm.cosym_analysis.target.dim),
        dim_curve=[float(v) for v in cm._dimension_optimisation_data.get("functional", [])],
        cpu_s=cpu,
        wall_s=wall,
    )
    if keep_arrays:
        out["_arrays"] = dict(
            rij=np.array(cm.cosym_analysis.target.rij_matrix),
            wij=np.array(cm.cosym_analysis.target.wij_matrix),
            coords=np.array(cm.cosym_analysis.coords),
            order=np.array(used),
        )
    out["_cm"] = cm
    return out
