"""Assemble the data of one trial: chosen crystals of a collection, each cut to its first
`w` degrees, from the jointly symmetrized (consistently indexed) per-crystal files, then
scaled together with dials.scale defaults, as xia2.multiplex does before clustering."""

from __future__ import annotations

import copy
import logging
import time

from dials.array_family import flex
from dials.util.slice import slice_experiments, slice_reflections
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory

from cf.paths import DATA

logger = logging.getLogger(__name__)


class Collection:
    """Per-crystal symmetrized data of one collection, loaded once and copied per trial."""

    def __init__(self, coll: str):
        self.coll = coll
        self.dir = DATA / "joint" / coll / "crystals"
        self._cache: dict[str, tuple] = {}

    def _load(self, name: str):
        if name not in self._cache:
            el = ExperimentListFactory.from_json_file(
                str(self.dir / f"{name}.expt"), check_format=False
            )
            refl = flex.reflection_table.from_file(str(self.dir / f"{name}.refl"))
            self._cache[name] = (el, refl)
        return self._cache[name]

    def trial(self, names: list[str], wedge_deg: float | None, segment: int = 0):
        """Experiments and per-crystal reflection tables for a trial, in `names` order.

        `wedge_deg=None` keeps the full sweep; otherwise images of wedge segment `segment`
        (0: the first `wedge_deg` degrees). Ids are 0..n-1 in the order of `names`."""
        expts = ExperimentList()
        tables = []
        for i, name in enumerate(names):
            el, refl = self._load(name)
            el = copy.deepcopy(el)
            refl = copy.deepcopy(refl)
            scan = el[0].scan
            if wedge_deg is not None:
                osc = scan.get_oscillation()[1]
                first, last = scan.get_image_range()
                n_img = max(1, int(round(wedge_deg / osc)))
                start = first + segment * n_img
                if start > last:
                    raise ValueError(f"{name}: segment {segment} of {wedge_deg} deg is beyond the sweep")
                rng = (start, min(start + n_img - 1, last))
                el = slice_experiments(el, [rng])
                refl = slice_reflections(refl, [rng])
            ident = el[0].identifier
            refl["id"] = flex.int(len(refl), i)
            for k in list(refl.experiment_identifiers().keys()):
                del refl.experiment_identifiers()[k]
            refl.experiment_identifiers()[i] = ident
            expts.append(el[0])
            tables.append(refl)
        return expts, tables


def scale(expts, tables, nproc: int = 1):
    """dials.scale with every default (as `dials.scale` on the command line)."""
    from dials.command_line.scale import phil_scope, run_scaling

    params = phil_scope.extract()
    params.scaling_options.nproc = nproc
    # dials.scale seeds flex's generator once, at import (scaler.py); reseed so that every
    # trial is scaled exactly as a fresh dials.scale process would scale it.
    flex.set_random_seed(42)
    t0 = time.process_time()
    scaled_expts, joint = run_scaling(params, copy.deepcopy(expts), copy.deepcopy(tables))
    return scaled_expts, joint, time.process_time() - t0
