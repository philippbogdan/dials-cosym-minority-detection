"""Split the jointly symmetrized (dials.cosym) data of a collection back into one
experiment/reflection file pair per crystal, keyed by crystal name, and record the
reindexing applied to each crystal."""

from __future__ import annotations

import json
import sys

from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory

from cf.paths import DATA, RESULTS, manifest


def main(coll: str) -> None:
    j = DATA / "joint" / coll
    expts = ExperimentListFactory.from_json_file(str(j / "symmetrized.expt"), check_format=False)
    refl = flex.reflection_table.from_file(str(j / "symmetrized.refl"))
    ident_to_name = {}
    for row in manifest(coll):
        p = DATA / "proc" / coll / row["crystal"] / "integrated.expt"
        if p.exists():
            el = ExperimentListFactory.from_json_file(str(p), check_format=False)
            ident_to_name[el[0].identifier] = row["crystal"]
    out = j / "crystals"
    out.mkdir(exist_ok=True)
    summary = {}
    for e in expts:
        name = ident_to_name[e.identifier]
        sub = refl.select_on_experiment_identifiers([e.identifier])
        sub.reset_ids()
        el = ExperimentList()
        el.append(e)
        el.as_file(str(out / f"{name}.expt"))
        sub.as_file(str(out / f"{name}.refl"))
        summary[name] = dict(
            identifier=e.identifier,
            space_group=str(e.crystal.get_space_group().info()),
            unit_cell=list(e.crystal.get_unit_cell().parameters()),
            n_reflections=len(sub),
        )
    log = (j / "dials.cosym.log").read_text()
    summary["_cosym_log_tail"] = log[log.find("Best solution"):][:3000]
    (RESULTS / "joint" / coll).mkdir(parents=True, exist_ok=True)
    (RESULTS / "joint" / coll / "crystals.json").write_text(json.dumps(summary, indent=1))
    print(len(expts), "crystals split;", "missing:", sorted(set(ident_to_name.values()) - set(summary)))


if __name__ == "__main__":
    main(sys.argv[1])
