"""Summarise one crystal's DIALS processing (run inside the processing work directory)."""

from __future__ import annotations

import json
import os
import re
import sys

from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory

coll, crystal, species = sys.argv[1:4]
el = ExperimentListFactory.from_json_file("integrated.expt", check_format=False)
refl = flex.reflection_table.from_file("integrated.refl")
strong = flex.reflection_table.from_file("strong.refl")
indexed = flex.reflection_table.from_file("indexed.refl")
e = el[0]
scan = e.scan
n_int = int(refl.get_flags(refl.flags.integrated_prf).count(True))
n_sum = int(refl.get_flags(refl.flags.integrated_sum).count(True))
log = open("dials.refine.log").read()
rmsd = re.findall(r"RMSD_X.*?\n(?:.*\n){2}\s*\|?\s*\d+\s*\|?\s*\d+\s*\|?\s*([\d.]+)\s*\|?\s*([\d.]+)\s*\|?\s*([\d.]+)", log)
out = dict(
    collection=coll,
    crystal=crystal,
    species=species,
    dials_version=os.popen("dials.version").read().splitlines()[0],
    n_images=scan.get_num_images(),
    oscillation=list(scan.get_oscillation()),
    image_range=list(scan.get_image_range()),
    unit_cell=list(e.crystal.get_unit_cell().parameters()),
    space_group=str(e.crystal.get_space_group().info()),
    n_strong=len(strong),
    n_indexed=int(indexed.get_flags(indexed.flags.indexed).count(True)),
    n_integrated_prf=n_int,
    n_integrated_sum=n_sum,
    n_reflections=len(refl),
    wavelength=e.beam.get_wavelength(),
    index_retry=os.path.exists("index_note.txt"),
    rmsd_xy_phi=[float(x) for x in rmsd[-1]] if rmsd else None,
)
print(json.dumps(out, indent=1))
