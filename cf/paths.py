from __future__ import annotations

from cf.result_io import open_result

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
import os

# Processed data and trial arrays are regenerated in .venv-data.
# CF_DATA can select a different data directory.
DATA = Path(os.environ["CF_DATA"]) if os.environ.get("CF_DATA") else ROOT / ".venv-data"
RESULTS = ROOT / "results"

# A design name with the suffix "_amp" (cryo_amp, rt_amp) means the trials of the base design
# with the pairwise matrix computed from amplitudes instead of intensities (cf.fastcc.trial_arrays);
# everything derived from the matrix (cmat, noise replicates, noise calibration, limits) is kept
# under the suffixed name.
AMP = "_amp"


# The suffix "_xcal" (cryo_xcal, rt_xcal): the limits of the design computed with a noise calibration
# learned on the other collection (one factor, its mean over wedges and species), as a procedure with
# a frozen calibration could use it (cf.limits3.calibration).
XCAL = "_xcal"


def matrix_design(design: str) -> str:
    """The design whose matrices and noise replicates a design uses (without "_xcal")."""
    return design[: -len(XCAL)] if design.endswith(XCAL) else design


def base_design(design: str) -> str:
    design = matrix_design(design)
    return design[: -len(AMP)] if design.endswith(AMP) else design


def representation(design: str) -> str:
    return "amplitude" if matrix_design(design).endswith(AMP) else "intensity"


def manifest(collection: str | None = None) -> list[dict]:
    with open_result(ROOT / "data" / "manifest.csv") as fh:
        rows = list(csv.DictReader(fh))
    return [r for r in rows if collection is None or r["collection"] == collection]


def species_of(collection: str) -> dict[str, str]:
    return {r["crystal"]: r["species"] for r in manifest(collection)}
