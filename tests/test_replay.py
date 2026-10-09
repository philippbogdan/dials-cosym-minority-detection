"""cf.replay reproduces dials.correlation_matrix's dimensions and clusters from saved rij/wij."""

from cf.result_io import open_result, glob_results


import glob
import json
import os

import numpy as np
import pytest

from cf.paths import DATA, RESULTS
from cf.replay import replay


def _cases():
    out = []
    for f in sorted(glob_results(str(RESULTS / "trials" / "*" / "today" / "shard_*.jsonl"))):
        design = f.split(os.sep)[-3]
        for line in open_result(f):
            r = json.loads(line)
            npz = DATA / "trials" / design / f"{r['id']}.npz"
            if "error" not in r and npz.exists() and "eps1" in np.load(npz).files:
                out.append((r, npz))
    return out[:8]


CASES = _cases()


def partition(labels_by_name):
    groups = {}
    for name, lab in labels_by_name.items():
        if lab >= 0:
            groups.setdefault(lab, set()).add(name)
    return sorted(sorted(g) for g in groups.values())


@pytest.mark.skipif(not CASES, reason="no saved trials")
@pytest.mark.parametrize("rec,npz", CASES)
def test_replay_matches_full_run(rec, npz):
    d = np.load(npz)
    out = replay(d["rij"], d["wij"])
    names = list(d["names"])
    assert out["dims"] == rec["dims"]
    got = partition({n: int(l) for n, l in zip(names, out["labels"])})
    want = partition(rec["labels"])
    assert got == want
