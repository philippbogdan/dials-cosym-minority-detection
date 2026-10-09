"""Fidelity against DIALS's own reference data: the cosym coordinates of the paper's data
sets shipped in DIALS's tests (tests/algorithms/correlation/test_coords: data_2 = cryo 36,
data_4 = room temperature 215) with their expected labels (the paper's clusters). Runs
today's OPTICS optimisation (CorrelationMatrix.optimise_clustering, DIALS 3.30.1) with the
initial min_samples DIALS's own test uses and with the heuristic of the code."""

from __future__ import annotations

import json
import warnings
from collections import Counter

import numpy as np
from sklearn.metrics import adjusted_rand_score, v_measure_score

from cf.paths import RESULTS, ROOT


def main():
    from dials.algorithms.correlation.analysis import CorrelationMatrix

    base = ROOT / "external" / "dials" / "tests" / "algorithms" / "correlation" / "test_coords"
    cases = {"data_2 (cryo, 36)": (2, 6), "data_4 (room temperature, 215)": (4, 27)}
    out = {}
    for name, (i, ms_test) in cases.items():
        x = np.loadtxt(base / f"data_{i}.txt")
        expected = np.loadtxt(base / f"labels_{i}.txt").astype(int)
        dims = x.shape[1]
        ms_heuristic = max(5, int(len(x) / dims * 0.5))
        rec = dict(n=len(x), dims=dims, expected=dict(zip(*map(lambda a: a.tolist(), np.unique(expected, return_counts=True)))))
        for label, ms in (("min_samples as in DIALS's test", ms_test), ("code heuristic", ms_heuristic)):
            runs = []
            for rep in range(5):  # optuna's TPE stage is unseeded: repeat
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    _, _, labels, _ = CorrelationMatrix.optimise_clustering(x, list(range(len(x))), initial_min_samples=ms)
                labels = np.asarray(labels)
                pairs = Counter(zip(expected.tolist(), labels.tolist()))
                runs.append(dict(counts=dict(zip(*map(lambda a: a.tolist(), np.unique(labels, return_counts=True)))),
                                 contingency={f"{a}->{b}": c for (a, b), c in sorted(pairs.items())},
                                 identical=bool(np.array_equal(labels, expected)),
                                 v_measure=float(v_measure_score(expected, labels)),
                                 ari=float(adjusted_rand_score(expected, labels))))
            rec[label] = dict(initial_min_samples=ms, runs=runs)
        out[name] = rec
    (RESULTS / "fidelity").mkdir(parents=True, exist_ok=True)
    (RESULTS / "fidelity" / "dials_test_coords.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:3000])


if __name__ == "__main__":
    main()
