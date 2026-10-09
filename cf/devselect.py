"""Selection check of the limit test (amendment 14): in single-population trials, random
pseudo-minorities (k crystals; the rest as majority) go through the same search as real
minorities (the group itself, its qualifying sub-groups and mixed extensions, both statistics);
the fraction whose best candidate stands out at the limit level should be at most that level."""

from cf.result_io import open_result

import json, math, sys
import numpy as np
from cf.analyze import load_design
from cf.bootstrap import load as bload
from cf.cmatrices import load
from cf.limits3 import (Trial, calibration, mixed_candidates, subgroup_candidates, LIMIT_ALPHA)
from cf.paths import RESULTS


def main(design, n_draws=10, ks=(2, 3, 4, 5, 6), seed=0):
    rng = np.random.default_rng(seed)
    trials = load_design(design)
    mats = load(design)
    cal = calibration(design)
    ids = sorted(t for t in trials if trials[t]["kind"] == "single" and t in mats)
    boots = bload(design, "analytic", ids=set(ids))
    out = []
    for tid in ids:
        m = mats[tid]
        T = Trial(m["rij"].astype(float), m["wij"].astype(float), boots[tid].astype(float),
                  cal.get(float(trials[tid]["w"]), 1.0))
        n = T.z.shape[0]
        for k in ks:
            if n - k < 4 or k > n // 2:
                continue
            for _ in range(n_draws):
                M = np.sort(rng.choice(n, k, replace=False))
                A = np.setdiff1d(np.arange(n), M)
                pure = subgroup_candidates(T, M, int(math.ceil(max(2, k / 2))))
                mixed = [G for G, _ in mixed_candidates(T, [M] + pure, A)]
                ps = {}
                for name, groups in (("whole", [M]), ("pure", [M] + pure), ("mixed", [M] + pure + mixed)):
                    ps[name] = min(min(T.p_noise_any_size(G), T.p_noise_any_size_iv(G)) for G in groups)
                out.append(dict(id=tid, w=trials[tid]["w"], k=k, **{f"p_{a}": float(b) for a, b in ps.items()}))
    (RESULTS / "dev" / base(design)).mkdir(parents=True, exist_ok=True)
    json.dump(out, open_result(RESULTS / "dev" / base(design) / f"selection_{design}.json", "w"))
    for name in ("whole", "pure", "mixed"):
        p = np.array([o[f"p_{name}"] for o in out], float)
        print(design, name, "draws", len(p), "p<%g: %d (%.2f%%)" % (LIMIT_ALPHA, (p < LIMIT_ALPHA).sum(), 100 * (p < LIMIT_ALPHA).mean()),
              "p<0.01: %.1f%%" % (100 * (p < 0.01).mean()))


def base(design):
    from cf.paths import base_design
    return base_design(design)


if __name__ == "__main__":
    main(sys.argv[1])
