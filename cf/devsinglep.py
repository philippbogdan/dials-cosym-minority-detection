"""Development check (amendment 14): p-value (look-elsewhere over subsets and sizes, calibrated noise) of the strongest group in
each single-population trial: uniform if the limit test's noise model is right and the trial
holds no structure; piled at 1 if the model overstates the noise."""

from cf.result_io import open_result

import json, sys
import numpy as np
from cf.analyze import load_design
from cf.bootstrap import load as bload
from cf.cmatrices import load
from cf.limits2 import gumbel_p
from cf.limits3 import Trial, calibration, max_stouffer, stouffer
from cf.paths import RESULTS
from cf.scan import grow_all

design = sys.argv[1]
trials = load_design(design)
mats = load(design)
cal = calibration(design)
ids = sorted(t for t in trials if trials[t]["kind"] == "single" and t in mats)
boots = bload(design, "analytic", ids=set(ids))
out = []
for tid in ids:
    t = trials[tid]; m = mats[tid]
    for model_noise in (False, True):
        factor = cal.get(float(t["w"]), 1.0) if not model_noise else 1.0
        T = Trial(m["rij"].astype(float), m["wij"].astype(float), boots[tid].astype(float), factor)
        n = T.z.shape[0]
        sims = T.null_all_sizes()
        mu, sd = sims.mean(0), sims.std(0) + 1e-12
        best = -np.inf
        xx = np.where(T.valid, T.z, 0.0)
        order, within, _ = grow_all(xx, n // 2)
        for k in range(2, n // 2 + 1):
            sums = within[:, k] * (k * (k - 1) / 2)
            zk = max(stouffer(T.z, T.valid, order[b, :k]) for b in np.argsort(sums)[::-1][:5])
            best = max(best, (zk - mu[k - 2]) / sd[k - 2])
        p = gumbel_p(best, ((sims - mu) / sd).max(1))
        out.append(dict(id=tid, w=t["w"], model_noise=model_noise, p=p))
(RESULTS / "dev" / "cryo").mkdir(parents=True, exist_ok=True)
json.dump(out, open_result(RESULTS / "dev" / "cryo" / f"single_pany_{design}.json", "w"))
for mn in (False, True):
    for w in sorted({o["w"] for o in out}):
        ps = np.array([o["p"] for o in out if o["w"] == w and o["model_noise"] == mn])
        print(design, "model noise" if mn else "calibrated ", f"w={w:g}", "n", len(ps), "p<0.001", int((ps < 0.001).sum()), "p<0.01", int((ps < 0.01).sum()), "p<0.1", int((ps < 0.1).sum()), "median", round(float(np.median(ps)), 3))
