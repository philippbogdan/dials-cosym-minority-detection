"""Development check (F5): Can one experiment know its counting noise? Jackknife (leave one reflection out) standard error of
every pair's weighted CC, against the scatter of the same pair's residual across disjoint repeat
segments (as cf.noisecal), cryo."""

from cf.result_io import open_result

import json
import numpy as np
from cf.analyze import load_design
from cf.fastcc import trial_arrays
from cf.noisecal import residual
from cf.bootstrap import analytic_replicates, neff_matrix
from cf.paths import DATA

def jack_se(y, v, m):
    n = y.shape[0]
    se = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            mm = m[i] & m[j]
            c = mm.sum()
            if c < 5:
                continue
            a, b = y[i, mm], y[j, mm]
            w = 1.0 / (v[i, mm] + v[j, mm])
            def wcc(sel):
                ww = w[sel] / w[sel].sum()
                da, db = a[sel] - (ww * a[sel]).sum(), b[sel] - (ww * b[sel]).sum()
                return (ww * da * db).sum() / np.sqrt((ww * da * da).sum() * (ww * db * db).sum())
            idx = np.arange(c)
            rs = np.array([wcc(idx != h) for h in range(c)])
            rs = rs[np.isfinite(rs)]
            if len(rs) < 4:
                continue
            se[i, j] = se[j, i] = np.sqrt((len(rs) - 1) / len(rs) * ((rs - rs.mean()) ** 2).sum())
    return se

trials = load_design("cryo")
out = {}
for w in (0.2, 0.3, 0.4, 0.5, 1):
    for sp in ("bovine", "porcine", "human"):
        es, sj, sa, names0 = [], [], [], None
        for seg in range(4):
            tid = f"cryo-single-{sp[0]}-all-w{w:g}" + (f"-s{seg}" if seg else "")
            f = DATA / "trials" / "cryo" / f"{tid}.npz"
            if not f.exists():
                continue
            a = trial_arrays(np.load(f))
            names = a["names"]
            names0 = names0 or names
            if sorted(names) != sorted(names0):
                continue
            idx = [names.index(x) for x in names0]
            r = a["rij"].copy(); np.fill_diagonal(r, 0)
            e = residual(r, a["wij"])
            ana = np.array([residual(r + x, a["wij"]) for x in analytic_replicates(r, neff_matrix(a["y"], a["v"], a["m"]), B=50).astype(float)])
            sub = np.ix_(idx, idx)
            es.append(e[sub]); sj.append(jack_se(a["y"], a["v"], a["m"])[sub]); sa.append(np.nanstd(ana, axis=0)[sub])
        es, sj, sa = map(np.array, (es, sj, sa))
        iu = np.triu_indices(es.shape[1], 1)
        seg_sd = np.nanstd(es[:, iu[0], iu[1]], axis=0, ddof=1)
        c4 = {3: 0.886, 4: 0.921}[len(es)]
        sa = np.where(sa > 0, sa, np.nan)
        rj = np.sqrt(np.nanmean(sj[:, iu[0], iu[1]] ** 2, axis=0)); ra = np.sqrt(np.nanmean(sa[:, iu[0], iu[1]] ** 2, axis=0))
        ok = np.isfinite(seg_sd) & np.isfinite(rj) & np.isfinite(ra)
        out[f"{w:g}|{sp}"] = dict(obs_over_jack=float(np.nanmean(seg_sd[ok]) / c4 / np.nanmean(rj[ok])),
                                  obs_over_analytic=float(np.nanmean(seg_sd[ok]) / c4 / np.nanmean(ra[ok])), n=int(ok.sum()))
        print(w, sp, {k: round(v, 2) for k, v in out[f"{w:g}|{sp}"].items()}, flush=True)
from cf.paths import RESULTS
json.dump(out, open_result(RESULTS / "dev" / "cryo" / "jackknife_cal.json", "w"), indent=1)
