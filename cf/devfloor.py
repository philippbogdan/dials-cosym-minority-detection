"""After the room-temperature judgement of v6 (amendment 17): v6 with a minimum group size m, derived
exactly from v6's own detector groups: groups smaller than m are merged into the largest population
before the take-out, as v7-v9 do (v6_m5 = v7); the others are taken out as the change does. Writes results/trials/<design>/v6_m<m>/shard_000.jsonl.
Usage: python -m cf.devfloor DESIGN M [M ...]"""

from cf.result_io import open_result

import json
import sys

from cf.analyze import load_design, load_results
from cf.newproc import carve
from cf.paths import RESULTS


def main(design, ms):
    trials = load_design(design)
    today = load_results(design, "today")
    v6 = load_results(design, "v6")
    for m in ms:
        od = RESULTS / "trials" / design / f"v6_m{m}"
        od.mkdir(parents=True, exist_ok=True)
        with open_result(od / "shard_000.jsonl", "w") as fh:
            for tid in sorted(v6):
                r = v6[tid]
                rec = dict(id=tid, procedure=f"v6_m{m}")
                if "error" in r or tid not in today or "error" in today[tid]:
                    rec["error"] = r.get("error") or "today's DIALS failed"
                    fh.write(json.dumps(rec) + "\n")
                    continue
                t = trials[tid]
                names = t["crystals"] if "crystals" in t else t["majority"] + t["minority"]
                groups = sorted(r["detector_groups"], key=len, reverse=True)
                kept = [g for g in groups[1:] if len(g) >= m]
                small = [x for g in groups[1:] if len(g) < m for x in g]
                rest = (groups[0] if groups else []) + small
                rec.update(clusters=carve(today[tid].get("clusters", []), [rest] + kept if kept else [rest], names, "half"),
                           detector_groups=[rest] + kept, detector_cpu_s=r.get("detector_cpu_s"),
                           noise_cpu_s=r.get("noise_cpu_s"), optics_cpu_s=r.get("optics_cpu_s"))
                fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2:]])
