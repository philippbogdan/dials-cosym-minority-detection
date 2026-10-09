"""Build data/manifest.csv: one row per crystal with its species label and Zenodo files.

Labels come only from the Zenodo records (file prefixes for 13890874, record titles for
the room-temperature records); nothing here is inferred from the diffraction data. For the third
collection (rtbh, Zenodo 15077303, PREREGISTRATION.md amendment 19) the paper's supporting
information (Figure S1: data sets 0-7 bovine, 8-15 human); which image numbers they are is settled by
the label verification (amendment 25): the first batch of image numbers human, the second bovine.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
Z = HERE / "data" / "zenodo"

CRYO = "13890874"
RT = {"15062310": "bovine", "15062327": "porcine", "15062343": "human", "15062350": "unclassified"}
PREFIX = {"CIX": "bovine", "PIX": "porcine", "X": "human"}
RTBH = "15077303"


def files_of(record):
    d = json.loads((Z / f"record_{record}.json").read_text())
    return {f["key"]: (f["size"], f["checksum"], f["links"]["self"]) for f in d["files"]}


def main():
    rows = []
    for key, (size, md5, url) in sorted(files_of(CRYO).items()):
        m = re.fullmatch(r"(CIX|PIX|X)(\d+)_1\.tar", key)
        assert m, key
        rows.append(
            dict(collection="cryo", crystal=key[:-4], species=PREFIX[m.group(1)],
                 record=CRYO, files=key, md5=md5, bytes=size, urls=url)
        )
    for record, species in RT.items():
        fs = files_of(record)
        stems = sorted({re.match(r"(image_\d+)", k).group(1) for k in fs})
        for stem in stems:
            keys = sorted(k for k in fs if k.startswith(stem + ".") or k.startswith(stem + "_"))
            assert len(keys) == 4, (stem, keys)
            rows.append(
                dict(collection="rt", crystal=stem, species=species, record=record,
                     files=";".join(keys), md5=";".join(fs[k][1] for k in keys),
                     bytes=sum(fs[k][0] for k in keys), urls=";".join(fs[k][2] for k in keys))
            )
    fs = files_of(RTBH)
    stems = sorted({re.match(r"(image_\d+)", k).group(1) for k in fs})
    assert len(stems) == 16, stems
    for i, stem in enumerate(stems):  # first batch human, second bovine (amendment 25)
        keys = sorted(k for k in fs if k.startswith(stem + ".") or k.startswith(stem + "_"))
        assert len(keys) == 4, (stem, keys)
        rows.append(
            dict(collection="rtbh", crystal=stem, species="human" if i < 8 else "bovine", record=RTBH,
                 files=";".join(keys), md5=";".join(fs[k][1] for k in keys),
                 bytes=sum(fs[k][0] for k in keys), urls=";".join(fs[k][2] for k in keys))
        )
    with open(HERE / "data" / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(Counter((r["collection"], r["species"]) for r in rows))
    print(sum(int(r["bytes"]) for r in rows) / 1e9, "GB")


if __name__ == "__main__":
    main()
