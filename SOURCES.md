# Sources and licences

The original analysis code, candidate detector, patch additions and reproduction scripts
are copyright 2026 Philipp Bogdan and released under the MIT licence in `LICENSE`.
The Zenodo data, metadata and data-derived result records retain their CC BY 4.0
attribution. Upstream source context in patches remains under its upstream licence.

All six diffraction-data records have the [Creative Commons Attribution 4.0
International licence](https://creativecommons.org/licenses/by/4.0/legalcode).
The public record metadata, file URLs and MD5 checksums are preserved in
`data/zenodo/record_*.json` and `data/manifest.csv`. The raw diffraction files are fetched
by the pipeline, rather than included in this repository.

| Record | Collection and attribution | Licence source |
|---|---|---|
| [10.5281/zenodo.13890874](https://doi.org/10.5281/zenodo.13890874) | Thompson, Amy and Winter, Graeme, Diamond Light Source (2024). Cows, Pigs and People: cubic insulin from three species recorded on I24. 36 cryogenic crystals, 12 per species. | [Record rights](https://zenodo.org/records/13890874) |
| [10.5281/zenodo.15062310](https://doi.org/10.5281/zenodo.15062310) | Diamond Light Source (2025). VMXi in situ bovine insulin. 25 room-temperature crystals. | [Record rights](https://zenodo.org/records/15062310) |
| [10.5281/zenodo.15062327](https://doi.org/10.5281/zenodo.15062327) | Diamond Light Source (2025). VMXi in situ porcine insulin. 25 room-temperature crystals. | [Record rights](https://zenodo.org/records/15062327) |
| [10.5281/zenodo.15062343](https://doi.org/10.5281/zenodo.15062343) | Diamond Light Source (2025). VMXi in situ human insulin. 25 room-temperature crystals. | [Record rights](https://zenodo.org/records/15062343) |
| [10.5281/zenodo.15062350](https://doi.org/10.5281/zenodo.15062350) | Diamond Light Source (2025). VMXi in situ unclassified insulin. 19 crystals included in the all-crystal fidelity check. | [Record rights](https://zenodo.org/records/15062350) |
| [10.5281/zenodo.15077303](https://doi.org/10.5281/zenodo.15077303) | Diamond Light Source (2025). VMXi bovine and human insulin. The held-out third collection, eight crystals per species. Labels follow the supporting information and the verification recorded in PREREGISTRATION.md amendment 25. | [Record rights](https://zenodo.org/records/15077303) |

Thompson, Beilsten-Edmands, Tam, Sanchez-Weatherby et al. (2025), *Enhanced
intensity-based clustering of isomorphous multi-crystal data sets in the presence of
subtle variations*, Acta Crystallographica D 81, 278-290,
[doi:10.1107/S2059798325004589](https://doi.org/10.1107/S2059798325004589), is the published
comparison. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12128884/) and supporting
information are available from the publisher. The downloaded full-text XML is not
vendored here.

| Software | Pinned source | Licence |
|---|---|---|
| DIALS | [v3.30.1](https://github.com/dials/dials/tree/v3.30.1), base source commit `17cfc5af87c711969b11e2848b1bbea45118fc27` | [BSD-3-Clause](https://github.com/dials/dials/blob/v3.30.1/LICENSE); notice and terms in `licences/DIALS-BSD-3-Clause.txt` |
| xia2 | [v3.30.0](https://github.com/xia2/xia2/tree/v3.30.0) | [BSD-3-Clause](https://github.com/xia2/xia2/blob/v3.30.0/LICENSE); notice and terms in `licences/XIA2-BSD-3-Clause.txt` |

Neither upstream source tree is vendored. The environment script retrieves pinned DIALS
source for its reference coordinates and regression tests; conda-forge supplies both
installed packages and their dependencies. `cf/frozen/groups_*.py` are the candidate
detector versions introduced by this research, retained to reproduce the frozen
configurations, rather than copies of a pre-existing third-party detector.

`dials-patch/manifest.json` maps the original patch MD5 checksums in the historic test
log to the public files' SHA-256 checksums. Only the author email headers have changed;
the patch diff bodies match the original files byte for byte.

`results/matrix_names/*.jsonl.gz` are lossless projections of the saved DIALS correlation
matrix records, containing only trial identifiers and crystal names in matrix order.
They retain the inputs for FINDINGS.md F9 without distributing regenerable matrices.
Local paths and operational details in the public copy have been neutralised; these are
the only editorial changes to the research documents. Numeric measurements are unchanged.
