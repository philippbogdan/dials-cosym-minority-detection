# Commands and inputs

`bash reproduce.sh report` reads the committed records, writes `results/report/numbers.json`
and renders `COMPLETE.md` and `FINDINGS.md`. It needs Python 3 with NumPy, SciPy
and Jinja2. If these are unavailable, the script creates `.venv-report` and installs
`requirements-report.txt`. No DIALS installation, diffraction images or external source
tree is needed for this step. `REPORT_PYTHON` can select a Python executable.

The JSON Lines files in `results/trials`, `results/limits`, `results/cause` and
`results/dev` are losslessly compressed as `.jsonl.gz`. `cf.result_io` reads them directly.
Fresh uncompressed records take precedence over compressed records of the same name.
For an independent analysis, use `gzip.open(path, "rt")` and decode each line with
`json.loads`. Some limit records contain Python JSON's `NaN` and `Infinity` values,
retained from the original measurement.

`PROVENANCE.json` keeps the original run's final and verdict commit identifiers. The
documents continue to identify that evidence, rather than the new public repository's
history. Result `code_commit` fields refer to the original research history. The
preparation has changed paths and operational wording, not measured values.

To regenerate the experiment from diffraction images, use a separate clone and retain
the committed records in the first clone for comparison:

```bash
python3 -m pip install -r requirements-report.txt
bash reproduce.sh env
mv results results-recorded
mkdir results
cp -R results-recorded/design results/design
mkdir -p results/dev results/fidelity
cp results-recorded/dev/configurations.txt results/dev/
cp results-recorded/fidelity/cryo_all.txt results/fidelity/
PARALLEL=1 bash reproduce.sh process joint fidelity pilot design trials analysis heldout limits dev e2e tests report
```

The `env` stage installs DIALS 3.30.1 and xia2 3.30.0 from conda-forge and downloads the
DIALS v3.30.1 source into ignored `external/dials`. Patches are applied to an environment
overlay and to a separate test source tree. `PATCH_MAX`, `PATCHES`, `OVERLAY` and
`FROZEN_PHIL` select the historical configurations, as recorded in `dials-patch/README.md`.
The default patch set is 0001-0015, the proposed v11 procedure.

The `process` stage downloads the Zenodo files listed in `data/manifest.csv` and checks
each MD5 checksum. About 71 GB of raw files are fetched across all collections, one
crystal at a time. Integrated data, consistently indexed data, scaled trial arrays and
noise replicates are regenerated under `.venv-data`. Plan additional storage and memory
for these, particularly the 94-crystal room-temperature indexing step. `PARALLEL`
controls concurrency; start with 1 and increase only with sufficient memory. Original
CPU times are measurements of the recorded run, not predictions for another computer.
Unseeded stages of today's DIALS can change clustering between full reruns.

`bash reproduce.sh tests` runs the new and existing DIALS correlation tests, their
regression tests, and this repository's tests. It requires the `env` stage and fetches
regression data through `dials-data`. New test output goes into ignored
`.venv-data/verification`, preserving the original test records in `results/tests`.
Tests that require saved trial arrays skip if those arrays have not been regenerated.

Inputs for plotting detection against data per crystal, false groups and k90 are in
`results/design/{cryo,rt,rtbh}.jsonl` and the corresponding
`results/trials/<collection>/{today,v6,v7,v8,v9,v10,v11}/shard_*.jsonl.gz`.
The minimum-cluster-size replay and coordinate separation measurements are in
`results/cause/{cryo,rt}/shard_*.jsonl.gz`. The report's aggregate detection cells retain
the draw counts needed for exact confidence intervals. No figures are included.
