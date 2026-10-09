# Work log (state of the work, for resuming)

## Layout

- `data/manifest.csv`: the 130 crystals (Zenodo records, species labels from the records).
- `pipeline/`: shell steps (processing, joint cosym, multiplex, runners,
  `analysis_stage.sh DESIGN P STAGE...` = everything after today's DIALS for one collection).
- `cf/`: the analysis package. Trial design (`design`), today's DIALS runner (`run_trials`,
  `today`), noise replicates (`bootstrap`: reflection bootstrap, parametric, analytic Fisher-z),
  noise calibration against repeat segments (`noisecal`), the frozen procedures from saved arrays
  (`newproc`, frozen copies of earlier module versions in `cf/frozen/`), limits (`limits`:
  pre-registered oracle/decoy; `limits3`: refined, amendments 4-14), cause ablations (`cause`,
  `replay`), prototypes (`proc5`, `proc6`, `devrun`, `devalt` = amplitude correlations,
  `devunion`, `devsinglep`, `devjack`), report and rendering (`report`, `render`).
- Design names with the suffix `_amp` (cryo_amp, rt_amp): the same trials with the pairwise matrix
  computed from amplitudes (cf.paths.base_design / representation; cf.fastcc.trial_arrays).
- `external/dials` (git branch `calibrated-groups` on DIALS v3.30.1, not committed here; the
  change is exported to `dials-patch/0001-0007`): the proposed DIALS change with its tests.
- `.venv-dials/` and `.venv-data/`: regenerated environments, processed data and trial arrays.
- `results/`: everything measured. Derived caches not in git: `results/cmat`, `results/boot*`.

## Procedures (PREREGISTRATION.md)

- today: DIALS 3.30.1 defaults (OPTICS significant clusters).
- v1 (amendment 3) ... v5 (amendment 9): see PREREGISTRATION.md; frozen tags frozen-v1 ... frozen-v5.
- **v6 (amendment 12, tag frozen-v6) = dials-patch/0001-0007, the proposed procedure**: v5's detector
  with the analytic noise and with the reflection bootstrap, groups combined. Its rt numbers are the
  held-out judgement.
- Development after v6 on cryo only, not adopted: amplitude correlations (results/dev/cryo/alt_*,
  amendment 14).

## Limits (cf/limits3.py)

Per mixture trial: counting noise (repeat-calibrated analytic noise; look-elsewhere over subsets and
sizes; alpha 0.001; standardised block sum or inverse-variance sum; the minority or any sub-group of
at least max(2, k/2) members, amendment 14), structure inside one species (single-population
comparators of the same size at the same or larger wedge, dominance in S and ratio within 2 SE; and
the trial's own majority), noise known only from repeats (model noise at face value), else
procedure. A minority found by an FG-safe frozen procedure is never a data limit.
Robustness (reported): amendment 13's rule; the amplitude matrix (`_amp`).

## State (2026-10-09)

- **The change as proposed: v11** (dials-patch/0001-0015, frozen.phil, tag frozen-v11; amendment 26). Tests
  48 passed (DIALS correlation + regression) and 24 (this repository); ruff clean; end to end on 72 sampled
  trials no failure beyond today's DIALS, the same detector groups as from the arrays in 72, the same
  outcome in 71.
- In-sample: cryo 252/720 (today 106), 0/189 false groups; rt 811/1548 (today 416), 11/180 (today's, the
  same trials). Across the collections (amendment 18, with amendment 24's comparators) its k90 is the
  data's in every cell of both collections; per minority 23 cryo + 6 rt separable but missed (open).
- **Held out, the third collection** (amendments 19-26; labels per amendment 25, verified 16/16): today
  26/110, 0/30; v11 58/110, 0/30, none added -> passes. v6 85/110, 0/30 (so the test of false groups is
  weak there: these crystals hold no sub-group v6 reports; 0/30 allows up to 12%).
- COMPLETE.md: "Verdict: holds" (held-out judgement passed; every cell at the data's limit across the
  collections; per-collection cells and individual minorities listed as open).
- Judged configurations: v6 (12), v7 (17), v8 (20), v9 (22), v10 (23), v11 (26); report.CONFIGURATIONS.
- Today's mishaps (all corrected, see REPORT "Defects"): newproc name collision (first v11 runs), an e2e
  run reading a phil file changed mid-run (rerun pinned; mixed records removed), a cause replay mixing two
  runs (rerun), a zsh word-splitting slip that deleted the e2e results (rerun from pinned
  commands; local copies kept).

## Next

- Nothing required for the goal. Possible further work (open in REPORT): the individual minorities, a
  third collection with sub-group structure to test false groups more strongly.
