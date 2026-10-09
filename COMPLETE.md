# COMPLETE

Status of the goal (GOAL.md, "Done"), measured from the results files of commit `458cb6bcd63306d74e83e5fa45f492897c64f4a1`
(`python -m cf.report && python -m cf.render`; `bash reproduce.sh` regenerates every number from the
raw data). Each line gives the measurement and the reference it is checked against; details in
`FINDINGS.md` (claims and evidence) and `REPORT.md` (criteria before and after, failures, open items).

## 1. Today's DIALS reproduces the published separation

- 36 cryo crystals, `xia2.multiplex` with every default: 3 clusters
  (12 bovine; 12 porcine; 12 human), 0
  outliers, 3 dimensions, resolution filter 1.26 A;
  the paper: three pure clusters, no outliers, three dimensions, 1.26 A. **Reproduced.**
- Cut to 1.0 and 0.8 degrees per crystal: 3 and
  3 species-pure clusters with 0
  and 1 outlier; the paper: separated down to 0.8 degrees. Reproduced, one
  porcine crystal an outlier at 0.8 degrees.
- The paper's own cosym coordinates shipped with DIALS: cryo labels reproduced exactly; the 215
  room-temperature data sets give 83, 82, 28 and 22 outliers against the published 83, 85, 28 and 19
  (three data sets become outliers under today's defaults, which changed after the paper).

## 2. Deliverable 1: the limit of today's DIALS on both public collections

k90 of today's DIALS (crystals; fraction of the trial), per species pair and data per crystal:

Cryo:

| species pair | 10 deg | 1 deg | 0.5 deg | 0.4 deg | 0.3 deg | 0.2 deg |
|---|---|---|---|---|---|---|
| bovine-human | 5 (29%) | 6 (33%) | 8 (40%) | 12 (50%) | 12 (50%) | > 12 |
| bovine-porcine | 5 (29%) | 8 (40%) | 10 (45%) | > 12 | > 12 | > 12 |
| human-porcine | 5 (29%) | 8 (40%) | 8 (40%) | 12 (50%) | > 12 | > 12 |

Room temperature:

| species pair | 60 deg | 20 deg | 10 deg | 5 deg | 2 deg | 1 deg |
|---|---|---|---|---|---|---|
| bovine-human | 5 (17%) | 5 (17%) | 5 (17%) | 6 (19%) | 12 (32%) | > 25 |
| bovine-porcine | 5 (17%) | 5 (17%) | 6 (19%) | 6 (19%) | 8 (24%) | 8 (24%) |
| human-porcine | 8 (24%) | 6 (19%) | 12 (32%) | > 25 | > 25 | > 25 |

False-group rate in single-population data: cryo 0/189 [0.00, 0.02], room temperature
11/180 [0.03, 0.11]. Minorities of 2-4 crystals are never reported, at any data amount.

## 3. Deliverable 2: the cause, by evidence

Replaying today's DIALS from every trial's own correlation matrix with one stage changed: lowering the
OPTICS minimum-samples floor (5) to 2 makes it report minorities of 2-4 at full data
(cryo k=2-4 at 10 degrees:
0 found with the floor,
73
without it, of 78),
and groups in single-population data (104 of
184 trials against 0); changing the number of
dimensions, the min_samples search or the scaling of the coordinates changes nothing. The minority is
already separated in the coordinates OPTICS receives (cryo, full data: median silhouette of the true
split 0.89-0.98,
every minority crystal's nearest neighbour a minority crystal in 96-100% of draws).
The cause is OPTICS's minimum cluster size: density clustering has no measure of how far a group stands
out of the noise, so its floor is what keeps noise and structure inside one species from being reported,
and it loses every group below it. Room temperature, the same (60 degrees, `FINDINGS.md` F3):
0 of
156 minorities of 2-4 found with the floor,
72 without it, and groups in
5 of 9 single-population trials against 0.

## 4. Deliverable 3: a better procedure, mergeable, with tests

The change: `dials-patch/0001-0015` against DIALS v3.30.1 (`dials.algorithms.correlation.groups`, a
`significant_clusters.find_small_groups` option of `dials.correlation_matrix` and so of `xia2.multiplex`),
tests in `tests/algorithms/correlation/test_groups.py`: 48 passed with DIALS's existing
correlation and regression tests; it passes DIALS's own pre-commit checks (ruff). Run end to end on 72 sampled trials: 0
failures beyond today's DIALS'; the same clusters as the module on the saved arrays in 62 of 72, the
same outcome in 71. 6 configurations were judged, each on data it was not tuned on:

- **v6** (amendment 12, `dials-patch/0001-0007`): after OPTICS, the calibrated detector's groups down to two data sets taken out as clusters (frozen before any room-temperature result).
- **v7** (amendment 17, `dials-patch/0001-0010`): v6, with groups smaller than five returned to the background population.
- **v8** (amendment 20, `dials-patch/0001-0011`): v7, and a group smaller than five when its coherence is precisely measured (standard error at most 0.002), two standard errors above 0.04 x 2/k and twice the rest's structure.
- **v9** (amendment 22, `dials-patch/0001-0013`): v7, and a group smaller than five when its coherence stands 30 standard errors out of the counting noise, two above 0.01 x 2/k and 1.5 times above the rest's structure.
- **v10** (amendment 23, `dials-patch/0001-0014`): v9, with the groups to take out of the OPTICS clusters decided on all of the detector's groups before the minimum size.
- **v11**, the change's defaults (amendment 26, `dials-patch/0001-0015`): v10, and a group smaller than five admitted by v8's condition as well as v9's (small_groups and small_groups.precise).

v6 was judged on room temperature. The others were chosen after its judgement, on both collections; their
held-out judgement is the third collection (Zenodo 15077303, amendments 19-26): none of its trials ran
before all of them were frozen.

| | cryo: found | cryo: false groups | room temperature: found | room temperature: false groups | third collection: found | third collection: false groups | median CPU s per trial (cryo / room temperature) |
|---|---|---|---|---|---|---|---|
| today | 106/720 [0.12, 0.18] | 0/189 [0.00, 0.02] | 416/1548 [0.25, 0.29] | 11/180 [0.03, 0.11] | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | 2.31 / 13.58 |
| v6 | 416/720 [0.54, 0.61] | 0/189 [0.00, 0.02] | 971/1548 [0.60, 0.65] | 33/180 [0.13, 0.25] | 85/110 [0.68, 0.85] | 0/30 [0.00, 0.12] | 2.96 / 46.76 |
| v7 | 174/720 [0.21, 0.27] | 0/189 [0.00, 0.02] | 490/1548 [0.29, 0.34] | 11/180 [0.03, 0.11] | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | 2.49 / 16.47 |
| v8 | 241/720 [0.30, 0.37] | 0/189 [0.00, 0.02] | 718/1548 [0.44, 0.49] | 11/180 [0.03, 0.11] | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | 2.78 / 16.44 |
| v9 | 251/720 [0.31, 0.38] | 0/189 [0.00, 0.02] | 796/1548 [0.49, 0.54] | 11/180 [0.03, 0.11] | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 2.70 / 15.91 |
| v10 | 252/720 [0.32, 0.39] | 0/189 [0.00, 0.02] | 797/1548 [0.49, 0.54] | 11/180 [0.03, 0.11] | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 2.65 / 16.04 |
| v11 | 252/720 [0.32, 0.39] | 0/189 [0.00, 0.02] | 811/1548 [0.50, 0.55] | 11/180 [0.03, 0.11] | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 2.60 / 15.73 |

**v6, held out on room temperature: fails.** It reports a group in 22 single-population trials
where today's DIALS reports none (none the other way round): it raises the false-group rate.
**v7, held out on the third collection: does not raise the false-group rate**
(0 single-population trial(s) with a group where today's DIALS has none, 0 the other way round) and finds
0 more of its 110 minorities than today's DIALS.
**v8, held out on the third collection: does not raise the false-group rate**
(0 single-population trial(s) with a group where today's DIALS has none, 0 the other way round) and finds
0 more of its 110 minorities than today's DIALS.
**v9, held out on the third collection: does not raise the false-group rate**
(0 single-population trial(s) with a group where today's DIALS has none, 0 the other way round) and finds
32 more of its 110 minorities than today's DIALS.
**v10, held out on the third collection: does not raise the false-group rate**
(0 single-population trial(s) with a group where today's DIALS has none, 0 the other way round) and finds
32 more of its 110 minorities than today's DIALS.
**v11, held out on the third collection: does not raise the false-group rate**
(0 single-population trial(s) with a group where today's DIALS has none, 0 the other way round) and finds
32 more of its 110 minorities than today's DIALS. v6 there: 0 trial(s) with a group where today's DIALS has
none, 59 more minorities.

## 5. Deliverable 4

`FINDINGS.md`, `REPORT.md`, this file; `bash reproduce.sh` (one command, staged; `bash reproduce.sh report`
re-derives every number from the results files).

## 6. The smallest detectable group and its cause, for today's DIALS and v11, both collections

k90 (crystals; fraction of the trial), today's DIALS and the 6 judged configurations:

Cryo:

| species pair | wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|---|
| bovine-human | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 2 (14%) | 2 (14%) | 2 (14%) | 2 (14%) |
| bovine-human | 1 deg | 6 (33%) | 2 (14%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| bovine-human | 0.5 deg | 8 (40%) | 3 (20%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| bovine-human | 0.4 deg | 12 (50%) | 4 (25%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) |
| bovine-human | 0.3 deg | 12 (50%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) |
| bovine-human | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| bovine-porcine | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 3 (20%) | 2 (14%) | 2 (14%) | 2 (14%) |
| bovine-porcine | 1 deg | 8 (40%) | 4 (25%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| bovine-porcine | 0.5 deg | 10 (45%) | 5 (29%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) |
| bovine-porcine | 0.4 deg | > 12 | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) |
| bovine-porcine | 0.3 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| bovine-porcine | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| human-porcine | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 3 (20%) | 2 (14%) | 2 (14%) | 2 (14%) |
| human-porcine | 1 deg | 8 (40%) | 2 (14%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| human-porcine | 0.5 deg | 8 (40%) | 4 (25%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| human-porcine | 0.4 deg | 12 (50%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) |
| human-porcine | 0.3 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| human-porcine | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |

Room temperature:

| species pair | wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|---|
| bovine-human | 60 deg | 5 (17%) | 2 (7%) | 5 (17%) | 3 (11%) | 2 (7%) | 2 (7%) | 2 (7%) |
| bovine-human | 20 deg | 5 (17%) | 2 (7%) | 5 (17%) | 4 (14%) | 3 (11%) | 3 (11%) | 3 (11%) |
| bovine-human | 10 deg | 5 (17%) | 4 (14%) | 5 (17%) | 4 (14%) | 5 (17%) | 5 (17%) | 4 (14%) |
| bovine-human | 5 deg | 6 (19%) | 4 (14%) | 6 (19%) | 5 (17%) | 6 (19%) | 6 (19%) | 5 (17%) |
| bovine-human | 2 deg | 12 (32%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) |
| bovine-human | 1 deg | > 25 | 25 (50%) | > 25 | > 25 | > 25 | 25 (50%) | 25 (50%) |
| bovine-porcine | 60 deg | 5 (17%) | 2 (7%) | 5 (17%) | 2 (7%) | 2 (7%) | 2 (7%) | 2 (7%) |
| bovine-porcine | 20 deg | 5 (17%) | 2 (7%) | 5 (17%) | 3 (11%) | 3 (11%) | 3 (11%) | 3 (11%) |
| bovine-porcine | 10 deg | 6 (19%) | 3 (11%) | 6 (19%) | 3 (11%) | 4 (14%) | 4 (14%) | 3 (11%) |
| bovine-porcine | 5 deg | 6 (19%) | 3 (11%) | 6 (19%) | 6 (19%) | 5 (17%) | 5 (17%) | 5 (17%) |
| bovine-porcine | 2 deg | 8 (24%) | 3 (11%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| bovine-porcine | 1 deg | 8 (24%) | 4 (14%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| human-porcine | 60 deg | 8 (24%) | 4 (14%) | 5 (17%) | 5 (17%) | 4 (14%) | 4 (14%) | 4 (14%) |
| human-porcine | 20 deg | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| human-porcine | 10 deg | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) |
| human-porcine | 5 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |
| human-porcine | 2 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |
| human-porcine | 1 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |

The third collection (eight crystals per species; held out):

| species pair | wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|---|
| bovine-human | 20 deg | 5 (38%) | 2 (20%) | 5 (38%) | 5 (38%) | 2 (20%) | 2 (20%) | 2 (20%) |
| bovine-human | 10 deg | 5 (38%) | 2 (20%) | 5 (38%) | 5 (38%) | 3 (27%) | 3 (27%) | 3 (27%) |
| bovine-human | 5 deg | 5 (38%) | 2 (20%) | 5 (38%) | 5 (38%) | 5 (38%) | 5 (38%) | 5 (38%) |
| bovine-human | 2 deg | 5 (38%) | 4 (33%) | 5 (38%) | 5 (38%) | 5 (38%) | 5 (38%) | 5 (38%) |
| bovine-human | 1 deg | 8 (50%) | 8 (50%) | 8 (50%) | 8 (50%) | 8 (50%) | 8 (50%) | 8 (50%) |

Cause, for today's DIALS: section 3. For the judged configurations: every minority they miss, classified by the
measurement that settles it (`cf/limits3.py`; `FINDINGS.md` F6), per collection:

| procedure | collection | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | not classified | procedure |
|---|---|---|---|---|---|---|---|---|
| v6 | cryo | 416 | 252 | 12 | 8 | 26 | 0 | 6 |
| v7 | cryo | 174 | 252 | 12 | 8 | 26 | 0 | 248 |
| v8 | cryo | 241 | 252 | 12 | 8 | 26 | 0 | 181 |
| v9 | cryo | 251 | 252 | 12 | 8 | 26 | 0 | 171 |
| v10 | cryo | 252 | 252 | 12 | 8 | 26 | 0 | 170 |
| v11 | cryo | 252 | 252 | 12 | 8 | 26 | 0 | 170 |
| v6 | rt | 971 | 248 | 306 | 7 | 3 | 0 | 13 |
| v7 | rt | 490 | 252 | 450 | 7 | 3 | 0 | 346 |
| v8 | rt | 718 | 252 | 450 | 7 | 3 | 0 | 118 |
| v9 | rt | 796 | 252 | 450 | 7 | 3 | 0 | 40 |
| v10 | rt | 797 | 252 | 450 | 7 | 3 | 0 | 39 |
| v11 | rt | 811 | 252 | 450 | 7 | 3 | 0 | 25 |

and across the collections (amendment 18: one configuration for every sample, so a group of one species
in either collection's single-population data counts as structure):

| procedure | collection | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | not classified | procedure |
|---|---|---|---|---|---|---|---|---|
| v6 | cryo | 416 | 252 | 31 | 8 | 7 | 0 | 6 |
| v7 | cryo | 174 | 264 | 161 | 8 | 12 | 0 | 101 |
| v8 | cryo | 241 | 264 | 161 | 8 | 12 | 0 | 34 |
| v9 | cryo | 251 | 264 | 161 | 8 | 12 | 0 | 24 |
| v10 | cryo | 252 | 264 | 161 | 8 | 12 | 0 | 23 |
| v11 | cryo | 252 | 264 | 161 | 8 | 12 | 0 | 23 |
| v6 | rt | 971 | 248 | 318 | 7 | 3 | 0 | 1 |
| v7 | rt | 490 | 252 | 469 | 7 | 3 | 0 | 327 |
| v8 | rt | 718 | 252 | 469 | 7 | 3 | 0 | 99 |
| v9 | rt | 796 | 252 | 469 | 7 | 3 | 0 | 21 |
| v10 | rt | 797 | 252 | 469 | 7 | 3 | 0 | 20 |
| v11 | rt | 811 | 252 | 469 | 7 | 3 | 0 | 6 |

## 7. Every remaining limit belongs to the problem itself

**Per cell** (the smallest detectable group, k90, of each species pair and data amount), v11's k90 against
the k90 the data allow with the noise one experiment can know (`FINDINGS.md` F6):

- across the collections (one configuration for every sample, amendment 18): cells where v11's is
  higher: cryo none; room temperature none;
- per collection (a procedure told which kind of sample it has): cryo bovine-porcine at 1 deg (5 against 4); bovine-porcine at 0.5 deg (6 against 5); bovine-human at 1 deg (5 against 2); bovine-human at 0.5 deg (5 against 3); bovine-human at 0.4 deg (6 against 4); human-porcine at 1 deg (5 against 2); human-porcine at 0.5 deg (5 against 4); room
  temperature human-porcine at 60 deg (4 against 3).

How much this depends on the measurement (`FINDINGS.md` F6): under amendment 13's rule, without the
two standard errors' tolerance in the structure test, on amplitude correlations, or on either matrix,
the k90 the data allow (per collection) moves in these cells only, each on minorities that stand out at p
close to the test's level:

- cryo: bovine-porcine 0.4 deg (the rule 10, amendment 13 10, no tolerance 10, amplitudes 8, either matrix 8); bovine-porcine 1 deg (the rule 4, amendment 13 4, no tolerance 2, amplitudes 4, either matrix 4); human-porcine 0.2 deg (the rule > 12, amendment 13 > 12, no tolerance > 12, amplitudes > 12, either matrix 10); human-porcine 0.3 deg (the rule > 12, amendment 13 > 12, no tolerance > 12, amplitudes 10, either matrix 12)
- rt: bovine-porcine 2 deg (the rule 6, amendment 13 6, no tolerance 5, amplitudes 6, either matrix 6); bovine-porcine 5 deg (the rule 5, amendment 13 5, no tolerance 4, amplitudes 5, either matrix 5); human-porcine 10 deg (the rule 12, amendment 13 12, no tolerance 10, amplitudes 12, either matrix 12); human-porcine 60 deg (the rule 3, amendment 13 3, no tolerance 3, amplitudes 4, either matrix 3)

**Per minority**: the minorities v11 misses that no measurement puts in the data, across the
collections (with the p-value of the group that makes each separable, noise from repeats / noise one
experiment can know; whether its detector examined anything; the frozen procedures that do find it):

- cryo: 23 of 720
  - `cryo-mix-bp-k04-d02-w1`: 4 porcine among bovine, 1 degree; separable as a group of 4 holding 4 of them, p 2e-09 / 8e-13; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-bp-k04-d02-w0.5`: 4 porcine among bovine, 0.5 degrees; separable as a group of 4 holding 4 of them, p 6e-05 / 0.0004; detector active; found by v6, v6_1, v5, v4, v3, v2, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-bh-k04-d00-w1`: 4 human among bovine, 1 degree; separable as a group of 4 holding 4 of them, p 1e-32 / 1e-21; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-bh-k04-d01-w1`: 4 human among bovine, 1 degree; separable as a group of 4 holding 4 of them, p 2e-24 / 4e-16; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-bh-k04-d02-w0.5`: 4 human among bovine, 0.5 degrees; separable as a group of 4 holding 4 of them, p 4e-11 / 2e-05; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-ph-k04-d00-w0.5`: 4 human among porcine, 0.5 degrees; separable as a group of 4 holding 4 of them, p 5e-09 / 0.0002; detector active; found by v6, v6_1, v5, v4, v3, v2, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-ph-k04-d02-w0.5`: 4 human among porcine, 0.5 degrees; separable as a group of 4 holding 4 of them, p 1e-11 / 6e-06; detector active; found by v6, v6_1, v5, v4, v3, v2, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-ph-k08-d00-w0.2`: 8 human among porcine, 0.2 degrees; separable as a group of 8 holding 8 of them, p 6e-13 / 3e-10; detector idle; found by no procedure tried
  - `cryo-mix-hb-k04-d00-w1`: 4 bovine among human, 1 degree; separable as a group of 4 holding 4 of them, p 2e-20 / 3e-17; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hb-k04-d01-w0.5`: 4 bovine among human, 0.5 degrees; separable as a group of 4 holding 4 of them, p 2e-13 / 2e-07; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hb-k04-d02-w1`: 4 bovine among human, 1 degree; separable as a group of 4 holding 4 of them, p 2e-26 / 2e-22; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hb-k04-d02-w0.5`: 4 bovine among human, 0.5 degrees; separable as a group of 4 holding 4 of them, p 7e-13 / 3e-07; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hb-k05-d00-w0.2`: 5 bovine among human, 0.2 degrees; separable as a group of 8 holding 5 of them, p 2e-05 / 0.0003; detector idle; found by no procedure tried
  - `cryo-mix-hb-k08-d00-w0.2`: 8 bovine among human, 0.2 degrees; separable as a group of 10 holding 8 of them, p 0.0005 / 0.004; detector idle; found by no procedure tried
  - `cryo-mix-hp-k03-d01-w0.2`: 3 porcine among human, 0.2 degrees; separable as a group of 5 holding 3 of them, p 0.0001 / 0.0003; detector idle; found by no procedure tried
  - `cryo-mix-hp-k04-d00-w1`: 4 porcine among human, 1 degree; separable as a group of 4 holding 4 of them, p 1e-08 / 1e-11; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k04-d00-w0.5`: 4 porcine among human, 0.5 degrees; separable as a group of 4 holding 4 of them, p 8e-07 / 1e-05; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k04-d00-w0.4`: 4 porcine among human, 0.4 degrees; separable as a group of 4 holding 4 of them, p 0.0003 / 0.0008; detector active; found by v6, v6_1, v5, v4, v3, v2, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k04-d01-w1`: 4 porcine among human, 1 degree; separable as a group of 4 holding 4 of them, p 3e-12 / 3e-16; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k04-d01-w0.5`: 4 porcine among human, 0.5 degrees; separable as a group of 4 holding 4 of them, p 6e-06 / 5e-05; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k04-d02-w1`: 4 porcine among human, 1 degree; separable as a group of 4 holding 4 of them, p 1e-05 / 1e-07; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `cryo-mix-hp-k10-d00-w0.2`: 10 porcine among human, 0.2 degrees; separable as a group of 9 holding 5 of them, p 0.0001 / 0.0002; detector idle; found by v3_analytic
  - `cryo-mix-hp-k12-d00-w0.2`: 12 porcine among human, 0.2 degrees; separable as a group of 11 holding 6 of them, p 0.0001 / 0.0002; detector idle; found by no procedure tried
- rt: 6 of 1548
  - `rt-mix-bp-k04-d00-w2`: 4 porcine among bovine, 2 degrees; separable as a group of 4 holding 4 of them, p 5e-37 / 4e-33; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `rt-mix-bp-k04-d05-w2`: 4 porcine among bovine, 2 degrees; separable as a group of 4 holding 4 of them, p 2e-45 / 1e-40; detector active; found by v6, v6_1, v5, v4, v3, v2, v1, v4_ratio_only, v4_coherence_only, v4_bootstrap, v3_ratio_only, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `rt-mix-bp-k05-d02-w5`: 5 porcine among bovine, 5 degrees; separable as a group of 4 holding 4 of them, p 3e-251 / 1e-152; detector active; found by v6, v6_1, v5, v4, v3, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `rt-mix-bp-k05-d02-w2`: 5 porcine among bovine, 2 degrees; separable as a group of 4 holding 4 of them, p 1e-47 / 1e-42; detector active; found by v6, v6_1, v5, v4, v3, v2, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `rt-mix-pb-k03-d00-w1`: 3 bovine among porcine, 1 degree; separable as a group of 4 holding 3 of them, p 3e-08 / 1e-06; detector active; found by v6, v6_1, v5, v4, v3, v4_coherence_only, v4_bootstrap, v3_coherence_only, v3_analytic, v6_m3, v6_m4
  - `rt-mix-hb-k05-d04-w1`: 5 bovine among human, 1 degree; separable as a group of 4 holding 3 of them, p 7e-05 / 0.0008; detector active; found by no procedure tried


**Verdict: holds.** v11 keeps today's false-group rate on the collection it was not tuned on (a weak test of
false groups: v6, which raised the rate on room temperature, reports no group there either, `FINDINGS.md`
F10), and across the collections, in every cell, its detection limit is the limit the data allow one
configuration for every sample. Per collection (a procedure told which kind of sample it has), it gives up the cells listed above:
smaller groups separable from the structure of that collection's own crystals but not from that of the
other's. The individual minorities listed above are separable in the data but missed by v11; they are open
(`REPORT.md`).
