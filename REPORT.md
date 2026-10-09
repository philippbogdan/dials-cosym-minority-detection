# Report

Numbers from the results files of commit `458cb6bcd63306d74e83e5fa45f492897c64f4a1` (`python -m cf.report && python -m
cf.render`; `bash reproduce.sh` regenerates everything). "Today" is DIALS 3.30.1 with every default
(`dials.correlation_matrix`, the clustering `xia2.multiplex` runs). 6 configurations of the proposed
change (`dials-patch/`) are judged, each on data it was not tuned on:

- **v6** (amendment 12, `dials-patch/0001-0007`): after OPTICS, the calibrated detector's groups down to two data sets taken out as clusters (frozen before any room-temperature result). Cryo numbers in-sample; **room temperature is its held-out judgement, and it fails it: it
  raises the false-group rate there** (`FINDINGS.md` F4, F8).
- **v7** (amendment 17, `dials-patch/0001-0010`): v6, with groups smaller than five returned to the background population.
- **v8** (amendment 20, `dials-patch/0001-0011`): v7, and a group smaller than five when its coherence is precisely measured (standard error at most 0.002), two standard errors above 0.04 x 2/k and twice the rest's structure.
- **v9** (amendment 22, `dials-patch/0001-0013`): v7, and a group smaller than five when its coherence stands 30 standard errors out of the counting noise, two above 0.01 x 2/k and 1.5 times above the rest's structure.
- **v10** (amendment 23, `dials-patch/0001-0014`): v9, with the groups to take out of the OPTICS clusters decided on all of the detector's groups before the minimum size.
- **v11**: the change as proposed, its defaults (amendment 26, `dials-patch/0001-0015`): v10, and a group smaller than five admitted by v8's condition as well as v9's (small_groups and small_groups.precise).

Cryo and room-temperature numbers of v7, v8, v9, v10, v11 are in-sample; **their held-out judgement is the third
collection** (Zenodo 15077303, amendments 19-26), none of whose trials ran before all were frozen.

Claims and their evidence: `FINDINGS.md`.

## Each criterion, before and after

### Fidelity (reference: the paper)

| setting | paper | today's DIALS here |
|---|---|---|
| 36 cryo crystals, full data (`xia2.multiplex`) | 3 species-pure clusters, no outliers, 3 dimensions, 1.26 A | 3 clusters, species-pure: True, 0 outliers, 3 dimensions, 1.26 A |
| cryo cut to 1 degree | separated (down to 0.8 degrees) | 3 clusters, species-pure: True, 0 outliers |
| cryo cut to 0.8 degrees | separated; below it data sets become outliers | 3 clusters, species-pure: True, 1 outlier |
| the paper's 215 room-temperature data sets (DIALS's reference coordinates) | 83, 85, 28 and 19 outliers | 83, 82, 28 and 22 outliers: 3 data sets of one cluster become outliers, every other label as published |

Difference explained: on the room-temperature reference set three data sets of one cluster are outliers
under today's defaults (the defaults changed after the paper); at 0.8 degrees one porcine crystal is an
outlier where the paper reports none until below 0.8 degrees.

### Detection (reference: species labels of the Zenodo records)

k90 (crystals; fraction of the trial), before (today) and after (v6, v7, v8, v9, v10, v11):

| collection | species pair | wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|
| cryo | bovine-human | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 2 (14%) | 2 (14%) | 2 (14%) | 2 (14%) |
| cryo | bovine-human | 1 deg | 6 (33%) | 2 (14%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| cryo | bovine-human | 0.5 deg | 8 (40%) | 3 (20%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| cryo | bovine-human | 0.4 deg | 12 (50%) | 4 (25%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) |
| cryo | bovine-human | 0.3 deg | 12 (50%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) |
| cryo | bovine-human | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| cryo | bovine-porcine | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 3 (20%) | 2 (14%) | 2 (14%) | 2 (14%) |
| cryo | bovine-porcine | 1 deg | 8 (40%) | 4 (25%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| cryo | bovine-porcine | 0.5 deg | 10 (45%) | 5 (29%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) | 6 (33%) |
| cryo | bovine-porcine | 0.4 deg | > 12 | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) | 10 (45%) |
| cryo | bovine-porcine | 0.3 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| cryo | bovine-porcine | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| cryo | human-porcine | 10 deg | 5 (29%) | 2 (14%) | 5 (29%) | 3 (20%) | 2 (14%) | 2 (14%) | 2 (14%) |
| cryo | human-porcine | 1 deg | 8 (40%) | 2 (14%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| cryo | human-porcine | 0.5 deg | 8 (40%) | 4 (25%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) | 5 (29%) |
| cryo | human-porcine | 0.4 deg | 12 (50%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) | 8 (40%) |
| cryo | human-porcine | 0.3 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| cryo | human-porcine | 0.2 deg | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 | > 12 |
| rt | bovine-human | 60 deg | 5 (17%) | 2 (7%) | 5 (17%) | 3 (11%) | 2 (7%) | 2 (7%) | 2 (7%) |
| rt | bovine-human | 20 deg | 5 (17%) | 2 (7%) | 5 (17%) | 4 (14%) | 3 (11%) | 3 (11%) | 3 (11%) |
| rt | bovine-human | 10 deg | 5 (17%) | 4 (14%) | 5 (17%) | 4 (14%) | 5 (17%) | 5 (17%) | 4 (14%) |
| rt | bovine-human | 5 deg | 6 (19%) | 4 (14%) | 6 (19%) | 5 (17%) | 6 (19%) | 6 (19%) | 5 (17%) |
| rt | bovine-human | 2 deg | 12 (32%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) | 8 (24%) |
| rt | bovine-human | 1 deg | > 25 | 25 (50%) | > 25 | > 25 | > 25 | 25 (50%) | 25 (50%) |
| rt | bovine-porcine | 60 deg | 5 (17%) | 2 (7%) | 5 (17%) | 2 (7%) | 2 (7%) | 2 (7%) | 2 (7%) |
| rt | bovine-porcine | 20 deg | 5 (17%) | 2 (7%) | 5 (17%) | 3 (11%) | 3 (11%) | 3 (11%) | 3 (11%) |
| rt | bovine-porcine | 10 deg | 6 (19%) | 3 (11%) | 6 (19%) | 3 (11%) | 4 (14%) | 4 (14%) | 3 (11%) |
| rt | bovine-porcine | 5 deg | 6 (19%) | 3 (11%) | 6 (19%) | 6 (19%) | 5 (17%) | 5 (17%) | 5 (17%) |
| rt | bovine-porcine | 2 deg | 8 (24%) | 3 (11%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| rt | bovine-porcine | 1 deg | 8 (24%) | 4 (14%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| rt | human-porcine | 60 deg | 8 (24%) | 4 (14%) | 5 (17%) | 5 (17%) | 4 (14%) | 4 (14%) | 4 (14%) |
| rt | human-porcine | 20 deg | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) | 6 (19%) |
| rt | human-porcine | 10 deg | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) | 12 (32%) |
| rt | human-porcine | 5 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |
| rt | human-porcine | 2 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |
| rt | human-porcine | 1 deg | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 | > 25 |

Minorities found (k >= 2), all cells: cryo today 106/720 [0.12, 0.18], v6 416/720 [0.54, 0.61], v7 174/720 [0.21, 0.27], v8 241/720 [0.30, 0.37], v9 251/720 [0.31, 0.38], v10 252/720 [0.32, 0.39], v11 252/720 [0.32, 0.39]; room temperature today 416/1548 [0.25, 0.29], v6 971/1548 [0.60, 0.65], v7 490/1548 [0.29, 0.34], v8 718/1548 [0.44, 0.49], v9 796/1548 [0.49, 0.54], v10 797/1548 [0.49, 0.54], v11 811/1548 [0.50, 0.55]. Third collection (held out): today 26/110 [0.16, 0.33], v6 85/110 [0.68, 0.85], v7 26/110 [0.16, 0.33], v8 26/110 [0.16, 0.33], v9 58/110 [0.43, 0.62], v10 58/110 [0.43, 0.62], v11 58/110 [0.43, 0.62].

### False groups (reference: single-population trials)

| wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|
| 10 deg | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] |
| 1 deg | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] |
| 0.5 deg | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] |
| 0.4 deg | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] |
| 0.3 deg | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] |
| 0.2 deg | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] | 0/36 [0.00, 0.10] |
| all | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] | 0/189 [0.00, 0.02] |

| wedge | today | v6 | v7 | v8 | v9 | v10 | v11 |
|---|---|---|---|---|---|---|---|
| 60 deg | 0/9 [0.00, 0.34] | 2/9 [0.03, 0.60] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] | 0/9 [0.00, 0.34] |
| 20 deg | 2/27 [0.01, 0.24] | 7/27 [0.11, 0.46] | 2/27 [0.01, 0.24] | 2/27 [0.01, 0.24] | 2/27 [0.01, 0.24] | 2/27 [0.01, 0.24] | 2/27 [0.01, 0.24] |
| 10 deg | 2/36 [0.01, 0.19] | 7/36 [0.08, 0.36] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] |
| 5 deg | 2/36 [0.01, 0.19] | 4/36 [0.03, 0.26] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] | 2/36 [0.01, 0.19] |
| 2 deg | 4/36 [0.03, 0.26] | 5/36 [0.05, 0.29] | 4/36 [0.03, 0.26] | 4/36 [0.03, 0.26] | 4/36 [0.03, 0.26] | 4/36 [0.03, 0.26] | 4/36 [0.03, 0.26] |
| 1 deg | 1/36 [0.00, 0.15] | 8/36 [0.10, 0.39] | 1/36 [0.00, 0.15] | 1/36 [0.00, 0.15] | 1/36 [0.00, 0.15] | 1/36 [0.00, 0.15] | 1/36 [0.00, 0.15] |
| all | 11/180 [0.03, 0.11] | 33/180 [0.13, 0.25] | 11/180 [0.03, 0.11] | 11/180 [0.03, 0.11] | 11/180 [0.03, 0.11] | 11/180 [0.03, 0.11] | 11/180 [0.03, 0.11] |

Third collection (held out; paired against today's DIALS on the same single-population trials):

| procedure | minorities found (k >= 2) | single-population trials with a false group | of them where today's DIALS reports none | where only today's does |
|---|---|---|---|---|
| today | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | - | - |
| v6 | 85/110 [0.68, 0.85] | 0/30 [0.00, 0.12] | 0 | 0 |
| v7 | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | 0 | 0 |
| v8 | 26/110 [0.16, 0.33] | 0/30 [0.00, 0.12] | 0 | 0 |
| v9 | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 0 | 0 |
| v10 | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 0 | 0 |
| v11 | 58/110 [0.43, 0.62] | 0/30 [0.00, 0.12] | 0 | 0 |


### The frontier (detection, false groups, cost; every procedure run on every trial)

Cryo:

| procedure | kind | minorities found | false groups | median CPU s per trial | non-dominated |
|---|---|---|---|---|---|
| today | frozen | 106/720 | 0/189 | 2.31 | yes |
| v11 | frozen | 252/720 | 0/189 | 2.60 |  |
| v10 | frozen | 252/720 | 0/189 | 2.65 |  |
| v9 | frozen | 251/720 | 0/189 | 2.70 |  |
| v8 | frozen | 241/720 | 0/189 | 2.78 |  |
| v7 | frozen | 174/720 | 0/189 | 2.49 |  |
| v6 | frozen | 416/720 | 0/189 | 2.96 |  |
| v6_1 | frozen | 416/720 | 0/189 | 2.52 | yes |
| v5 | frozen | 407/720 | 0/189 | 2.36 | yes |
| v4 | frozen | 404/720 | 0/189 | 2.37 |  |
| v3 | frozen | 378/720 | 0/189 | 2.92 |  |
| v2 | frozen | 351/720 | 0/189 | 2.91 |  |
| v1 | frozen | 271/720 | 0/189 | 3.26 |  |
| v4_ratio_only | frozen | 264/720 | 0/189 | 2.36 | yes |
| v4_coherence_only | frozen | 403/720 | 0/189 | 2.37 |  |
| v4_bootstrap | frozen | 369/720 | 0/189 | 2.92 |  |
| v3_ratio_only | frozen | 270/720 | 0/189 | 2.92 |  |
| v3_coherence_only | frozen | 376/720 | 0/189 | 2.92 |  |
| v3_analytic | frozen | 431/720 | 5/189 | 2.45 | yes |
| v6_m3 | frozen | 315/720 | 0/189 | 2.96 |  |
| v6_m4 | frozen | 239/720 | 0/189 | 2.96 |  |
| v6_m6 | frozen | 146/720 | 0/189 | 2.96 |  |
| v6_m8 | frozen | 118/720 | 0/189 | 2.96 |  |
| today, floor2 (replay) | ablation | 261/716 | 104/184 | - |  |
| today, heuristic (replay) | ablation | 106/716 | 0/184 | - |  |
| today, unit (replay) | ablation | 110/716 | 0/184 | - |  |
| today, dims+1 (replay) | ablation | 95/716 | 0/184 | - |  |
| today, dims5 (replay) | ablation | 100/716 | 0/184 | - |  |
| today, floor2_dims+1 (replay) | ablation | 227/716 | 76/184 | - |  |

Room temperature:

| procedure | kind | minorities found | false groups | median CPU s per trial | non-dominated |
|---|---|---|---|---|---|
| today | frozen | 416/1548 | 11/180 | 13.58 | yes |
| v11 | frozen | 811/1548 | 11/180 | 15.73 | yes |
| v10 | frozen | 797/1548 | 11/180 | 16.04 |  |
| v9 | frozen | 796/1548 | 11/180 | 15.91 |  |
| v8 | frozen | 718/1548 | 11/180 | 16.44 |  |
| v7 | frozen | 490/1548 | 11/180 | 16.47 |  |
| v6 | frozen | 971/1548 | 33/180 | 46.76 |  |
| v6_1 | frozen | 971/1548 | 33/180 | 16.41 | yes |
| v5 | frozen | 960/1548 | 25/180 | 15.82 | yes |
| v4 | frozen | 959/1548 | 25/180 | 15.92 |  |
| v3 | frozen | 921/1548 | 25/180 | 45.52 |  |
| v2 | frozen | 794/1548 | 24/180 | 45.52 |  |
| v1 | frozen | 779/1548 | 22/180 | 51.48 |  |
| v4_ratio_only | frozen | 787/1548 | 19/180 | 15.86 |  |
| v4_coherence_only | frozen | 759/1548 | 20/180 | 15.50 | yes |
| v4_bootstrap | frozen | 918/1548 | 24/180 | 45.53 | yes |
| v3_ratio_only | frozen | 763/1548 | 23/180 | 45.57 |  |
| v3_coherence_only | frozen | 742/1548 | 17/180 | 45.53 |  |
| v3_analytic | frozen | 962/1548 | 28/180 | 15.79 | yes |
| v6_m3 | frozen | 792/1548 | 23/180 | 46.76 |  |
| v6_m4 | frozen | 633/1548 | 17/180 | 46.76 |  |
| v6_m6 | frozen | 446/1548 | 11/180 | 46.76 |  |
| v6_m8 | frozen | 427/1548 | 11/180 | 46.76 |  |

Development prototypes on cryo (detector alone, without OPTICS; cost not comparable):

| procedure | kind | minorities found | false groups | median CPU s per trial | non-dominated |
|---|---|---|---|---|---|
| prototype alt_amp (detector alone) | development | 437/720 | 7/189 | - |  |
| prototype alt_amp_c0.08 (detector alone) | development | 388/720 | 3/189 | - |  |
| prototype alt_amp_c0.1 (detector alone) | development | 343/720 | 1/189 | - |  |
| prototype alt_amp_ratio (detector alone) | development | 257/720 | 1/189 | - |  |
| prototype alt_ampw (detector alone) | development | 334/720 | 0/189 | - |  |
| prototype alt_nscore (detector alone) | development | 423/720 | 18/189 | - |  |
| prototype p5ucb_c4 (detector alone) | development | 281/720 | 4/189 | - |  |
| prototype p5ucb_c6 (detector alone) | development | 222/720 | 2/189 | - |  |
| prototype p6_r3 (detector alone) | development | 163/720 | 0/189 | - |  |
| prototype p6abs_0.03 (detector alone) | development | 356/720 | 0/189 | - |  |
| prototype p6b_r3 (detector alone) | development | 319/720 | 6/189 | - |  |
| prototype p6c_ana_r3.5 (detector alone) | development | 204/720 | 0/189 | - |  |
| prototype p6c_param_r3 (detector alone) | development | 218/720 | 2/189 | - |  |
| prototype p6c_param_r3.5 (detector alone) | development | 187/720 | 1/189 | - |  |
| prototype p6c_r3 (detector alone) | development | 234/720 | 1/189 | - |  |
| prototype p6c_r3.5 (detector alone) | development | 207/720 | 0/189 | - |  |
| prototype p6c_within_r3 (detector alone) | development | 205/720 | 1/189 | - |  |
| prototype p6d_gap_r3.5 (detector alone) | development | 217/720 | 3/189 | - |  |
| prototype p6e_ana_r3.5_a0.03 (detector alone) | development | 381/720 | 5/189 | - |  |
| prototype p6e_ana_r3.5_a0.03_q2.33 (detector alone) | development | 344/720 | 3/189 | - |  |
| prototype p6e_ana_r3.5_a0.03_t0.2 (detector alone) | development | 320/720 | 3/189 | - |  |
| prototype p6e_ana_r3.5_a0.03_t0.5 (detector alone) | development | 342/720 | 3/189 | - |  |
| prototype p6e_ana_r3.5_a0.03_test (detector alone) | development | 272/720 | 0/189 | - |  |
| prototype p6e_cal_r3.5_a0.03 (detector alone) | development | 417/720 | 4/189 | - |  |
| prototype p6e_r3.5_a0.03 (detector alone) | development | 342/720 | 0/189 | - |  |
| prototype p6e_r3.5_a0.035 (detector alone) | development | 337/720 | 0/189 | - |  |
| prototype p6e_r3.5_a0.03_fix (detector alone) | development | 342/720 | 0/189 | - |  |
| prototype p6e_r3.5_a0.05 (detector alone) | development | 301/720 | 0/189 | - |  |
| prototype p6e_r3.5_a0.05_invk (detector alone) | development | 376/720 | 0/189 | - |  |
| prototype p6e_r3.5_a0.06_invk (detector alone) | development | 373/720 | 0/189 | - |  |
| prototype p6e_rob0.5_r3.5_a0.03 (detector alone) | development | 286/720 | 1/189 | - |  |
| prototype p6e_rob_r3.5_a0.03 (detector alone) | development | 286/720 | 2/189 | - |  |
| prototype p6u_ana_r3.5_a0.03 (detector alone) | development | 375/720 | 2/189 | - |  |
| prototype p6u_ana_r3.5_a0.03_q3 (detector alone) | development | 316/720 | 1/189 | - |  |
| prototype p6u_ana_r3.5_a0.03_q3.5 (detector alone) | development | 292/720 | 1/189 | - |  |
| prototype p6u_ana_r3.5_a0.04_p0.001 (detector alone) | development | 314/720 | 0/189 | - |  |
| prototype p6u_ana_r3.5_a0.05 (detector alone) | development | 331/720 | 1/189 | - |  |
| prototype p6u_ana_r3.5_a0.05_p0.001 (detector alone) | development | 293/720 | 0/189 | - |  |
| prototype p6u_ana_r3.5_a0.05_q2 (detector alone) | development | 289/720 | 1/189 | - |  |
| prototype p6u_ana_r3.5_a0.05_q3 (detector alone) | development | 258/720 | 0/189 | - |  |
| prototype p6u_ana_v3rule (detector alone) | development | 415/720 | 1/189 | - |  |
| prototype p6u_ana_v3rule_p0.001 (detector alone) | development | 368/720 | 0/189 | - |  |
| prototype p6u_ana_v3rule_pair0.001 (detector alone) | development | 400/720 | 0/189 | - |  |
| prototype p6u_ana_v3rule_small3_0.001 (detector alone) | development | 388/720 | 0/189 | - |  |
| prototype p6u_cal_v3rule (detector alone) | development | 459/720 | 4/189 | - | yes |
| prototype p7_sparse (detector alone) | development | 416/720 | 10/189 | - |  |
| prototype p7_sparse_iter (detector alone) | development | 416/720 | 0/189 | - |  |
| prototype p7_sparse_local (detector alone) | development | 417/720 | 0/189 | - | yes |

### Limits (reference: the measurements of `cf/limits3.py`)

Misses of each frozen procedure by what limits them (minorities with k >= 2), across the collections
(amendment 18: one configuration for every sample, so a group of one species in either collection's
single-population data counts as structure):

Cryo:

| procedure | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | procedure |
|---|---|---|---|---|---|---|
| today | 106 | 264 | 161 | 8 | 12 | 169 |
| v1 | 271 | 257 | 103 | 8 | 8 | 73 |
| v2 | 351 | 255 | 60 | 8 | 8 | 38 |
| v3 | 378 | 254 | 51 | 8 | 8 | 21 |
| v4 | 404 | 254 | 36 | 8 | 7 | 11 |
| v6 | 416 | 252 | 31 | 8 | 7 | 6 |
| v7 | 174 | 264 | 161 | 8 | 12 | 101 |
| v8 | 241 | 264 | 161 | 8 | 12 | 34 |
| v9 | 251 | 264 | 161 | 8 | 12 | 24 |
| v10 | 252 | 264 | 161 | 8 | 12 | 23 |
| v11 | 252 | 264 | 161 | 8 | 12 | 23 |

Room temperature:

| procedure | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | procedure |
|---|---|---|---|---|---|---|
| today | 416 | 252 | 469 | 7 | 3 | 401 |
| v1 | 779 | 251 | 409 | 7 | 3 | 99 |
| v2 | 794 | 251 | 396 | 7 | 3 | 97 |
| v3 | 921 | 250 | 336 | 6 | 3 | 32 |
| v4 | 959 | 249 | 326 | 7 | 3 | 4 |
| v6 | 971 | 248 | 318 | 7 | 3 | 1 |
| v7 | 490 | 252 | 469 | 7 | 3 | 327 |
| v8 | 718 | 252 | 469 | 7 | 3 | 99 |
| v9 | 796 | 252 | 469 | 7 | 3 | 21 |
| v10 | 797 | 252 | 469 | 7 | 3 | 20 |
| v11 | 811 | 252 | 469 | 7 | 3 | 6 |

Per collection (a procedure told which kind of sample it has):

Cryo:

| procedure | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | procedure |
|---|---|---|---|---|---|---|
| today | 106 | 252 | 12 | 8 | 26 | 316 |
| v1 | 271 | 252 | 12 | 8 | 26 | 151 |
| v2 | 351 | 252 | 12 | 8 | 26 | 71 |
| v3 | 378 | 252 | 12 | 8 | 26 | 44 |
| v4 | 404 | 252 | 12 | 8 | 26 | 18 |
| v6 | 416 | 252 | 12 | 8 | 26 | 6 |
| v7 | 174 | 252 | 12 | 8 | 26 | 248 |
| v8 | 241 | 252 | 12 | 8 | 26 | 181 |
| v9 | 251 | 252 | 12 | 8 | 26 | 171 |
| v10 | 252 | 252 | 12 | 8 | 26 | 170 |
| v11 | 252 | 252 | 12 | 8 | 26 | 170 |

Room temperature:

| procedure | found | counting noise | structure inside one species | membership not resolved | noise known only from repeats | procedure |
|---|---|---|---|---|---|---|
| today | 416 | 252 | 450 | 7 | 3 | 420 |
| v1 | 779 | 251 | 392 | 7 | 3 | 116 |
| v2 | 794 | 251 | 379 | 7 | 3 | 114 |
| v3 | 921 | 250 | 319 | 6 | 3 | 49 |
| v4 | 959 | 249 | 313 | 7 | 3 | 17 |
| v6 | 971 | 248 | 306 | 7 | 3 | 13 |
| v7 | 490 | 252 | 450 | 7 | 3 | 346 |
| v8 | 718 | 252 | 450 | 7 | 3 | 118 |
| v9 | 796 | 252 | 450 | 7 | 3 | 40 |
| v10 | 797 | 252 | 450 | 7 | 3 | 39 |
| v11 | 811 | 252 | 450 | 7 | 3 | 25 |

The k90 the data allow against the procedures' (`FINDINGS.md` F6 for the full tables). Cells where v11's k90
is above what the data allow with the noise one experiment can know, across the collections: cryo
none; room temperature none;
per collection: cryo bovine-porcine at 1 deg (5 against 4); bovine-porcine at 0.5 deg (6 against 5); bovine-human at 1 deg (5 against 2); bovine-human at 0.5 deg (5 against 3); bovine-human at 0.4 deg (6 against 4); human-porcine at 1 deg (5 against 2); human-porcine at 0.5 deg (5 against 4); room temperature
human-porcine at 60 deg (4 against 3).

### Evidence

Every number above comes from commit `458cb6bcd63306d74e83e5fa45f492897c64f4a1`; the trial runners record the code commit of
every result line. Each configuration is judged on data it was not tuned on. v6: every procedure and
frontier point run on room temperature was frozen before any room-temperature result of any of them
existed (PREREGISTRATION.md amendments 3-13), and v6 failed that judgement on false groups.
v7, v8, v9, v10, v11 (amendments 17-26) were chosen after it, on both collections;
their judgement is the third collection (amendment 19: its labels verified against the labelled
room-temperature crystals before any trial; amendment 21: its pipeline stopped during processing so that
no trial ran before the last of them was frozen; amendment 25: its labels in the order the verification
settles; nothing tuned on it).

## What failed

- **v6 on room temperature, its held-out judgement**: it finds 971 of 1548 minorities (today's
  DIALS 416) but reports a group in 33 of 180 single-population trials (today's 11), 22 of
  them where today's DIALS reports none. Room-temperature crystals of one species hold sub-groups of
  two to four as distinct as a species minority there (`FINDINGS.md` F8); its thresholds, set on cryo where
  the species differ several times more, report them. Every other procedure frozen before room
  temperature raised it as well (the frontier above).
- **v1** (separation ratio only, amendment 3): missed minorities next to a majority with strong internal
  structure (porcine majorities), because relative to the majority's own structure they are no more
  distinct than the strongest group inside one species.
- **A flat coherence threshold** (v2, 0.03 then 0.05): 0.03 reported a porcine pair as a group in a
  single-population trial once the patch's own bootstrap seed was used (0.038, p = 0.004); 0.05 for
  every size missed groups of four to six next to structured majorities (coherence falls as 1/k).
- **The bootstrap over reflections** as counting noise (v1-v3): it overstates the noise of correlations
  over few shared reflections up to five-fold (F5), so small groups at small wedges never reached
  significance.
- **Calibrated noise without guards** (several prototypes, v3 with analytic noise): pairs of crystals of
  one species at small wedges stood out as strongly as species pairs and were reported (1-6 false groups
  on cryo).
- **Self-calibration of the noise from the residuals' spread** (scaled parametric replicates, robust
  scaling of the analytic model): wrong where real heterogeneity is present (full data) or the matrix is
  sparse (0.2 degrees).
- **Noise-scaled or look-elsewhere thresholds on coherence**: removed the false groups but gave back
  more detections than the stricter pair significance of v4.
- **A combined scan statistic** (gap or mean within) and a weaker guard on sparse data (development v6):
  no gain on cryo, one false group with the weaker guard; not adopted.
- **Amplitude correlations instead of DIALS's intensity correlations** (development after v6, cryo;
  amendment 14): the detector on Pearson's r of the square roots of the intensities, with v6's
  parameters, found 437 of 720 minorities but reported groups in 7 of
  189 single-population trials (Fisher's z understates the noise of amplitude
  correlations up to 1.7-fold, F5); with its coherence threshold doubled, 343 found and
  1 false group(s). Calibrated against the repeats, amplitudes separate almost exactly the
  minorities DIALS's matrix separates (F6), so the gain was the noise model's error, not information.
- **A detector for the sparse regime** (development after v6, cryo; `cf/proc7.py`): where v6's detector
  examines nothing (every trial at 0.2 degrees, nearly every one at 0.3), a scan of the inverse-variance
  statistic over the pairs that have an error, groups of three or more, look-elsewhere at 0.001, v6's
  rules: 416 minorities found (as v6) and groups in 10 of 189
  single-population trials; restricted to trials of ten or more and repeated up to three times,
  416 found and 0 false groups; with the search improved by single add, remove and swap moves (in
  the noise replicates alike), 417 found and 0 false groups. Searched blind,
  the most significant group there holds part of the minority with as many majority crystals, or is a
  strongly correlated pair elsewhere. One minority in 720 did not justify a new procedure untested on
  room temperature; not adopted.
- **A small-group rule with an absolute precision gate and coherence floor** (v8, amendment 20): the
  bovine-human and human-bovine pairs at 60 degrees on room temperature stand 60-100 standard errors out of
  their noise but have a coherence of 0.03-0.04, under its floor of 0.04 x 2/k; it leaves
  5 room-temperature cell(s) above the data's limit across the collections (bovine-porcine at 5 deg (6 against 5); bovine-human at 60 deg (3 against 2); bovine-human at 20 deg (4 against 3); bovine-human at 1 deg (> 25 against 25); human-porcine at 60 deg (5 against 4)). Replaced by a condition on how far the coherence stands
  out of its own noise (v9, amendment 22).
- **Merging the groups too small into the largest population before the take-out** (v7-v9): when the
  detector splits a trial into two populations of about half each, the merged one can come to hold more than
  half and stay as the background, losing the minority (`FINDINGS.md` F4); v10 decides the take-out first
  (amendment 23), as v11 does.
- **Comparators of the design's minority sizes only** (the limits, amendments 14-23): sub-groups of a minority
  of 7, 9 or 11 crystals had no single-population group to be compared with and counted as separable
  untested. Amendment 24 compares every size from 2 to n/2; on room temperature this moves
  7 minorities across the collections (7 per collection) from separable to a limit of the
  data, 7 of them misses of v11.
- **Defects found and fixed on the way**: the rank-1 fit stopped short of convergence (simultaneous
  update), a rejected candidate was merged back by residual (sending half of a varying population to a
  minority), correlations without an analytic error were kept with a near-zero spread, and a minority
  that was the detector's largest group was left inside the OPTICS cluster; and in the measurement of the limits, the test of a
  minority judged only the whole minority, while "found" accepts a cluster holding half of it
  (amendment 14; on cryo it moved individual verdicts, no cell's limit). In the harness that runs the
  frozen procedures from the saved arrays: a helper named like the take-out rule replaced it in the first
  v11 runs (they disagreed with v11's exact derivation from v6's groups, which is how it was found; rerun);
  an end-to-end run read the parameter file per trial while a later configuration's replaced it (rerun
  with the file pinned; for a while records of both runs sat in one directory, which showed two failures
  that the rerun does not have); one cause replay mixed the records of two runs (rerun clean).

## What remains open

- **Minorities attributed to the procedure** after all limits of the data, across the collections (one
  configuration for every sample): cryo 23 of 720 (8 at 1 deg, 8 at 0.5 deg, 1 at 0.4 deg, 6 at 0.2 deg); room temperature 6 of 1548 (1 at 5 deg, 3 at 2 deg, 2 at 1 deg);
  per collection: cryo 170, room temperature 25
  (`COMPLETE.md` section 7 lists them). Cryo: 17 of the 23 are
  groups of two to four that v6 reports, measured as v11 measures them standing 3-17 standard errors out
  of their noise with a standard error of 0.0029-0.0368; each fails both of v11's conditions
  (30 standard errors out and 1.5 times the rest's own structure, or a standard error of at most 0.002 and
  twice that structure), which the sub-groups of one species set (F8). 6 are where its detector
  examines nothing. Room temperature: 4 of the 6 are
  groups of two to four that v6 reports, measured as v11 measures them standing 23-64 standard errors out
  of their noise with a standard error of 0.0008-0.0024; each fails both of v11's conditions
  (30 standard errors out and 1.5 times the rest's own structure, or a standard error of at most 0.002 and
  twice that structure), which the sub-groups of one species set (F8). At 0.2 degrees per crystal (cryo), v11's
  detector examines nothing: no data set keeps three correlations with an error, the guard that keeps
  pairs of one species from being reported at small wedges. With their members known a group of them
  stands out of the noise; searched blind (two prototypes, What failed) the most significant group mixes
  them with majority crystals. Whether any blind procedure separates them is open.
- **How much the limits depend on the measurement** (`FINDINGS.md` F6): in a few cells the k90 the
  data allow moves by one to a few crystals under another statistic, without the structure test's
  tolerance, on amplitude correlations or on either matrix, always on minorities at the test's edge
  (p close to 0.001, two to four draws per size). A finer limit there needs more draws than the public
  collections hold, and a procedure using several matrices at once that keeps the false-group rate.
- **Noise known only from repeat measurements**: minorities that would be separable if the counting
  noise of correlations over 10-20 shared reflections were known as precisely as repeat measurements of
  the same crystals show it. Neither the analytic error nor a jackknife over reflections, which one
  experiment can compute, reproduces it (F5); whether a calibration learned on other crystals transfers
  is checked against room temperature's own repeats (F5).
- **Structure inside one species** is a property of the crystals, not of the procedure: groups of one
  species (cell size, single crystals with something in common at small wedges) are as coherent as a
  species minority; separating them would need information beyond the intensities' correlations.
- The size of the public collections (12 and 25 per species) bounds the largest minority fraction
  tested, and the single-population sample (189 cryo, 180 room-temperature trials, 30 of the third
  collection) how small a false-group rate can be shown: none in 189 is compatible with a true rate up to
  1.9% (95%), none in 30 with up to 12%. The third collection's single-population trials hold no sub-group
  that even v6 reports, so its judgement tests the false-group rate weakly (`FINDINGS.md` F10).
