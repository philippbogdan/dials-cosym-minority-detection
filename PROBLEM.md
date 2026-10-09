# The problem

Multi-crystal crystallography merges data from tens to hundreds of crystals, which are not always one population: a few may have bound a ligand, changed conformation or grown in another form, and merged blindly they are averaged away. DIALS, the open diffraction-processing package behind Diamond Light Source's automatic multi-crystal pipelines, separates populations by intensity: cosym places each dataset as a point fitted to the pairwise intensity correlations, and OPTICS density clustering labels clusters and outliers (`dials.correlation_matrix`, `xia2.multiplex`). It separates groups of comparable size. Nobody has measured how small a group it can find.

## To find

The smallest minority, as a number of crystals and as a fraction, that today's DIALS reports as a group of its own on real multi-crystal data, and how that depends on how different the minority is and how much data each crystal gives; which part of the procedure loses smaller groups; and a procedure that finds smaller groups without reporting groups in data holding one population.

## Given

- DIALS: https://github.com/dials/dials (BSD-3; v3.30.1, 17 September 2026).
- The method paper: Acta Cryst. D 81, 278-290 (2025), doi:10.1107/S2059798325004589, in data/. On 12 + 12 + 12 cryogenic crystals of bovine, porcine and human insulin it found three pure clusters with no outliers.
- Labelled raw data, CC BY 4.0, retrieved separately: Zenodo 13890874, those 36 crystals (6.6 GB); Zenodo 15062310, 15062327 and 15062343, 25 room-temperature crystals of each species, and 15062350, 19 more the paper left unclassified (about 60 GB together).

## Why it is not routine

The public collections cap how rare a real group can be, and subsets of the same crystals are not independent trials. The populations differ by one to three residues, so their signal sits close to the noise of the pairwise correlations. Dimension choice, outlier handling and minimum cluster size interact, and the defaults have changed since the paper. A group visible in the plots is not a group the procedure reports.
