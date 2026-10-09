# Pre-registration: definitions fixed before any detection was measured

Committed before any mixture or single-population trial was run with any procedure.
Later documents refer to these definitions; they are not changed after measuring. Anything
added later is marked as an amendment with its date and commit, and says whether any
result had been seen when it was written.

## Units

- A **crystal** is one rotation sweep from one Zenodo entry (`data/manifest.csv`).
- A **population** is one insulin species (bovine B, porcine P, human H) within one
  collection: `cryo` (Zenodo 13890874, 12 per species, 100 K, 10 degree wedges) or `rt`
  (Zenodo 15062310/27/43, 25 per species, 293 K, 60 degree wedges). The 19 crystals of
  15062350 ("unclassified") have no species label and are not used as labelled data.
- **Data per crystal** is the first `w` degrees of each sweep (images 1 to w/0.1), as in
  section 3.3.3 of the paper. Every crystal in a trial is cut to the same `w`.

## Trials

- A **mixture trial** is a set of crystals from one collection: all crystals of a majority
  species A (12 cryo, 25 rt) plus `k` crystals of a minority species B, all cut to `w`.
- A **single-population trial** is a set of crystals of one species of one collection, cut
  to `w`.
- Every procedure is run on exactly the same trials (paired design).
- The indexing ambiguity (I213 in an m-3m lattice) is resolved once per collection by
  `dials.cosym` on the full data of all crystals of that collection; trials select crystals
  and image ranges from that consistently indexed data. Each trial is then scaled on its
  own with `dials.scale` defaults, as `xia2.multiplex` does, and passed to the procedure.

## Draws

For an ordered pair (A, B) and size `k`, the crystals of B are shuffled with a fixed seed
and cut into floor(N_B / k) disjoint sets of `k`; each set is one draw. Minority sets never
share a crystal within a cell (collection, A, B, k, w). The same draws are used for every
`w` and every procedure. Both directions of each species pair are run; a cell's
"difference level" is the unordered pair. Intervals are exact (Clopper-Pearson, 95%) over
the draws of a cell. Draws in a cell share the majority, so intervals are conditional on
that majority; pooled numbers say so.

## Found (the one definition)

A minority M (|M| = k >= 2) is **found** in a mixture trial if the procedure reports a
cluster C (a reported group, not crystals labelled as outliers or noise) with

    |C ∩ M| >= max(2, k/2)   and   |C ∩ M| > |C| / 2.

That is: the reported group holds at least two and at least half of the minority, and
minority crystals are the majority of the group. Nothing else about the output counts
toward "found". A group visible in a plot but not reported is not found.

k = 1 cannot form a group; it is reported separately as **isolated** if the crystal is not
in any reported cluster that contains a majority crystal.

## False group (single-population data)

In a single-population trial the procedure **reports a false group** if it reports two or
more clusters. The **false-group rate** is the fraction of single-population trials with a
false group. Reports of zero clusters, and the number of crystals labelled as outliers, are
recorded separately.

## Detection limit

For a cell (collection, species pair, w) the detection probability at size k is the
fraction of draws found. The **detection limit** k90 is the smallest k >= 2 at which the
detection probability is at least 0.9 and stays at least 0.9 at every larger k tested. It
is also given as the fraction k90 / (m + k90), with m the majority size. k50 (same rule at
0.5) is a secondary summary. If no tested k reaches 0.9, the limit is "above the largest k
tested".

## Procedures and tuning

- **Today's DIALS** is the OPTICS "significant clusters" of `dials.correlation_matrix` from
  DIALS 3.30.1 with every default (the same code `xia2.multiplex` runs), on data scaled
  with `dials.scale` defaults.
- A new procedure is developed on the cryo collection only. Its configuration is frozen in
  a named commit before it is run on any rt trial; rt numbers for it are the held-out
  judgement. Cryo numbers for it are labelled in-sample.
- Cost is CPU seconds of the clustering stage per trial (scaled data in, labels out), on the
  same machine type, reported as the median over the trials of a cell.

## Limits

A failure to find a minority counts as a limit of the data only if an oracle that is told
which crystals form the minority cannot distinguish them from the same number of other
crystals in the same trial: a likelihood-ratio statistic on the merged intensities the
procedures use, calibrated by permuting crystal labels within the trial, does not reject at
the 5% level. Otherwise the failure is attributed to the procedure.

## Amendment 1 (before any trial; no procedure had been scored on any mixture or single-population trial)

Written after the data-quality pilot (`cf/pilot.py`, `results/pilot/`), which scaled all
labelled crystals of each collection at a range of wedge widths and summarised the
pairwise correlation coefficients; no clustering procedure was run in it. Two pipeline
checks had been run on three test trials (`results/design/test_cryo*.jsonl`) to debug the
runner; their outcomes are not used anywhere.

1. **Wedge grid.** cryo: w = 10 (full), 1, 0.5, 0.4, 0.3, 0.2 degrees (within-species CC
   holds at about 0.98 down to 0.5 and collapses between 0.3 and 0.2). rt: w = 60 (full), 20,
   10, 5, 2, 1 degrees (most crystals hold, a few have no diffraction in their first degrees).
   Minority sizes: cryo k = 1, 2, 3, 4, 5, 6, 8, 10, 12; rt k = 1, 2, 3, 4, 5, 6, 8, 10, 12, 25.
2. **More single-population trials.** Besides the first w degrees, single-population trials
   are also run on the later disjoint wedge segments 2, 3 and 4 (images 1+s*w/0.1 onwards,
   s = 1, 2, 3) where the sweep holds them. They add independent counting noise on the same
   crystals; false-group rates are reported for the first segment alone and pooled.
3. **Oracle statistic (replaces the chi-square of the original "Limits" section, whose
   self-terms let one crystal without diffraction inflate the null).** For crystals i, j of
   a trial, d_i and d_j are their merged intensities minus the inverse-variance weighted mean
   of all other crystals of the trial (excluding i and j), each divided by its standard
   deviation; c_ij is the correlation of d_i and d_j over the reflections both observe. For a
   group G, C(G) is the mean of c_ij over the pairs of G: how much the members deviate from the
   rest in the same way. The oracle is told M and compares C(M) with C(S) for the size-k
   subsets S of the majority (all of them, or 20 000 at random): p = (1 + #{C(S) >= C(M)}) /
   (1 + #S), using the same trial's data.
4. **Limits (replaces the original rule).** A failure to find M counts as a limit of the data
   if (a) the oracle does not reject at the 5% level, or (b) the majority of the same trial
   contains a decoy: a size-k subset S with C(S) >= C(M) (maximum over all subsets, or a
   lower bound found by local search when there are more than 2 million). (b) means the
   single population alone holds a group at least as coherent as the minority, so a procedure
   that reported M would report the decoy in single-population data. Otherwise the failure is
   attributed to the procedure.

## Amendment 2 (before any oracle, decoy or new-procedure result)

By this point today's DIALS had been scored on the cryo trials (its results do not depend on
this amendment) and the deviation-correlation matrix of one single-population cryo trial
(porcine, 10 degrees) had been inspected; no oracle p-value, decoy or new procedure had been
computed. That matrix showed the statistic of amendment 1 to be unfit for its purpose: being a
correlation of deviations it ignores their size, so at full data the small, continuous
crystal-to-crystal variation within one species gives pairwise values up to 0.88 and a
minority would be judged against decoys that are tiny in absolute terms; and two crystals
without diffraction share the same deviation (minus the true intensities), so they would form
a spurious decoy. Items 3 and 4 of amendment 1 therefore use, instead of c, the residual of
DIALS's own pairwise correlation matrix from a single population:

- rho is the weighted least-squares rank-1 fit to the trial's rij matrix off the diagonal,
  with DIALS's own weights wij (as cosym in one dimension); e_ij = r_ij - rho_i rho_j.
- S(G) is the mean of e_ij over the pairs of G (pairs DIALS could not calculate are left out):
  how much more the members of G agree with each other than one population predicts, in
  correlation units. A crystal without signal has rho near 0 and residuals near 0, so it
  neither inflates the null nor forms a decoy with another such crystal.
- The oracle p-value and the decoy are as in amendment 1, with S in place of C.

## Amendment 3 (procedure frozen; before any room-temperature result of it)

**Frozen procedure** (git tag `frozen-v1`): DIALS 3.30.1 with
`dials-patch/0001-correlation-find-small-groups.patch`, configured by
`dials-patch/frozen.phil`: today's OPTICS significant clusters, then the calibrated
small-group detector (bootstrap over reflections with 100 replicates, significance 0.01,
separation ratio 3.5, smallest judged population 5); groups other than the largest are taken
out of the OPTICS clusters. It was developed on the cryo collection only. No rt trial had
been run with it or any of its prototypes; rt had been looked at only through the
data-quality pilot and today's DIALS. Its rt numbers are the held-out judgement; its cryo
numbers are in-sample.

What was seen during development (cryo): today's DIALS scores, the oracle and decoy of
amendment 2, the stage ablations, and the scores of the prototypes P1-P7 (`results/dev/cryo`).
The configuration was chosen as the one with no false group in the 189 cryo single-population
trials and the most minorities found.

**Post-hoc refinement of "limits" (secondary; the pre-registered rules of amendments 1-2 stay
primary).** The oracle results showed that the supervised oracle "detects" a pair whose
coherence is 2.6 counting-noise standard deviations, only because it ranks first among the
66 majority pairs (p = 1/67 < 0.05); no unsupervised procedure with a low false-group rate can
report such a pair. A secondary limit is therefore reported: a minority lies within the
counting noise for any unsupervised procedure with false-group rate alpha = 0.01 if its
standardised block sum (residuals divided by bootstrap standard errors, summed over its pairs,
over sqrt(#pairs)) does not exceed the 99th percentile of the largest such sum over all
subsets of the same size in the trial's counting-noise replicates (`cf/limits2.py`). Results
under both rules are reported.

## Amendment 4 (procedure v2 frozen; refined limits; before any room-temperature result of a new procedure or of any limit analysis)

**What had been seen of rt.** The data-quality pilot (amendment 1; its within- and between-species
correlation summaries were tabulated again for this amendment: at 60 degrees the species differ by
about 0.03 (bovine-porcine), 0.02 (bovine-human) and 0.005 (human-porcine) in correlation, against
0.05-0.09 at cryo). Today's DIALS on rt was then partly done (about 2000 of 2628 trials); seen from it:
minorities of 2-4 never found, k90 of 5-6 at full data, and about 6 false groups in the first 112
single-population trials scored. No prototype, oracle, decoy, noise calibration or new procedure had
been run on any rt trial.

**What was learned on cryo since amendment 3** (all in `results/noisecal/cryo.json`,
`results/dev/cryo`, `results/limits/cryo_refined`):

1. Counting noise. The same crystals measured in disjoint wedge segments show that the bootstrap over
   reflections overstates the noise of the rij residuals where pairs share few reflections (median 9
   shared reflections per pair at 0.3 degrees; observed / predicted 0.2-0.4 at 0.2-0.3 degrees,
   0.5-0.9 at 0.4-0.5, about 1 at 1 degree). The large-sample error of Fisher's z with DIALS's own
   effective number of reflections of each pair (`cf.bootstrap.analytic_replicates`) agrees with the
   repeats within a factor 0.56-1.43 without any tuning; measurement errors alone (parametric) understate
   the noise up to sevenfold.
2. The separation ratio alone misses minorities next to a majority with strong internal structure
   (porcine majorities at 10 and 1 degrees): their coherence is 2-3 times that of the strongest subset
   of the majority, as large in relative terms as the strongest structure inside a single species. In
   absolute terms (correlation units) they differ: species minorities 0.03-0.10, structure inside one
   species at most 0.027 in the 189 single-population cryo trials.
3. Detector variants that used a calibrated noise model (analytic, or the bootstrap rescaled by the
   segment calibration) found more minorities at 0.3-0.5 degrees but reported 1-5 false groups (pairs
   of crystals of one species at 0.2-1 degree).
4. Two defects found with synthetic tests and fixed: the rank-1 fit (simultaneous update) stopped short
   of convergence when reliabilities differ widely (residuals off by up to 0.012 at 0.2 degrees); and a
   rejected candidate was merged back by mean residual, which can send half of a continuously varying
   population to a minority (now by correlation affinity).

**Frozen procedure v2** (git tag `frozen-v2`): DIALS 3.30.1 with `dials-patch/0001` + `0002`,
configured by `dials-patch/frozen.phil`: today's OPTICS significant clusters, then the small-group
detector of v1 (bootstrap over reflections, 100 replicates, significance 0.01, min population 5) in
which a significant candidate is a population of its own when its coherence exceeds 3.5 times the
strongest structure inside the rest **or is at least 0.03** (`min_coherence`), with the converged fit and
the affinity merge-back. Chosen on cryo as the configuration with no false group in the 189 cryo
single-population trials and the most minorities found (382 of 720 with today's OPTICS, against 268
for v1 and 106 for today's DIALS) among those tried (`results/dev/configurations.txt`). Its rt numbers
are the held-out judgement.

**Also run on rt, frozen now** (`cf/newproc.py` CONFIGS), as frontier points: v1 exactly as frozen in
amendment 3 (`cf/frozen/groups_v1.py`); v2 without `min_coherence`; v2's coherence floor alone (the
absolute criterion; expected not to carry over, since rt species differences are up to ten times smaller
in correlation); v2's rule with the analytic noise model (cheaper).

**Refined limits** (secondary; the rules of amendments 1-2 stay primary; this replaces the post-hoc rule
of amendment 3, whose noise came from the bootstrap shown above to overstate it), `cf/limits3.py`:
counting noise from the analytic model multiplied by the factor measured against repeat segments for the
collection and wedge (1 where no second segment exists). A minority not found is a limit of the data if
(a) its standardised block sum does not exceed counting noise at 0.01 with look-elsewhere over all
subsets of its size in the trial, or (b) some single-population trial of the same collection and wedge
holds a group of the same size, standing out of counting noise at 0.01 and leaving a rest of at least 5
crystals, that is at least as coherent and at least as separated (coherence over the strongest structure
in its rest) as the minority is in its trial. Otherwise the miss is attributed to the procedure. Per cell
(collection, species pair, wedge) the data's k90 is computed from the minorities that are not data
limits, and compared with each procedure's k90.

## Amendment 5 (v2 re-frozen with min_coherence 0.05; still before any room-temperature result of a new procedure or of any limit analysis)

Run with the DIALS module itself and the patch's own bootstrap seed (`cf/newproc.py`, all 1341
cryo trials), v2 as frozen in amendment 4 reported one false group in the 189 cryo
single-population trials (the prototype, on differently seeded replicates, had none): two porcine
crystals of the small-cell sub-population at 1 degree, third segment (PIX15, PIX16), coherence
0.038, p = 0.004 against counting noise. That is the strongest significant group of one species
seen on cryo (the next are 0.027 and 0.025). With the module and its seed, min_coherence 0.035
still gave that false group, 0.04, 0.045 and 0.05 gave none and found 366, 359 and 351 of the 720
cryo minorities. **min_coherence is set to 0.05** (a margin above 0.038; the false-group rate is
the criterion that must not rise), in `dials-patch/0003`, `dials-patch/frozen.phil` and
`cf/newproc.py`; everything else of amendment 4 stands (git tag `frozen-v2.1`). The frontier
points of amendment 4 that use v2's floor use 0.05. No rt trial had been run with any new
procedure; nothing more of rt had been seen than amendment 4 lists.

## Amendment 6 (procedure v3 frozen; still before any room-temperature result of a new procedure or of any limit analysis)

The refined limits of amendment 4 (`cf/limits3.py`) on cryo, with v2 (amendment 5) and every
variant tried since, attributed most remaining misses to the procedure at 1 degree and below. Two
causes were measured on cryo only:

1. **Size.** A group's coherence falls about as 1/k with its size k (the rank-1 fit absorbs more of a
   larger group), and so does the strongest structure inside one species: 0.038 for two crystals,
   0.025 for three, 0.012 for five (cryo single-population trials, significant against counting
   noise). A threshold of 0.05 for every size therefore missed groups of four to six crystals next to a
   structured majority (bovine-porcine at 1 degree: coherence 0.027 for five crystals).
2. **Counting noise at small wedges.** Where pairs share few reflections (median effective number 7-23
   at 0.3-0.5 degrees) both noise models available to a procedure without repeat measurements
   overstate the noise (bootstrap 1.5-3 times, analytic 1.2-1.5 times, `results/noisecal/cryo.json`);
   rules that used the analytic model found more minorities there but reported crystals of one species
   that stand out of the noise as strongly as species minorities (pairs involving porcine PIX5 at 0.5
   degrees, coherence 0.07-0.095) as groups, and every guard against that (noise-scaled floors,
   stricter significance for the floor, look-elsewhere floor tests) gave the gain back.

**Frozen procedure v3** (git tag `frozen-v3`): v2 with the threshold min_coherence * 2/k for a group of
k (0.05 for a pair, `dials-patch/0004`), everything else unchanged. On cryo it found 383 of 720 minorities with no false group on the development replicates (355 for the
same rule with a flat threshold), and 378 of 720 with no false group run as the DIALS module with the
patch's own bootstrap seed (`cf/newproc.py`; v2: 351). Its rt numbers are the held-out judgement; v2 (amendment 5) and v1
(amendment 3) and the frontier points of amendment 4 (now on v3's rule: separation alone, threshold
alone, analytic noise) are also run on rt, all frozen now (`cf/newproc.py` CONFIGS).

**Refined limits, two corrections (before any rt limit result).** (a) Coherence and separation are
measured with counting noise; a minority counts as within the structure of one population when some
qualifying single-population group (amendment 4) is at least as coherent and as separated **within two
standard errors** of the difference (counting noise of both, from the calibrated replicates); at full
data this is the strict comparison of amendment 4, at small wedges it stops treating differences
smaller than the noise as real. (b) The counting-noise test scores only the pairs inside the minority
and so is weak for large minorities (it called minorities of 10-12 at 0.3 degrees noise that today's
DIALS finds); a minority that one of the frozen procedures with no more false groups than today's DIALS
finds in that trial is not counted as a counting-noise limit. Both rules are reported with the strict
version of amendment 4 alongside.

## Amendment 7 (refined limits only; before any room-temperature limit result)

Two completions of the refined limits of amendments 4 and 6, found while checking them on cryo:
(a) the single-population comparators are the 20 most coherent groups of every size in every
single-population trial (amendment 4 took only the most coherent one; a procedure is just as liable
to report the others); (b) the counting-noise test of a minority is made against the noise
replicates' maxima over all sizes as well as all subsets (standardised per size, Gumbel tail over
sizes, as the procedures' own scan), because an unsupervised procedure knows neither the members nor
the size of a group. Both are reported; the known-size test of amendment 4 is kept in the records
(`p_noise`) beside the new one (`p_noise_any_size`).

## Amendment 8 (procedure v4 frozen; before any room-temperature result of a new procedure or of any limit analysis)

With v3 (amendment 6) the refined limits of amendments 4, 6 and 7 still attributed misses on cryo to
the procedure at 0.3-0.5 degrees, mostly groups of three or more crystals that stand clear of the
counting noise (measured with repeat segments) while no group of one species of that size does at
that or any larger wedge. The cause is v3's noise model: the bootstrap over reflections overstates
the noise there 1.5-3 fold. A defect was found and fixed on the way: pairs whose correlation has too
few reflections for an analytic error got the small spread of the refit instead of being left out
(prototype scan and limit analysis; the bootstrap-based results are unaffected).

**Frozen procedure v4** (git tag `frozen-v4`; `dials-patch/0005`, `dials-patch/frozen.phil`): v3 with
the analytic counting noise (Fisher's z with DIALS's effective number of reflections per pair, no
recomputation of correlations), correlations without an error left out and data sets with fewer than
3 usable correlations not judged, and p < 0.001 for groups of two (the only false group of the
analytic variants on cryo after the other changes was a pair of porcine crystals at 0.5 degrees, p =
0.0011, coherence 0.095). On cryo, run as the DIALS module with its own seed: 404 of 720 minorities
found, no false group in the 189 single-population trials (v3: 378, today's DIALS: 106). Its rt
numbers are the held-out judgement. Also run on rt, frozen now (`cf/newproc.py` CONFIGS): v3, v2, v1,
and frontier points v4 by separation alone, v4 by coherence threshold alone, v4's rule with the
bootstrap noise, v3's frontier points. Nothing of rt beyond what amendment 4 lists had been seen.

## Amendment 9 (procedure v5 frozen; before any room-temperature result of a new procedure or of any limit analysis)

v4 on cryo found three minorities that it then failed to report: the detector had split the majority
into several groups, so the minority was its largest group (or tied for largest) and kept the OPTICS
labels. **Frozen procedure v5** (git tag `frozen-v5`, `dials-patch/0006`): v4 with
`take_out_groups`: groups other than the largest are taken out of the OPTICS clusters, and when no
group holds more than half of the data sets every group is. On cryo (DIALS module, its own seed): 407
of 720 minorities found, no false group in the 189 single-population trials. Its rt numbers are the
held-out judgement; v4 and the earlier procedures and frontier points of amendment 8 are run on rt
too. Nothing of rt beyond amendment 4's list had been seen.

## Amendment 10 (refined limits only; before any room-temperature limit result)

(a) Structure inside one population seen with at least as much data per crystal counts for smaller
wedges too (a difference between crystals does not vanish with less data; it only gets noisier): the
comparators of a minority at wedge w are the qualifying single-population groups at w and at every
larger wedge of the same collection. (b) The counting noise used for the limit is measured with
repeat measurements (amendment 4's calibration factor). A single experiment has no repeats; the noise
it can establish for itself is the analytic model at face value (calibration measured to be 0.56-1.43
of the truth at 0.2-1 degree on cryo, and no estimator from one experiment's own residuals did
better at 0.4-0.5 degrees). A minority that stands out of the repeat-calibrated noise but not of the
model noise (same test, `p_noise_model`) is reported in a category of its own, "beyond the noise only
with repeat measurements", distinct from both the data limits and the procedure's misses.

## Amendment 11 (refined limits only; before any room-temperature limit result)

The pre-registered rule (amendments 1-2) counts a minority as a limit of the data when the majority of
its own trial holds a group of the same size at least as coherent; the refined rule had dropped that
part. It is restored in the refined form: the most coherent group of the minority's size inside the
trial's majority, standing out of the (repeat-calibrated) counting noise at 0.01, at least as coherent
as the minority within two standard errors, makes the minority a limit of the data ("structure inside
one species"), unless a frozen procedure with no more false groups than today's found it.

## Amendment 12 (procedure v6 frozen; before any room-temperature result of a new procedure or of any limit analysis)

With the refined limits, the minorities v5 missed on cryo that no measurement puts in the data were,
apart from four that no procedure found, found by the bootstrap-based procedures (v1-v3): mostly
pairs at 0.4-0.5 degrees, which v5's pair threshold (p < 0.001 with the analytic noise) rejects. The
two noise models rank pairs differently: the bootstrap discounts correlations resting on a few strong
reflections (the porcine pairs that forced the pair threshold are not significant under it), the
analytic error is better calibrated where pairs share few reflections. **Frozen procedure v6** (git tag
`frozen-v6`, `dials-patch/0007`, `significant_clusters.groups.noise_model=both`): the detector of v5
run with the analytic noise (as v5) and with the bootstrap (no pair threshold, every data set judged);
the groups of the first, plus the groups of the second that share no data set with them, are taken out
of the OPTICS clusters. On cryo (DIALS module, its seeds): 416 of 720 minorities found, no false group
in the 189 single-population trials (v5: 407; the union of every frozen procedure without more false
groups than today's: 416). Its rt numbers are the held-out judgement; v5 and the procedures and frontier
points of amendments 8-9 are run on rt as well. Nothing of rt beyond amendment 4's list had been seen.

## Amendment 13 (refined limits only; before any room-temperature limit result)

The level of the refined limit's counting-noise tests is set to 0.001 (it was 0.01): a procedure that
reported a group whenever it stood out of the noise at 0.01 would be expected to report about two
groups from noise alone in the 189 cryo single-population trials, raising the false-group rate above
today's (0 of 189), which the problem forbids; at 0.001 the expectation is below 0.2. The comparators
of "structure inside one species" are held to the same level. The results with the tests at 0.01
and the strict comparison of amendment 4 are reported alongside ("strict").

## Amendment 14 (refined limits only; before any room-temperature refined-limit result)

Seen of room temperature before this amendment: today's DIALS on every trial, the pre-registered
oracle and decoy limits (amendment 2's rule), and the frozen procedures v1-v6 on 287-290 of the
2628 trials (mixtures only; how many minorities each found per wedge). No room-temperature
single-population result of a new procedure and no refined-limit result.

A development detector on amplitude correlations (Pearson's r of sqrt(I) over the shared reflections,
unweighted; `results/dev/cryo/alt_*`, not adopted: it reported groups in 1-18 of the 189 cryo
single-population trials depending on its thresholds) found human-porcine minorities at 0.3 degrees
that the refined limits put within the counting noise. Two defects of the limit test showed:

1. It judged only the whole minority, while "found" accepts a cluster holding at least max(2, k/2) of
   the minority's members as its majority (that detector had isolated three of four). The test now also
   judges, for every size from max(2, ceil(k/2)) to k - 1, the sub-group of the minority with the
   largest standardised block sum and the one with the largest inverse-variance statistic (all
   subsets when there are at most 5000, otherwise backward elimination). A minority is a data limit
   only if it and every such sub-group are, each by the same tests as before sized to the sub-group:
   counting noise if all are within the noise, structure inside one species otherwise.
2. Its statistic, the standardised block sum, weights every pair's z equally, which is optimal when
   the excess correlation scales with the pair's noise; for an excess common to the pairs the
   inverse-variance weighted sum, sum(z_ij / se_ij) / sqrt(sum(1 / se_ij^2)), is. Both are computed,
   each against its own look-elsewhere null over subsets and sizes; a group is within the noise only
   if neither stands out at 0.001.

On cryo neither change moves the k90 the data allow in any cell; the minorities v6 misses that the
data allow go from 1 to 4. The detector's finds at 0.3 degrees that remain within the noise stand out
of the repeat-calibrated noise at p = 0.01-0.1 (look-elsewhere): a procedure reporting groups that
strong reports groups from noise alone in several percent of single-population trials (about 1.4 of
the 36 at 0.3 degrees expected); it reported one in the 189.

Checked alongside, not part of the rule: the same analysis on the amplitude correlation matrix
(design suffix `_amp`; Fisher's z with the number of shared reflections, calibrated against repeat
measurements: observed / model 1.0-1.7). On the 600 cryo minorities with arrays here (10 degrees
excepted) the verdicts agree for 593; amplitudes alone separate five (four pairs at 0.5 degrees,
p 4e-5 to 0.001, which on DIALS's matrix are within one species' structure) and intensities alone two.
This moves one cell, bovine-human at 0.5 degrees (k90 2 with amplitudes, 3 on DIALS's matrix), on one
minority at p = 0.00099. Room-temperature limits are computed on both matrices.

Addendum to amendment 14 (same day, still before any room-temperature refined-limit result). An
independent review of the limit analysis (development notes, not part of the results) pointed out that
"found" also accepts a cluster holding majority crystals, fewer than its minority members. The test now
also judges such groups: from the minority and each of its qualifying sub-groups, majority crystals are
added one at a time (each time the one that raises the group's standardised statistic most, for each
statistic) while they stay fewer than the minority members; the extension that stands out most is
judged, by the same tests, when it stands out more than its base. Selection check (`cf/devselect.py`):
random pseudo-minorities of 2-6 crystals in the cryo single-population trials, searched the same way,
stand out at 0.001 in none of 1680 draws at 0.2-0.3 degrees (where those trials hold little structure
of their own); at 0.5-10 degrees they stand out through the species' own structure, which the
structure test judges. On cryo the change moves no cell's k90 the data allow; the minorities v6
misses that the data allow go from 4 to 16 (11 of them at 0.2 degrees). The review's other points
are reported with the results (the two-standard-error tolerance of the structure test: k90 under
tolerance 0 and 1 alongside; the amplitude matrix: alongside; a noise calibration learned on other
crystals: tested on room temperature with cryo's factors).

## Amendment 15 (refined limits only; before any room-temperature refined-limit result)

At 0.2 degrees v6's detector examines no candidate in any trial (no data set keeps three correlations
with an error), while the limit test finds, with the members known, qualifying groups of 16 cryo
minorities that stand out of the noise. Two prototypes for this regime (`cf/proc7.py`: the
inverse-variance scan over the pairs that have an error, groups of three or more, look-elsewhere at
0.001, v6's rules; once per trial, or up to three times; not adopted) found none of them: searched
blind, the most significant group holds three or four of the minority with as many or more majority
crystals, which "found" does not accept. The limit test now also records the group a search that does
not know the members finds (`blind_best`: for each statistic, greedy growth from every pair over all
sizes, standardised per size by the noise maxima). A minority that would otherwise be separable is
**membership not resolved** (a data limit) when, for both statistics, that group overlaps the minority,
would not count as finding it, and stands out more than every group that would: the data then favour a
grouping that mixes the minority with majority crystals. On cryo: 7 of the 16 (all at 0.2 degrees),
and one minority at 0.4 degrees that was separable only with the noise known from repeats; the other 9
stay attributed to the procedure; no cell's k90 moves.

## Amendment 16 (refined limits only; before any room-temperature refined-limit result)

The repeat measurements calibrate the counting noise per wedge and species (`cf/noisecal.py`), and the
species differ: at 0.2 degrees observed / analytic is 1.43 for porcine pairs and 0.76-0.79 for bovine
and human, at 0.4-0.5 degrees 0.56-0.62 for bovine and human and 0.87-0.92 for porcine. The limit test
used the mean over species, understating the noise of porcine pairs and overstating the others'. Each
pair's noise is now scaled by the geometric mean of its two crystals' species factors at that wedge (the
wedge's mean for a species without its own; 1 where the wedge has no repeat segments); the
inverse-variance weights use the calibrated error. The test knows the crystals' species, as a
measurement of the data may; no procedure is given them. On cryo: v6's misses that are counting noise
252 (unchanged), structure inside one species 12 (was 18), noise known only from repeats 26 (was 17),
membership not resolved 8, procedure 6 (was 9, all at 0.2 degrees); no cell's k90 with the noise one
experiment can know moves (with the noise from repeats, bovine-human 0.4 and human-porcine 0.5 degrees
would allow 3 instead of 4). `CF_CAL_BY_SPECIES=0` gives amendment 15's test.
The verdict with the noise one experiment can know applies every check with that noise, the blind
search of amendment 15 included (recorded as `blind_model`). The six minorities left to the procedure
are all at 0.2 degrees, where v6's detector examines nothing; in each, a group a blind search finds is
about as significant as the minority's best qualifying group (within 0.1 in the standardised statistic)
or a strongly correlated pair of majority crystals stands out first. A third prototype for that regime
(the blind search improved by single add, remove and swap moves, the same search in the noise replicates;
`cf/proc7.py` `sparse_groups_local`) found one minority more than v6 (417 of 720), none of these six, and
no false group.

## Amendment 17 (after the room-temperature judgement of v6; nothing after it is held out)

**The held-out judgement.** On the 2628 room-temperature trials the frozen v6 finds 971 of the 1548
minorities (today's DIALS 416) but reports a group in 33 of the 180 single-population trials (today's
11): in 22 trials where today's DIALS reports none, and in none the other way round. v6 raises the
false-group rate on the collection it was not tuned on, which the problem forbids. Every frozen
procedure and frontier point raised it too (6 to 22 trials more than today's).

**What was looked at afterwards (room temperature is no longer held out for anything below).** The
groups v6 adds in room-temperature single-population trials hold 2 (12 groups), 3 (7) or 4 (6)
crystals; they are sub-groups of one species as coherent as a species minority at room temperature
(porcine image_170462, image_170464, image_170482 and image_170540 together at p = 3e-11, two of which
join the human cluster in today's DIALS at full data; human groups of two to four at p down to 1e-16);
the paper notes "clear substructures within the porcine and human insulin groups" (section 3.2.2).
A coherence floor on the separation-ratio rule, a stricter significance for every size, requiring both
rules, and a minimum group size were replayed on v6's candidates (development notes, in-sample on both
collections). Only a minimum group size removed every added false group: with groups of fewer than 5
crystals returned to the background population (derived exactly from v6's detector groups,
`cf/devfloor.py`, results `v6_m3`-`v6_m8`), the false-group rate is today's on both collections (0 of
189, 11 of 180, the same trials) and 174 of 720 cryo and 490 of 1548 room-temperature minorities are
found (today's 106 and 416).

**v7** (`dials-patch/0010`, `groups.min_group_size = 5`, the default; `cf/newproc.py` config `v7`):
v6 with groups smaller than five returned to the background population. Five is OPTICS's own
smallest cluster here (min_samples >= 5): the change then reports no group smaller than DIALS reports
without it. The value was chosen after seeing room temperature; its room-temperature numbers are
in-sample. A check of the rule "the smallest size that adds no false group on the data seen":
chosen on cryo and the room-temperature trials of half the crystals (half0), it gives 4, which adds 2
false groups on the other half's (half1) single-population trials: a size learned on some crystals
does not protect against the sub-groups of others. v6 remains the pre-registered procedure and its
held-out judgement stands; with `min_group_size = 2` the change is v6 (`v6_1`, 0001-0009, the same
decisions).

## Amendment 18 (limits; after v6's room-temperature procedure results, before any room-temperature limit result)

A clustering procedure run automatically has one configuration for every sample and is not told what
kind of sample it has. For that question the structure test is also run **across the collections**:
a minority is within the structure inside one species when a qualifying single-population group of
either collection, of its size (any wedge of the other collection, as much or more data per crystal in
its own), dominates it as in amendments 6-16 (coherence and separation within two standard errors,
standing out of the noise at 0.001, a rest of at least five): a procedure deciding by coherence and
separation that reported the minority would report that group in that collection's single-population
data. Both classifications are reported: per collection (what a procedure tuned to that kind of sample
can reach) and across the collections (what one configuration for both can reach); COMPLETE.md states
the verdict for each.

## Amendment 19 (a third collection, untouched: the held-out judgement of v7 and a second one of v6)

Before any of its data were downloaded or processed. The paper's section 3.1 used a third public,
labelled collection that PROBLEM.md does not list: 16 crystals measured at room temperature on VMXi
(20 degree sweeps, 293 K), Zenodo 15077303 ("bovine and human", CC BY 4.0, 4.2 GB). Its supporting
information (Figure S1) states that data sets 0-7 are bovine and 8-15 human; the deposit names its
crystals by image number, in two batches of eight (image_140500 ... image_140818 and image_143421 ...
image_143867), taken in that order as bovine and human. **Verification of the labels before any trial**:
each crystal's mean correlation (DIALS's pairwise matrix of a joint run with the labelled room-temperature
crystals of 15062310 and 15062343) must be higher with its labelled species than with the other; a crystal
that fails is left out, and this is reported.

**Design** (collection `rtbh`, as `rt` with eight crystals per species): mixtures of all eight crystals of
one species and k = 1-8 of the other, draws of the minority disjoint within a cell, both directions;
single-population trials of all eight crystals of one species, on the first w degrees and on the later
disjoint segments (up to three more); w = 20, 10, 5, 2 and 1 degrees.

**Procedures judged**, frozen before any result on it: today's DIALS; **v7** (tag `frozen-v7`,
`dials-patch/0001-0010`, `cf/frozen/groups_v7.py`); **v6** (tag `frozen-v6`). No other procedure is run on
it before the judgement, and nothing is tuned on it. **Criteria**: the false-group rate in its
single-population trials must not be above today's (paired: no trial where the procedure reports a group
and today's DIALS does not); the minorities found, against today's, by cell (k90) and in total, with
exact intervals. With eight crystals per species the smallest group v7 can report is five, and the rest a
group is split from must hold five: the collection tests v7's false-group rate and its gain for groups of
five to eight.

## Amendment 20 (v8 frozen; before any result on the third collection)

The minimum group size of v7 gives up small minorities that the data separate from every
single-population group of either collection. Mapped on the limit records (oracle level, both
collections): no rule on a group's coherence and separation takes any group of 2-4 crystals without also
taking a single-population group, unless the coherence is precisely measured; with a precision
condition such rules exist. Measured as the procedure itself can (`cf/devsmall.py`: DIALS's analytic noise at
face value, the detector's coherence and decoy), on v6's groups of two to four in every cryo and
room-temperature trial, a search over 256 rules (standard error of the coherence at most 0.001-0.005;
coherence minus z = 0-3 standard errors at least 0.03-0.06 x 2/k; separation ratio at least 1.5-5) for
the most minorities found with no single-population trial reporting a group where today's DIALS does not:
the best takes 249 cryo and 741 room-temperature minorities. **Frozen procedure v8** (tag `frozen-v8`,
`dials-patch/0001-0011`, `groups.small_groups`: max_se 0.002, z 2, min_coherence 0.04, min_ratio 2, a
neighbour of the best with a margin; `cf/frozen/groups_v8.py`; `dials-patch/frozen.phil`): v7, and a group
smaller than five when its coherence is measured that precisely and stands that clearly above the floor
and the rest's structure. In-sample on both collections: 241 of 720 cryo and 718 of 1548
room-temperature minorities, no single-population trial with a group where today's DIALS has none (v7:
174 and 490). It is judged with v7 and v6 on the third collection (amendment 19), of which no result
exists yet: it is the held-out judgement of v8 and v7.

## Amendment 21 (the third collection's trials held; before any trial on it)

v8 measured on both collections: on cryo its run (241 of 720 minorities, 0 of 189 single-population
trials with a group); on room temperature derived exactly from v6's groups (`cf.newproc` run pending),
718 of 1548 and 11 of 178 trials, today's rate. Across the collections (amendment 18) its k90 is the
data's in every cryo cell but above it in four room-temperature cells (bovine-human at 60 and 1
degrees, bovine-porcine at 5, human-porcine at 10), with 37 room-temperature and 23 cryo minorities
separable but missed. At 60 degrees the pairs it misses fail only its coherence floor (0.04 x 2/k, the
margin taken in amendment 20 from the best rule of that search). A further configuration (v9) is
developed on cryo and room temperature and judged with v6, v7 and v8 on the third collection. So
that it is held out for v9 too, the third collection's pipeline (`pipeline/rtbh_all.sh`) was stopped
during processing (16:33, 8 of its 16 crystals processed, no trial run, no procedure run on it);
`pipeline/rtbh_prepare.sh` continues processing, joint cosym and the label check, and
`pipeline/rtbh_judge.sh` runs its trials and every judged procedure once v9 is frozen. Nothing of the
third collection is used to develop v9.

## Amendment 22 (v9 frozen; before any trial on the third collection)

Developed on cryo and room temperature only (amendment 21), on v6's groups of two to four in every
trial measured as the procedure measures them (`cf/devsmall.py`: coherence S, its standard error se
over DIALS's analytic counting-noise replicates, the rest's strongest structure of the size). Every
group of two to four v6 reports in a room-temperature single-population trial has S at most 24.9
standard errors (31 groups); the room-temperature minorities v8 misses at 60 degrees stand 60-100
standard errors out, and so do the cryo minorities v8 finds at 10 degrees. Two rule families were
searched, in-sample, for the most minorities found with no single-population trial reporting a group
where today's DIALS does not (development notes): S/se >= z0 with S - 2 se >= c0 x 2/k and
S/decoy >= r0 (140 rules; with no added group from z0 = 25, the best 251 cryo and 825 room-temperature
minorities), and S - z se >= c0 (x 2/k) with S/decoy >= r0 and no precision condition (210 rules; the
best, 262 and 806, is within 5% of a human sub-group at 20 degrees). Chosen with a margin: the first
form at z0 = 30, c0 = 0.01, r0 = 1.5 (every neighbour in the grid adds no group; groups appear at
z0 = 20).

**Frozen procedure v9** (tag `frozen-v9`, `dials-patch/0001-0013`, `cf/frozen/groups_v9.py`,
`dials-patch/frozen.phil`; v8's parameters in `frozen-v8.phil`): v7, and a group smaller than five
when its coherence stands 30 standard errors out of the counting noise, 2 standard errors above
0.01 x 2/k and 1.5 times above the strongest structure of its size inside the rest
(`groups.small_groups`, `small_group_admitted`; 0013 formats the change with DIALS's ruff
configuration, no change in behaviour). In-sample, derived exactly from v6's groups: 251 of 720 cryo
and 796 of 1548 room-temperature minorities (v8: 241 and 718), no single-population trial with a group
where today's DIALS has none (0 of 189; 11 of 178, today's). Across the collections its k90 is the
data's in every cryo cell; on room temperature v8's cells bovine-human at 60 and bovine-porcine at 5
degrees close, and two remain above: human-porcine at 10 degrees (12 against 8: one draw each at
k = 8 and 10 that v6's detector does not report; under amendment 13's rule the data's k90 there is 12)
and bovine-human at 1 degree (one of the two draws at k = 25).

**Judgement.** v9 is the change as proposed. It is judged with v8, v7 and v6 on the third
collection, by the criteria of amendment 19, on trials none of which had run when it was frozen. If
v9 raises the false-group rate there (a single-population trial with a group where today's DIALS has
none), it fails its judgement, which is reported, and the change is proposed with the latest of v8
and v7 that does not, stated as chosen after the judgement.

## Amendment 23 (v10 frozen; before any trial on the third collection)

Checked on cryo and room temperature only: of v9's two room-temperature cells above the data's limit
across the collections, bovine-human at 1 degree rests on one draw (`rt-mix-hb-k25-d00-w1`, k = 25 of
50) that v6 finds and v7, v8 and v9 miss. v6's detector splits it into 22 minority crystals, 20 of the
majority and three groups of two to four; v7-v9 merge the small groups into the largest population
before the take-out (`apply_min_group_size`), which then holds more than half of the data sets and so
becomes the background, and the minority is lost. `groups.min_group_size` is documented as "groups
smaller than this are not taken out of the OPTICS clusters"; v10 does that: which groups are taken out
is decided on all of the detector's groups (`take_out_groups`, as v6), and then a group smaller than
five that `small_group_admitted` does not admit keeps its OPTICS labels.

**Frozen procedure v10** (tag `frozen-v10`, `dials-patch/0001-0014`, `cf/frozen/groups_v10.py`,
`dials-patch/frozen.phil`, the same parameters as v9; v9 in `frozen-v9.phil` with `PATCH_MAX=13`). In-sample,
derived exactly from v6's groups (development notes): 252 of 720 cryo and 797 of 1548 room-temperature
minorities (v9: 251 and 796), the clusters of 17 trials changed and no outcome lost, no
single-population trial with a group where today's DIALS has none. Across the collections its k90 is
the data's in every cell but human-porcine at 10 degrees on room temperature (12 against 8: one draw
each at k = 8 and 10 that v6's detector itself rejects, separation 1.9 and 2.7 times the rest's
structure against its 3.5; under amendment 13's rule the data's k90 there is 12).

**Judgement.** v10 is the change as proposed, judged with v9, v8, v7 and v6 on the third collection by
the criteria of amendment 19. If v10 raises the false-group rate there, it fails its judgement, which is
reported, and the change is proposed with the latest of v9, v8 and v7 that does not, stated as chosen
after the judgement.

## Amendment 24 (limits: comparators of every size; after v10 was frozen, before any trial on the third collection)

Found while examining v10's one room-temperature cell above the data's limit across the collections
(human-porcine at 10 degrees): the single-population comparators of the structure test were computed
only for the design's minority sizes (room temperature 2-6, 8, 10, 12 and 25), while the sub-groups of a
minority recorded since amendment 14 take every size from max(2, k/2) up. A sub-group of 7, 9 or 11
crystals could never be matched by structure inside one species and counted as separable untested;
the two draws that keep v10 above the data's limit there are separable only through sub-groups of 7
and 9. The comparators are now computed for every size from 2 to n/2 (`cf/limits3.comparator_sizes`).
On room temperature this adds sizes 7, 9 and 11 to the 59 single-population trials of 24-25 crystals;
cryo already had every size. The records already written are supplemented by the same computation
(`python -m cf.limits3 supplement`, checked to reproduce the stored comparators of the old sizes), and
every classification is recomputed for every procedure. The change can only move a minority from
separable to a limit of the data; the classification with the design's sizes only is kept
(`classify(design_sizes_only=True)`) and how many minorities each procedure's classification moves is
reported (`REPORT.md`).

## Amendment 25 (the third collection's labels; before any trial on it)

The label verification of amendment 19 failed for all 16 crystals in the same direction: every
crystal of the first batch (image_140500 ... image_140818, taken as bovine) correlates more with the
labelled room-temperature human crystals (0.976-0.980 against 0.972-0.977), every one of the second
(image_143421 ... image_143867, taken as human) more with the bovine ones (0.983-0.985 against
0.957-0.960; DIALS's pairwise correlations on the first 20 degrees, `results/joint/rtbh_verify.json`).
The supporting information gives the species of data sets 0-7 and 8-15 but not which image numbers
they are; their order was an assumption of amendment 19. Three measurements settle it, none of them
a procedure under test: (1) the two batches are two populations by their own correlations (within
0.9947 and 0.9933, between 0.9770; the labelled bovine and human crystals: within 0.988 and 0.985,
between 0.970); (2) the second batch correlates with the labelled bovine crystals as they do with each
other (0.9845 against 0.9883) and far less with the human ones (0.9581), the first batch more with the
human (0.9784) than with the bovine (0.9743); (3) the first batch diffracts more strongly (median 3960
strong spots against 3153, every crystal of the second batch at most 4362), as the supporting
information's human crystals, larger (66 x 67 x 35 against 44 x 43 x 22 cubic micrometres) in the same
10 x 10 micrometre beam, would. The labels are therefore taken in the other order: the first batch
human, the second bovine (`pipeline/manifest.py`), the design regenerated and the verification rerun
with them before any trial; a crystal that still disagrees is left out.

## Amendment 26 (v11 frozen; before any trial on the third collection)

With amendment 24's comparators and v8's and v9's room-temperature runs complete (`cf.newproc`, both
false-group-safe there, so their finds count as separable), v10's k90 is above the data's across the
collections in three room-temperature cells (bovine-human at 10 and 5 degrees, bovine-porcine at 10):
pairs and triples precisely measured and strong (v8's condition) but under 30 standard errors (v9's),
which v8 finds and v10 does not. A group admitted by either condition: **Frozen procedure v11** (tag
`frozen-v11`, `dials-patch/0001-0015`, `cf/frozen/groups_v11.py`, `dials-patch/frozen.phil`; v10 in
`frozen-v10.phil` with `PATCH_MAX=14`): v10, and a group smaller than five also admitted when it meets
v8's condition (`groups.small_groups.precise`: standard error at most 0.002, two standard errors above
0.04 x 2/k, twice the rest's structure; `admitted_by_any`). In-sample, derived exactly from v6's groups
(`cf/devrules.py`, `results/dev/smallrules.json`): 252 of 720 cryo and 811 of 1548 room-temperature
minorities, no single-population trial with a group where today's DIALS has none; across the
collections its k90 is the data's in every cell of both collections. The within-species group nearest
its admission is v8's (7% short). **Judgement**: v11 is the change as proposed, judged with v10, v9, v8,
v7 and v6 on the third collection by the criteria of amendment 19; if it raises the false-group rate
there it fails, and the change is proposed with the latest of v10, v9, v8 and v7 that does not.
