# Goal

Establish how small a population of crystals automatic multi-crystal clustering can be trusted to find, why today's DIALS loses smaller groups, and the best procedure that could exist for finding them without inventing groups. Aspire to the best that could exist, not a good one. Today's DIALS is a floor, not a target.

## Specific

1. **The limit of today's DIALS** on the public labelled data: the smallest minority it reports as a group, against how different it is and how much data each crystal carries, and how often it reports groups in single-population data.
2. **The cause,** shown by evidence rather than inferred.
3. **A better procedure,** as a change DIALS's developers could merge, with tests.
4. **`FINDINGS.md`** (each claim with its evidence), **`REPORT.md`** (each criterion before and after, what failed, what remains open), and one command that regenerates every number.

## Measurable

Every number is checked against references this work does not produce.

- **Fidelity:** today's DIALS reproduces the paper's three pure clusters on the 36 cryogenic crystals, or the difference is explained.
- **Detection:** the probability that a minority is found, with intervals over independent draws, and the limit it implies, under one definition of "found" fixed before measuring.
- **False groups:** the rate at which groups are reported in single-population data, which a better procedure must not raise.
- **The frontier:** the non-dominated set of detection limit, false-group rate and cost across every procedure tried. After every success, push it on whichever criterion is furthest from its limit.
- **Limits:** only the data themselves: a minority shown, by the measurement that settles it, to lie within the noise for any procedure. A default, parameter or sample size chosen here and anything left untried are not limits.
- **Evidence:** every number comes from one named commit and says what it is compared against; no claim is stronger than its data (a group seen in a plot is not a group found); a procedure is judged with one configuration on data it was not tuned on.

## Relevant

Crystallographers screening fragments across many crystals need to know whether a handful of ligand-bound or differently packed crystals will appear in the automatic output or must be hunted by hand.

## Done

Done means `COMPLETE.md` shows, through the measurements above, that all four deliverables hold, today's DIALS reproduces the published separation, the smallest detectable group and its cause are measured for today's DIALS and the new procedure on both public collections, and every remaining limit belongs to the problem itself. There is no deadline and no other way to finish.

## Rules

- Every choice is yours.
- Never invent results or references; label estimates as estimates.
- Contact nobody and publish nothing: no email, posts, issue comments, pushes or pull requests.
- Work only inside this directory; do not read or search local files elsewhere. The internet is yours.
- Commit as you go, so the repository alone can resume the work. Do not edit PROBLEM.md or GOAL.md.
- Keep data and environments small, and choose compute resources appropriate to each stage.
