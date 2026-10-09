# The proposed DIALS change

Against DIALS v3.30.1 (commit 17cfc5a of https://github.com/dials/dials).

- `squashed/correlation-find-small-groups.patch`: the change as one patch (0001-0015):
  `dials/algorithms/correlation/groups.py` (new), a `significant_clusters.find_small_groups` option
  and its `groups` scope in `dials/algorithms/correlation/analysis.py` (used by
  `dials.correlation_matrix` and `xia2.multiplex`), tests in
  `tests/algorithms/correlation/test_groups.py`. It passes DIALS's pre-commit checks (ruff check and
  ruff format, the version DIALS pins).
- `0001-0015`: the same change as it was developed, one commit per frozen version
  (0001 = v1, +0002-0003 = v2, +0004 = v3, +0005 = v4, +0006 = v5, +0007 = v6), then:
  0008 fixes a failure found by running the patched DIALS end to end (DIALS logged OPTICS's
  silhouette as a number whenever there was more than one cluster, but when OPTICS finds one it is a
  message, and the change can add clusters to that one; 25 of 72 sampled trials failed), with an
  end-to-end test of that case; 0009 computes the bootstrap replicates by weighted sums, the same
  matrices to rounding and the same groups on every cryo trial, 10-15 times faster
  (`cf/devfast.py`); neither changes v6's decisions. 0010 adds `groups.min_group_size` (default 5,
  OPTICS's own smallest cluster): groups smaller than that are not taken out (v7, amendment 17,
  after v6 raised the false-group rate on room temperature, where crystals of one species hold
  sub-groups of two to four; `min_group_size = 2` gives v6). 0011 adds `groups.small_groups`: a
  smaller group is kept when its coherence is precisely measured (v8, amendment 20). 0012 makes that
  condition one of signal to noise, `small_group_admitted` (v9, amendment 22: the coherence 30
  standard errors out of the counting noise); 0013 formats the change with ruff (no change in
  behaviour); 0014 decides which groups are taken out on all of them before the minimum group size
  (v10, amendment 23: a group too small is no longer merged into the largest population first, which
  could make a minority the background population); 0015 admits a small group by either condition,
  far out of its noise (0012) or precisely measured and strong (0011's, `small_groups.precise`; v11,
  amendment 26).
- `frozen.phil`: the parameters of the change as proposed (v11); `frozen-v10.phil` (with 0001-0014,
  `PATCH_MAX=14`), `frozen-v9.phil` (0001-0013,
  `PATCH_MAX=13`), `frozen-v8.phil` (0001-0011, `PATCH_MAX=11`), `frozen-v7.phil` (0001-0010), `frozen-v6.phil` (0001-0007): the configurations
  judged with it; `frozen-v1.phil`, `frozen-v3.phil`, `frozen-v5.phil`: earlier versions.

Apply with `git apply` from the DIALS source root, or as an overlay on an installed DIALS 3.30.1
(`pipeline/make_overlay.sh`; `PATCH_MAX` keeps the patches up to a number); run its tests with
`pipeline/test_patch.sh`.
