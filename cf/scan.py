"""Fast greedy group scan on a symmetric pair matrix (vectorised over all seed pairs).

For every pair (a, b) a group is grown by repeatedly adding the crystal with the largest
summed value to the current members. Along every path, for every size k, the scan records a
score from the group's within-pair mean and the best outsider's mean to the group, and keeps
the best group per size and overall. Missing pairs must be given as 0 with weight 0.
"""

from __future__ import annotations

import numpy as np


def grow_all(x: np.ndarray, kmax: int, seeds: np.ndarray | None = None):
    """Greedy growth of all seeds. Returns (members order [s, kmax], within-sum [s, kmax+1],
    best-outsider mean [s, kmax+1]); index k of the last two arrays is for size k."""
    n = x.shape[0]
    if seeds is None:
        seeds = np.array(np.triu_indices(n, 1)).T
    s = len(seeds)
    order = np.zeros((s, kmax), dtype=np.int64)
    order[:, 0], order[:, 1] = seeds[:, 0], seeds[:, 1]
    inset = np.zeros((s, n), bool)
    inset[np.arange(s), seeds[:, 0]] = True
    inset[np.arange(s), seeds[:, 1]] = True
    tot = x[seeds[:, 0]] + x[seeds[:, 1]]  # summed value of each crystal to the members
    within = np.full((s, kmax + 1), np.nan)
    outside = np.full((s, kmax + 1), np.nan)
    ssum = x[seeds[:, 0], seeds[:, 1]].copy()
    for k in range(2, kmax + 1):
        within[:, k] = ssum / (k * (k - 1) / 2)
        cand = np.where(inset, -np.inf, tot)
        best = cand.argmax(1)
        outside[:, k] = cand[np.arange(s), best] / k
        if k == kmax:
            break
        ssum = ssum + tot[np.arange(s), best]
        inset[np.arange(s), best] = True
        tot = tot + x[best]
        order[:, k] = best
    return order, within, outside
