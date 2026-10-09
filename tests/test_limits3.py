import math

import numpy as np

from cf.limits3 import Trial, mixed_candidates, pair_factors, subgroup_candidates


def _trial(n=12, group=(0, 1, 2, 3), shift=0.05, seed=0):
    """A correlation matrix of one population with a coherent group, analytic-like noise."""
    rng = np.random.default_rng(seed)
    rho = np.full(n, 0.95)
    r = np.outer(rho, rho)
    for i in group:
        for j in group:
            if i != j:
                r[i, j] += shift
    se = 0.01
    noise = rng.normal(0, se, (n, n))
    r = np.clip(r + (noise + noise.T) / 2, -0.99, 0.99)
    np.fill_diagonal(r, 0.0)
    w = np.full((n, n), 50.0)
    np.fill_diagonal(w, 0.0)
    deltas = rng.normal(0, se, (60, n, n))
    deltas = (deltas + deltas.transpose(0, 2, 1)) / 2
    for b in range(len(deltas)):
        np.fill_diagonal(deltas[b], 0.0)
    return Trial(r, w, deltas, 1.0)


def test_pair_factors_are_geometric_means():
    f = pair_factors(["a", "b", "c"], {"a": "x", "b": "y", "c": "z"}, {"x": 0.64, "y": 1.44}, 1.0)
    assert np.allclose(f, [[0.64, 0.96, 0.8], [0.96, 1.44, 1.2], [0.8, 1.2, 1.0]])


def test_subgroups_are_qualifying_subsets_of_the_minority():
    T = _trial()
    M = np.array([0, 1, 2, 3, 4, 5])
    smin = int(math.ceil(max(2, len(M) / 2)))
    subs = subgroup_candidates(T, M, smin)
    assert subs
    for G in subs:
        assert set(G.tolist()) <= set(M.tolist())
        assert smin <= len(G) < len(M)


def test_mixed_groups_keep_the_minority_in_the_majority_of_the_group():
    T = _trial()
    M = np.array([0, 1, 2, 3])
    A = np.arange(4, 12)
    for G, n_min in mixed_candidates(T, [M], A):
        members = set(G.tolist())
        hit = len(members & set(M.tolist()))
        assert hit == n_min
        assert len(members) - hit < hit
