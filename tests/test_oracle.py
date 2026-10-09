"""Synthetic checks of the oracle statistic, its p-value and the decoy search."""

import numpy as np

from cf.oracle import decoy_max, deviation_correlations, group_means, oracle_p


def synthetic(n=14, nu=3000, k=2, delta=0.0, outlier=None, seed=1):
    rng = np.random.default_rng(seed)
    mu = rng.exponential(1.0, nu)
    sig = 0.3 * np.ones((n, nu))
    y = mu[None, :] + sig * rng.standard_normal((n, nu))
    d = delta * rng.standard_normal(nu)
    y[:k] += d[None, :]
    if outlier is not None:  # a crystal with no signal at all
        y[outlier] = 5 * rng.standard_normal(nu)
    m = np.ones((n, nu), bool)
    return y, sig**2, m


def test_null_is_uniformish():
    ps = []
    for s in range(30):
        y, v, m = synthetic(delta=0.0, seed=s)
        c, _ = deviation_correlations(y, v, m)
        ps.append(oracle_p(c, [0, 1], list(range(2, 14)))["p"])
    assert 0.25 < np.mean(ps) < 0.75


def test_planted_pair_found():
    y, v, m = synthetic(delta=0.15, seed=3)
    c, _ = deviation_correlations(y, v, m)
    assert oracle_p(c, [0, 1], list(range(2, 14)))["p"] < 0.02


def test_outlier_does_not_inflate_null():
    y, v, m = synthetic(delta=0.15, outlier=5, seed=4)
    c, _ = deviation_correlations(y, v, m)
    assert oracle_p(c, [0, 1], list(range(2, 14)))["p"] < 0.02
    assert abs(c[5, 6]) < 0.1


def test_decoy_local_search_matches_exhaustive():
    rng = np.random.default_rng(0)
    a = rng.standard_normal((12, 12))
    c = (a + a.T) / 2
    np.fill_diagonal(c, 0)
    exact = decoy_max(c, list(range(12)), 4)
    approx = decoy_max(c, list(range(12)), 4, exhaustive_limit=0)
    assert exact["exact"] and not approx["exact"]
    assert approx["max"] <= exact["max"] + 1e-12
    assert approx["max"] >= exact["max"] - 0.05
    assert np.isclose(group_means(c, np.array([exact["members"]]))[0], exact["max"])


def test_rank1_residual_oracle():
    from cf.oracle import decoy_max_resid, oracle_p_resid, residuals

    rng = np.random.default_rng(5)
    n = 14
    rho = rng.uniform(0.8, 0.99, n)
    r = np.outer(rho, rho)
    r[np.ix_([0, 1], [0, 1])] = rho[0] * rho[1] * 1.0
    r[0, 2:] *= 0.95  # crystals 0 and 1 differ from the rest
    r[1, 2:] *= 0.95
    r[2:, 0] = r[0, 2:]
    r[2:, 1] = r[1, 2:]
    r = r + 0.002 * rng.standard_normal((n, n))
    r = (r + r.T) / 2
    w = np.ones((n, n))
    np.fill_diagonal(w, 0)
    e, _ = residuals(r, w)
    o = oracle_p_resid(e, [0, 1], list(range(2, n)))
    assert o["p"] < 0.02
    assert decoy_max_resid(e, list(range(2, n)), 2)["max"] < o["S"]
