import numpy as np

from cf.bootstrap import weighted_cc_matrix_w
from cf.fastcc import weighted_cc_matrix


def test_unit_weights_match_fastcc():
    rng = np.random.default_rng(0)
    n, nu = 6, 400
    mu = rng.exponential(1.0, nu)
    y = mu[None, :] + 0.2 * rng.standard_normal((n, nu))
    v = (0.2 + 0.1 * rng.random((n, nu))) ** 2
    m = rng.random((n, nu)) < 0.7
    r1, _, _ = weighted_cc_matrix(y, v, m)
    r2 = weighted_cc_matrix_w(y, v, m, np.ones(nu))
    assert np.allclose(r1, r2)
