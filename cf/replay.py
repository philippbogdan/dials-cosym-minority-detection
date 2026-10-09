"""Replay today's DIALS clustering from a trial's saved rij and wij matrices, calling DIALS's
own functions (CosymAnalysis._determine_dimensions, _optimise, CorrelationMatrix.
optimise_clustering), so that single stages can be changed one at a time (cause ablations)
without re-scaling. `replay(...)` with no changes reproduces dials.correlation_matrix's
significant clusters for the same rij/wij (checked in tests/test_replay.py).

Changes (keyword arguments):
  dims            force the number of dimensions (None: DIALS's automatic choice)
  min_samples_floor  lower bound of the OPTICS min_samples search (DIALS: 5)
  heuristic_only  use the paper's heuristic min_samples, no optuna search
  unit_vectors    run OPTICS on coordinates scaled to unit length (cosine geometry)
"""

from __future__ import annotations

import copy
import random
import warnings

import numpy as np
from ordered_set import OrderedSet
from sklearn.cluster import OPTICS

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.algorithms.symmetry.cosym import CosymAnalysis
from dials.algorithms.symmetry.cosym.target import Target
from dials.array_family import flex
from dials.util.observer import Subject


class SavedTarget(Target):
    """A cosym Target built from saved rij/wij instead of from intensities."""

    def __init__(self, rij, wij):  # noqa: super().__init__ deliberately not called
        n = rij.shape[0]
        self.rij_matrix = np.array(rij, dtype=np.float64)
        np.fill_diagonal(self.rij_matrix, 0.0)  # cosym excludes self-correlations
        self.wij_matrix = np.array(wij, dtype=np.float64)
        self._lattices = np.arange(n)
        self.sym_ops = OrderedSet(["x,y,z"])
        self._weights = "count"
        self.dim = 2


def _analysis(rij, wij, params):
    ca = object.__new__(CosymAnalysis)
    Subject.__init__(ca, events=["optimised", "analysed_symmetry", "analysed_clusters"])
    ca.params = params
    ca.target = SavedTarget(rij, wij)
    ca.dataset_ids = flex.int(range(rij.shape[0]))
    return ca


def coordinates(rij, wij, params, dims=None):
    """DIALS's dimension choice and cosym coordinates for saved rij/wij."""
    ca = _analysis(rij, wij, params)
    n = rij.shape[0]
    dims_to_test = min(n, params.dimensionality_assessment.maximum_dimensions)
    curve = None
    if dims is None:
        _, curve = ca._determine_dimensions(
            dims_to_test,
            outlier_rejection=params.dimensionality_assessment.outlier_rejection,
        )
    else:
        ca.target.set_dimensions(dims)
    ca._optimise(
        params.minimization.engine,
        max_iterations=params.minimization.max_iterations,
        max_calls=params.minimization.max_calls,
    )
    return np.array(ca.coords), int(ca.target.dim), curve


def optics_labels(coords, dim, params, min_samples_floor=5, heuristic_only=False,
                  unit_vectors=False):
    """DIALS's cluster_cosine_coords on given coordinates, optionally changed."""
    n = coords.shape[0]
    x = coords
    if unit_vectors:
        x = coords / np.linalg.norm(coords, axis=1, keepdims=True)
    sc = params.significant_clusters
    initial = max(min_samples_floor, int((n / dim) * sc.min_points_buffer))
    if heuristic_only:
        model = OPTICS(min_samples=min(initial, n), xi=sc.xi, max_eps=sc.max_distance)
        model.fit(x)
        return np.array(model.labels_), initial
    if min_samples_floor == 5:
        fn = CorrelationMatrix.optimise_clustering
    else:
        fn = _optimise_clustering_with_floor(min_samples_floor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ms, _, labels, _ = fn(x, list(range(n)), initial, sc.xi, sc.max_distance,
                              sc.noise_penalty.alpha, sc.noise_penalty.gamma)
    return np.array(labels), ms


def _optimise_clustering_with_floor(floor):
    """CorrelationMatrix.optimise_clustering with the min_samples search starting at
    `floor` instead of 5 (its source is otherwise unchanged)."""
    import inspect
    import textwrap

    src = textwrap.dedent(inspect.getsource(CorrelationMatrix.optimise_clustering))
    src = src.replace("@staticmethod\n", "")
    import re

    pattern = r'trial\.suggest_int\(\s*"min_samples",\s*5,'
    assert re.search(pattern, src), "DIALS source changed; update the ablation"
    src = re.sub(pattern, f'trial.suggest_int("min_samples", {floor},', src)
    import dials.algorithms.correlation.analysis as mod

    ns = dict(mod.__dict__)
    exec(src, ns)
    return ns["optimise_clustering"]


def replay(rij, wij, params=None, dims=None, seed=True, **changes):
    if params is None:
        from cf.today import default_params

        params = default_params()
    params = copy.deepcopy(params)
    if seed and params.seed is not None:
        flex.set_random_seed(params.seed)
        np.random.seed(params.seed)
        random.seed(params.seed)
    coords, dim, curve = coordinates(rij, wij, params, dims=dims)
    labels, ms = optics_labels(coords, dim, params, **changes)
    return dict(labels=labels, dims=dim, min_samples=ms, coords=coords, curve=curve)
