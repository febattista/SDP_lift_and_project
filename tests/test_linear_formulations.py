# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for nodal coefficient computation, caching and validation."""

import json
import os

import networkx as nx
import pytest

from pyModules.LinearFormulations import (
    _CLIQUER, check_nodal_coefficients, compute_alpha, compute_theta,
)

needs_cliquer = pytest.mark.skipif(
    not os.path.exists(_CLIQUER),
    reason='cliquer binary not built (run install.py)')


class TestCheckNodalCoefficients:

    def test_valid_passes(self):
        check_nodal_coefficients({0: 3.0, 1: 2.0}, {0: 3, 1: 1}, 'g')

    def test_violation_raises_and_names_node(self):
        with pytest.raises(ValueError, match='node 1.*theta=1.*alpha=2'):
            check_nodal_coefficients({0: 3.0, 1: 1.0}, {0: 3, 1: 2}, 'g')

    def test_disjoint_keys_ignored(self):
        check_nodal_coefficients({0: 1.0}, {1: 5}, 'g')


class TestCoefficientCaches:
    """Cache semantics of compute_alpha / compute_theta.

    On C5 every neighborhood is two non-adjacent nodes: alpha = theta = 2.
    """

    @pytest.fixture
    def C5(self):
        return nx.cycle_graph(5)

    def _seed_bad_cache(self, path):
        with open(path, 'w') as f:
            json.dump({str(v): [99.0, 0.0] for v in range(5)}, f)

    @needs_cliquer
    def test_alpha_reads_cache_by_default(self, C5, tmp_path):
        self._seed_bad_cache(tmp_path / 'C5_alpha.json')
        assert compute_alpha(C5, 'C5', model_out_dir=str(tmp_path)) == \
            {v: 99.0 for v in range(5)}

    @needs_cliquer
    def test_alpha_refresh_recomputes_and_overwrites(self, C5, tmp_path):
        self._seed_bad_cache(tmp_path / 'C5_alpha.json')
        alpha = compute_alpha(C5, 'C5', model_out_dir=str(tmp_path),
                              refresh=True)
        assert alpha == {v: 2 for v in range(5)}
        with open(tmp_path / 'C5_alpha.json') as f:
            assert {int(k): v[0] for k, v in json.load(f).items()} == alpha

    def test_theta_reads_cache_by_default(self, C5, tmp_path):
        self._seed_bad_cache(tmp_path / 'C5_theta.json')
        assert compute_theta(C5, 'C5', model_out_dir=str(tmp_path)) == \
            {v: 99.0 for v in range(5)}

    def test_theta_refresh_recomputes_and_overwrites(self, C5, tmp_path):
        self._seed_bad_cache(tmp_path / 'C5_theta.json')
        theta = compute_theta(C5, 'C5', model_out_dir=str(tmp_path),
                              refresh=True)
        # each N(v) is edgeless on 2 nodes: theta = 2 exactly — this is the
        # integral case the safe-bound/rounding fix must get right
        assert theta == {v: 2.0 for v in range(5)}
