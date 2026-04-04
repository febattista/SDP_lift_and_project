# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Tests for LP parsing, mat_idx correctness, and M+ lifting matrix dimensions.

import os
import tempfile
import itertools
import numpy as np
import pytest

from pyModules.SDPLifting import (
    mat_idx,
    read_constr_from_lp,
    _compute_m_plus_model,
    write_sdpnal_mat,
)


# ---------------------------------------------------------------------------
# mat_idx
# ---------------------------------------------------------------------------

class TestMatIdx:

    def test_diagonal_entry(self):
        # mat_idx(i, i) should map to the i-th triangular number
        assert mat_idx(0, 0) == 0
        assert mat_idx(1, 1) == 2
        assert mat_idx(2, 2) == 5

    def test_symmetry(self):
        for i in range(5):
            for j in range(5):
                assert mat_idx(i, j) == mat_idx(j, i)

    def test_returns_int(self):
        idx = mat_idx(3, 2)
        assert isinstance(idx, int)

    def test_no_collision_for_small_matrix(self):
        n = 6
        indices = set()
        for i in range(n):
            for j in range(i + 1):
                indices.add(mat_idx(i, j))
        # All indices should be distinct
        assert len(indices) == (n * (n + 1)) // 2


# ---------------------------------------------------------------------------
# LP parsing
# ---------------------------------------------------------------------------

# Minimal LP that models a triangle (3 nodes, 3 edges):
#   max  x[1] + x[2] + x[3]
#   s.t. x[1] + x[2] <= 1
#        x[1] + x[3] <= 1
#        x[2] + x[3] <= 1
#   bounds 0 <= x[i] <= 1
_TRIANGLE_LP = """\
obj: x[1] + x[2] + x[3]
subject to
 edge1: x[1] + x[2] <= 1
 edge2: x[1] + x[3] <= 1
 edge3: x[2] + x[3] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
 0 <= x[3] <= 1
end
"""


@pytest.fixture
def triangle_lp_path(tmp_path):
    path = tmp_path / "triangle.lp"
    path.write_text(_TRIANGLE_LP)
    return str(path)


class TestReadConstrFromLp:

    def test_first_item_is_objective(self, triangle_lp_path):
        constrs = list(read_constr_from_lp(triangle_lp_path))
        obj = constrs[0]
        # Objective has 3 variables each with coefficient 1.0
        assert len(obj) == 3
        assert all(v == 1.0 for v in obj.values())

    def test_constraint_count(self, triangle_lp_path):
        constrs = list(read_constr_from_lp(triangle_lp_path))
        # 1 objective + 3 constraints
        assert len(constrs) == 4

    def test_constraint_rhs(self, triangle_lp_path):
        constrs = list(read_constr_from_lp(triangle_lp_path))
        for constr, b in constrs[1:]:
            assert b == 1.0

    def test_constraint_variables(self, triangle_lp_path):
        constrs = list(read_constr_from_lp(triangle_lp_path))
        for constr, b in constrs[1:]:
            assert len(constr) == 2
            assert all(v == 1.0 for v in constr.values())


# ---------------------------------------------------------------------------
# M+ lifting — matrix dimension checks
# ---------------------------------------------------------------------------

class TestMPlusModel:

    def test_matrix_variable_dimension(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        n = model['n']
        assert n == 3
        # Packed dimension of (n+1) x (n+1) symmetric matrix
        expected_dim = ((n + 1) * (n + 2)) // 2
        assert model['At'].shape[0] == expected_dim
        assert model['Bt'].shape[0] == expected_dim

    def test_equality_count(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        n = model['n']
        # Equalities: X_00 = 1 (1) + X_ii = X_0i for each node (n)
        assert model['At'].shape[1] == n + 1

    def test_cut_classes_length_matches_Bt(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        assert len(model['cut_classes']) == model['Bt'].shape[1]

    def test_cut_classes_values(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        # Classes are 1, 2, 3, or 4
        assert set(model['cut_classes']).issubset({1., 2., 3., 4.})

    def test_write_sdpnal_mat(self, triangle_lp_path, tmp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        out_path = str(tmp_path / "triangle_edge")
        write_sdpnal_mat(model, out_path, do_split=False)
        assert os.path.exists(out_path + '.mat')

    def test_write_sdpnal_mat_split(self, triangle_lp_path, tmp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        out_path = str(tmp_path / "triangle_edge")
        write_sdpnal_mat(model, out_path, do_split=True)
        for i in range(1, 6):
            assert os.path.exists(out_path + '_%d.mat' % i)
