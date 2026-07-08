# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Tests for LP parsing, mat_idx correctness, M+ lifting matrix dimensions,
# and SDPModel export correctness (to_sdpnal, to_adal).

import os
import tempfile
import itertools
import numpy as np
import pytest
from scipy import sparse

from pyModules.SDPLifting import (
    mat_idx,
    read_constr_from_lp,
    _compute_m_plus_model,
)
from pyModules.SDPModel import SDPModel


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

    def test_missing_file_raises(self, tmp_path):
        missing = str(tmp_path / "does_not_exist.lp")
        with pytest.raises(FileNotFoundError):
            list(read_constr_from_lp(missing))


# ---------------------------------------------------------------------------
# M+ lifting — SDPModel structure
# ---------------------------------------------------------------------------

class TestMPlusModel:

    def test_returns_sdpmodel(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        assert isinstance(model, SDPModel)

    def test_matrix_variable_dimension(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        n = model.n
        assert n == 3
        assert model.dim == n + 1
        assert model.A.shape[1] == model.packed_dim

    def test_equality_count(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        n = model.n
        eq_rows = int((model.row_senses == '=').sum())
        # X_00 = 1  plus  X_ii = X_0i for each node
        assert eq_rows == n + 1

    def test_cut_classes_length_matches_inequality_rows(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        ineq_rows = int((model.row_senses == '>').sum())
        assert len(model.cut_classes) == ineq_rows

    def test_cut_classes_values(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        assert set(model.cut_classes).issubset({1., 2., 3., 4.})

    def test_write_sdpnal_mat(self, triangle_lp_path, tmp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        out_path = str(tmp_path / "triangle_edge")
        model.write_sdpnal_mat(out_path, do_split=False)
        assert os.path.exists(out_path + '.mat')

    def test_write_sdpnal_mat_split(self, triangle_lp_path, tmp_path):
        model = _compute_m_plus_model(triangle_lp_path)
        out_path = str(tmp_path / "triangle_edge")
        model.write_sdpnal_mat(out_path, do_split=True)
        for i in range(1, 6):
            assert os.path.exists(out_path + '_%d.mat' % i)


# ---------------------------------------------------------------------------
# SDPModel export correctness
# ---------------------------------------------------------------------------

class TestSDPModelExports:

    @pytest.fixture
    def triangle_model(self, triangle_lp_path):
        return _compute_m_plus_model(triangle_lp_path)

    def test_to_sdpnal_off_diagonal_scaling(self, triangle_model):
        """Off-diagonal entries of At must be sqrt(2) times the canonical value."""
        model = triangle_model
        s = model._svec_scale()
        sdpnal = model.to_sdpnal()
        At = sdpnal['At']         # packed_dim x m_eq
        A_eq = model.A[model.row_senses == '=']  # m_eq x packed_dim canonical

        # At = diag(s) @ A_eq.T  =>  At.T = A_eq * s
        At_T_dense = np.asarray(At.T.todense())
        A_eq_dense = np.asarray(A_eq.todense())
        np.testing.assert_allclose(At_T_dense, A_eq_dense * s, atol=1e-12)

    def test_to_adal_full_matrix_symmetry(self, triangle_model):
        """The ADAL constraint matrix must produce a symmetric full matrix for each row."""
        model = triangle_model
        adal = model.to_adal()
        A_full = adal['A']   # (num_rows x dim^2)
        dim = model.dim
        for row_idx in range(min(5, A_full.shape[0])):
            row = np.asarray(A_full[row_idx].todense()).ravel()
            M = row.reshape(dim, dim)
            np.testing.assert_allclose(M, M.T, atol=1e-12)

    def test_inner_product_round_trip_sdpnal(self, triangle_model):
        """trace(A X) via canonical inner product == dot(svec(A), svec(X)) in SDPNAL format."""
        model = triangle_model
        dim = model.dim
        rng = np.random.default_rng(0)
        # Random symmetric positive-ish matrix X (doesn't need to be PSD for this check)
        X_lower = rng.standard_normal((dim, dim))
        X = X_lower + X_lower.T
        # Packed svec of X: diagonal entries as-is, off-diagonal * sqrt(2)
        s = model._svec_scale()
        x_packed = np.zeros(model.packed_dim)
        for i in range(dim):
            for j in range(i + 1):
                k = i * (i + 1) // 2 + j
                x_packed[k] = X[i, j] * s[k]  # svec scaling applied to X

        sdpnal = model.to_sdpnal()
        At = sdpnal['At']  # packed_dim x m_eq

        # canonical inner products: <A_i, X> = sum_j A_ii X_ii + 2*sum_{i>j} A_ij X_ij
        A_eq = model.A[model.row_senses == '='].toarray()  # (m_eq x packed_dim) canonical
        pi, pj = model._packed_ij()
        # Build the full inner-product by hand
        trace_results = []
        for row in range(A_eq.shape[0]):
            val = 0.0
            for k in range(model.packed_dim):
                i_k, j_k = pi[k], pj[k]
                if i_k == j_k:
                    val += A_eq[row, k] * X[i_k, j_k]
                else:
                    val += 2 * A_eq[row, k] * X[i_k, j_k]
            trace_results.append(val)

        # SDPNAL inner products: At.T @ x_packed (for each equality constraint)
        sdpnal_results = np.asarray(At.T @ x_packed).ravel()
        np.testing.assert_allclose(sdpnal_results, trace_results, atol=1e-10)

    def test_to_mosek_inner_product(self, triangle_model):
        """trace(A X) via canonical inner product == MOSEK triplet dot X for each row."""
        model = triangle_model
        dim = model.dim
        rng = np.random.default_rng(2)
        X_lower = rng.standard_normal((dim, dim))
        X = X_lower + X_lower.T   # symmetric; need not be PSD for this check

        d = model.to_mosek()
        pi, pj = model._packed_ij()

        # Canonical inner product for each constraint row
        A_arr = model.A.toarray()   # (num_rows x packed_dim)
        trace_canonical = []
        for row in range(model.num_rows):
            val = 0.0
            for k in range(model.packed_dim):
                i_k, j_k = pi[k], pj[k]
                val += (A_arr[row, k] * X[i_k, j_k]
                        if i_k == j_k
                        else 2 * A_arr[row, k] * X[i_k, j_k])
            trace_canonical.append(val)

        # MOSEK triplet inner product: same formula via (cons, row, col, val) tuples
        trace_mosek = np.zeros(model.num_rows)
        for con, r, c, v in zip(d['A_cons'], d['A_rows'], d['A_cols'], d['A_val']):
            trace_mosek[con] += v * X[r, c] if r == c else 2 * v * X[r, c]

        np.testing.assert_allclose(trace_mosek, trace_canonical, atol=1e-12)

    def test_to_mosek_bounds(self, triangle_model):
        """blc/buc must encode the constraint senses correctly."""
        model = triangle_model
        d = model.to_mosek()
        for i, sense in enumerate(model.row_senses):
            rhs = model.row_rhs[i]
            if sense == '=':
                assert d['blc'][i] == rhs and d['buc'][i] == rhs
            elif sense == '>':
                assert d['blc'][i] == rhs and d['buc'][i] == np.inf
            else:
                assert d['blc'][i] == -np.inf and d['buc'][i] == rhs

    def test_inner_product_round_trip_adal(self, triangle_model):
        """trace(A X) via canonical inner product == dot(A_full_row, X_flat) in ADAL format."""
        model = triangle_model
        dim = model.dim
        rng = np.random.default_rng(1)
        X_lower = rng.standard_normal((dim, dim))
        X = X_lower + X_lower.T
        X_flat = X.ravel()

        adal = model.to_adal()
        A_full = adal['A']   # (num_rows x dim^2)

        # canonical inner products for equality rows
        A_eq = model.A[model.row_senses == '='].toarray()
        pi, pj = model._packed_ij()

        trace_results = []
        for row in range(A_eq.shape[0]):
            val = 0.0
            for k in range(model.packed_dim):
                i_k, j_k = pi[k], pj[k]
                if i_k == j_k:
                    val += A_eq[row, k] * X[i_k, j_k]
                else:
                    val += 2 * A_eq[row, k] * X[i_k, j_k]
            trace_results.append(val)

        # In to_adal, inequalities come first; equality rows follow
        n_ineq = int((model.row_senses == '>').sum())
        adal_results = np.asarray(A_full[n_ineq:] @ X_flat).ravel()
        np.testing.assert_allclose(adal_results, trace_results, atol=1e-10)
