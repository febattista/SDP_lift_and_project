# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Tests for LP parsing, mat_idx correctness, M+ lifting matrix dimensions,
# and SDPModel export correctness (to_sdpnal, to_adal).

import os
import itertools
import numpy as np
import pytest

import networkx as nx

from pyModules.SDPLifting import (
    mat_idx,
    read_constr_from_lp,
    _compute_m_plus_model,
    lovasz_schrijver_filter,
    Theta_SDP,
    Theta_plus_SDP,
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
        ineq_rows = int((model.row_senses != '=').sum())
        assert len(model.cut_classes) == ineq_rows

    def test_lifted_rows_are_upper_bounds(self, triangle_lp_path):
        """Lifted inequality rows are <= constraints (see docs/plans/m-plus-kk-lifting.md)."""
        model = _compute_m_plus_model(triangle_lp_path)
        senses = set(model.row_senses)
        assert senses == {'=', '<'}

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

        # MOSEK triplet inner product: same formula via (cons, row, col, val)
        # tuples (d['num_rows'] > model.num_rows: L/U bounds become extra rows)
        trace_mosek = np.zeros(d['num_rows'])
        for con, r, c, v in zip(d['A_cons'], d['A_rows'], d['A_cols'], d['A_val']):
            trace_mosek[con] += v * X[r, c] if r == c else 2 * v * X[r, c]

        np.testing.assert_allclose(trace_mosek[:model.num_rows],
                                   trace_canonical, atol=1e-12)

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
        n_ineq = int((model.row_senses != '=').sum())
        adal_results = np.asarray(A_full[n_ineq:] @ X_flat).ravel()
        np.testing.assert_allclose(adal_results, trace_results, atol=1e-10)


# ---------------------------------------------------------------------------
# Elementwise bounds L <= X <= U
# ---------------------------------------------------------------------------
# SDPNAL+ takes L/U natively (write_sdpnal_mat); MOSEK bar-variables cannot
# carry entry bounds, so to_mosek converts finite bounds into explicit rows;
# to_lp maps them onto variable bounds; the internal ADMM solver already
# enforces X >= 0, so to_adal converts only bounds differing from that.

class TestElementwiseBounds:

    @pytest.fixture
    def triangle_model(self, triangle_lp_path):
        return _compute_m_plus_model(triangle_lp_path)

    def test_builders_set_L(self, triangle_lp_path):
        for mode in ('ls', 'kk'):
            model = _compute_m_plus_model(triangle_lp_path, lift_mode=mode)
            assert model.L == 0. and model.U is None
        G = nx.cycle_graph(5)
        assert Theta_plus_SDP(G).L == 0. and Theta_plus_SDP(G).U is None
        assert Theta_SDP(G).L is None and Theta_SDP(G).U is None

    def test_bounds_packed_scalar_and_matrix(self, triangle_model):
        lb, ub = triangle_model._bounds_packed()
        assert (lb == 0.).all() and (ub == np.inf).all()
        dim = triangle_model.dim
        L_mat = np.full((dim, dim), -np.inf)
        L_mat[1, 0] = L_mat[0, 1] = 0.25
        triangle_model.L = L_mat
        lb, _ = triangle_model._bounds_packed()
        assert lb[mat_idx(1, 0)] == 0.25
        assert np.isneginf(np.delete(lb, mat_idx(1, 0))).all()

    def test_to_mosek_bound_rows(self, triangle_model):
        model = triangle_model
        d = model.to_mosek()
        # L = 0 scalar: one X_ij >= 0 row per packed entry
        assert d['num_rows'] == model.num_rows + model.packed_dim
        assert (d['row_senses'][model.num_rows:] == '>').all()
        np.testing.assert_array_equal(d['blc'][model.num_rows:], 0.)
        assert (d['buc'][model.num_rows:] == np.inf).all()
        # each converted row's inner product is exactly X_ij
        rng = np.random.default_rng(4)
        X_half = rng.standard_normal((model.dim, model.dim))
        X = X_half + X_half.T
        vals = np.zeros(d['num_rows'])
        for con, r, c, v in zip(d['A_cons'], d['A_rows'],
                                d['A_cols'], d['A_val']):
            vals[con] += v * X[r, c] if r == c else 2 * v * X[r, c]
        pi, pj = model._packed_ij()
        np.testing.assert_allclose(vals[model.num_rows:], X[pi, pj],
                                   atol=1e-12)

    def test_to_mosek_no_bounds_no_extra_rows(self, triangle_model):
        triangle_model.L = None
        d = triangle_model.to_mosek()
        assert d['num_rows'] == triangle_model.num_rows

    def test_write_sdpnal_mat_bounds(self, triangle_model, tmp_path):
        from scipy.io import loadmat
        out = str(tmp_path / 'with_L')
        triangle_model.write_sdpnal_mat(out)
        m = loadmat(out + '.mat')
        assert m['L'].item() == 0.
        assert 'U' not in m
        # models without bounds write no L (Theta must not gain X >= 0)
        out2 = str(tmp_path / 'no_L')
        Theta_SDP(nx.cycle_graph(5)).write_sdpnal_mat(out2)
        assert 'L' not in loadmat(out2 + '.mat')

    def test_to_adal_builtin_nonneg_unconverted(self, triangle_model):
        """L = 0 matches ADMM_3b's built-in X >= 0: no rows are added."""
        d = triangle_model.to_adal()
        assert d['A'].shape[0] == triangle_model.num_rows

    def test_to_adal_other_bounds_become_rows(self, triangle_model):
        model = triangle_model
        model.U = 1.
        d = model.to_adal()
        # one '<' row per packed entry for U; L = 0 stays implicit
        assert d['A'].shape[0] == model.num_rows + model.packed_dim
        n_ineq = int((model.row_senses != '=').sum())
        assert d['mleq'] == n_ineq + model.packed_dim
        # converted rows sit after the model inequalities: row . X == X_ij
        rng = np.random.default_rng(5)
        X_half = rng.standard_normal((model.dim, model.dim))
        X = X_half + X_half.T
        block = d['A'][n_ineq:n_ineq + model.packed_dim]
        pi, pj = model._packed_ij()
        np.testing.assert_allclose(np.asarray(block @ X.ravel()).ravel(),
                                   X[pi, pj], atol=1e-12)
        np.testing.assert_array_equal(
            d['b'][n_ineq:n_ineq + model.packed_dim].ravel(), 1.)


# ---------------------------------------------------------------------------
# Integral-feasibility invariant
# ---------------------------------------------------------------------------
# Any valid relaxation must contain the moment matrix Y = y y^T, y = (1, chi_S),
# of every stable set S. Checking each row under its declared sense catches any
# sign or sense inversion in the lifting (see docs/plans/m-plus-kk-lifting.md).

# 5-cycle (nodes 1..5, edges of C5):
_C5_LP = """\
obj: x[1] + x[2] + x[3] + x[4] + x[5]
subject to
 edge1: x[1] + x[2] <= 1
 edge2: x[2] + x[3] <= 1
 edge3: x[3] + x[4] <= 1
 edge4: x[4] + x[5] <= 1
 edge5: x[1] + x[5] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
 0 <= x[3] <= 1
 0 <= x[4] <= 1
 0 <= x[5] <= 1
end
"""

_C5_EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
_TRIANGLE_EDGES = [(1, 2), (1, 3), (2, 3)]


@pytest.fixture
def c5_lp_path(tmp_path):
    path = tmp_path / "c5.lp"
    path.write_text(_C5_LP)
    return str(path)


def _stable_sets(n, edges):
    """All stable sets of the graph on nodes 1..n with the given edge list."""
    for r in range(n + 1):
        for S in itertools.combinations(range(1, n + 1), r):
            members = set(S)
            if not any(u in members and v in members for u, v in edges):
                yield members


def _canonical_row_values(model, Y):
    """<A_i, Y> for every row, under the canonical packed inner product."""
    pi, pj = model._packed_ij()
    weight = np.where(pi == pj, 1.0, 2.0)
    y_packed = weight * Y[pi, pj]
    return np.asarray(model.A @ y_packed).ravel()


def _assert_feasible(model, Y, label):
    vals = _canonical_row_values(model, Y)
    for i, (val, sense, rhs) in enumerate(
            zip(vals, model.row_senses, model.row_rhs)):
        if sense == '=':
            assert abs(val - rhs) <= 1e-9, \
                '%s: row %d: %g != %g' % (label, i, val, rhs)
        elif sense == '<':
            assert val <= rhs + 1e-9, \
                '%s: row %d: %g > %g' % (label, i, val, rhs)
        else:
            assert val >= rhs - 1e-9, \
                '%s: row %d: %g < %g' % (label, i, val, rhs)

    lb, ub = model._bounds_packed()
    pi, pj = model._packed_ij()
    entries = Y[pi, pj]
    assert (entries >= lb - 1e-9).all() and (entries <= ub + 1e-9).all(), \
        '%s: elementwise L <= Y <= U violated' % label


class TestIntegralFeasibility:

    @pytest.mark.parametrize('graph', ['triangle', 'c5'])
    @pytest.mark.parametrize('kwargs', [
        {'lift_mode': 'ls'},
        {'lift_mode': 'kk'},
        {'lift_mode': 'kk', 'include_squares': True},
    ], ids=['ls', 'kk', 'kk-squares'])
    def test_m_plus_contains_all_stable_sets(self, graph, kwargs,
                                             triangle_lp_path, c5_lp_path):
        if graph == 'triangle':
            path, n, edges = triangle_lp_path, 3, _TRIANGLE_EDGES
        else:
            path, n, edges = c5_lp_path, 5, _C5_EDGES

        model = _compute_m_plus_model(path, **kwargs)
        for S in _stable_sets(n, edges):
            y = np.zeros(model.dim)
            y[0] = 1.0
            for i in S:
                y[i] = 1.0
            _assert_feasible(model, np.outer(y, y), '%s S=%s' % (graph, S))

    @pytest.mark.parametrize('builder', [Theta_SDP, Theta_plus_SDP])
    def test_theta_models_contain_all_stable_sets(self, builder):
        G = nx.cycle_graph(5)   # 0-based nodes; SDP index = node + 1
        model = builder(G)
        edges = [(u + 1, v + 1) for u, v in G.edges()]
        for S in _stable_sets(len(G.nodes()), edges):
            y = np.zeros(model.dim)
            y[0] = 1.0
            for i in S:
                y[i] = 1.0
            _assert_feasible(model, np.outer(y, y),
                             '%s S=%s' % (builder.__name__, S))


# ---------------------------------------------------------------------------
# Regression: vectorized builder reproduces the original loop builder
# ---------------------------------------------------------------------------
# The fixtures in tests/fixtures/ were generated by the pre-vectorization
# loop builder (July 2026, see docs/plans/m-plus-kk-lifting.md).  The 'ls'
# mode must reproduce them exactly — same rows, same order, same values,
# same cut_classes — which is what keeps the MATLAB cut-class pipeline safe.

_FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')

_REGRESSION_CASES = [
    # (fixture stem, lp stem, builder kwargs)
    ('triangle_ls', 'triangle', {}),
    ('triangle_ls_skip', 'triangle', {'skip_func': lovasz_schrijver_filter}),
    ('triangle_ls_nobounds', 'triangle', {'lift_bounds': False}),
    ('clique_ls', 'clique', {}),
    ('nod_ls', 'nod', {}),
]


class TestRegressionAgainstLoopBuilder:

    @pytest.mark.parametrize('fixture,lp,kwargs', _REGRESSION_CASES,
                             ids=[c[0] for c in _REGRESSION_CASES])
    def test_ls_mode_reproduces_loop_builder(self, fixture, lp, kwargs):
        ref = np.load(os.path.join(_FIXTURE_DIR, fixture + '.npz'))
        model = _compute_m_plus_model(
            os.path.join(_FIXTURE_DIR, lp + '.lp'), **kwargs)

        assert model.dim == ref['dim']
        assert model.n == ref['n']
        np.testing.assert_array_equal(model.A.toarray(), ref['A'])
        np.testing.assert_array_equal(model.row_rhs, ref['row_rhs'])
        assert (model.row_senses == ref['row_senses']).all()
        np.testing.assert_array_equal(model.cut_classes, ref['cut_classes'])

    def test_generic_skip_func_matches_stock_filter(self):
        """A user-supplied callable takes the generic path but must select
        the same pairs as the special-cased stock filter."""
        path = os.path.join(_FIXTURE_DIR, 'triangle.lp')
        stock = _compute_m_plus_model(path, skip_func=lovasz_schrijver_filter)
        generic = _compute_m_plus_model(
            path, skip_func=lambda constr, i: lovasz_schrijver_filter(constr, i))
        np.testing.assert_array_equal(stock.A.toarray(), generic.A.toarray())
        np.testing.assert_array_equal(stock.row_rhs, generic.row_rhs)
        np.testing.assert_array_equal(stock.cut_classes, generic.cut_classes)


# ---------------------------------------------------------------------------
# LP parser sense guard
# ---------------------------------------------------------------------------

class TestParserSenseGuard:

    def test_geq_row_raises(self, tmp_path):
        path = tmp_path / "bad.lp"
        path.write_text(
            "obj: x[1] + x[2]\n"
            "subject to\n"
            " c1: x[1] + x[2] >= 1\n"
            "bounds\n"
            " 0 <= x[1] <= 1\n"
            " 0 <= x[2] <= 1\n"
            "end\n")
        with pytest.raises(ValueError, match='sense'):
            list(read_constr_from_lp(str(path)))


# ---------------------------------------------------------------------------
# M_+(K,K) mode
# ---------------------------------------------------------------------------

# 2-variable LP whose first K x K product row is hand-checkable:
#   c1 x c2 = (1 - x1 - x2)(1 - x1) >= 0, which with x_i^2 = x_i reduces to
#   X_11 + X_22 - X_12 <= 1.
_TWO_VAR_LP = """\
obj: x[1] + x[2]
subject to
 c1: x[1] + x[2] <= 1
 c2: x[1] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
end
"""


class TestKKMode:

    @pytest.fixture
    def two_var_lp_path(self, tmp_path):
        path = tmp_path / "two_var.lp"
        path.write_text(_TWO_VAR_LP)
        return str(path)

    def test_cut_classes_is_none(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path, lift_mode='kk')
        assert model.cut_classes is None

    def test_row_count_triangle(self, triangle_lp_path):
        # extended system: 3 LP rows + 2*3 bounds = 9 forms.
        # C(9,2) = 36 distinct pairs, minus the 3 trivial x_i(1 - x_i)
        # products, plus n + 1 = 4 homogenization equalities.
        model = _compute_m_plus_model(triangle_lp_path, lift_mode='kk')
        assert model.num_rows == 4 + 36 - 3
        model_sq = _compute_m_plus_model(triangle_lp_path, lift_mode='kk',
                                         include_squares=True)
        assert model_sq.num_rows == 4 + 36 - 3 + 9

    def test_no_zero_rows(self, triangle_lp_path):
        model = _compute_m_plus_model(triangle_lp_path, lift_mode='kk')
        nnz_per_row = np.diff(model.A.tocsr().indptr)
        assert (nnz_per_row > 0).all()

    def test_lp_pair_product_hand_check(self, two_var_lp_path):
        """First 'kk' inequality row is the c1 x c2 product, expanded by hand."""
        model = _compute_m_plus_model(two_var_lp_path, lift_mode='kk')
        # rows 0..2 are the equalities; row 3 is the pair (c1, c2)
        row = model.A[3].toarray().ravel()
        expected = np.zeros(model.packed_dim)
        expected[mat_idx(1, 1)] = 1.0    # X_11
        expected[mat_idx(2, 2)] = 1.0    # X_22
        expected[mat_idx(2, 1)] = -0.5   # -X_12 (canonical implicit x2)
        np.testing.assert_array_equal(row, expected)
        assert model.row_rhs[3] == 1.0
        assert model.row_senses[3] == '<'

    def test_x_nonneg_products_present(self, two_var_lp_path):
        """The e_i x e_j pairs must appear as explicit X_ij >= 0 rows
        (canonically -X_ij <= 0): MOSEK bar-variables have no elementwise
        bounds, so these rows are the only way to impose X >= 0 there."""
        model = _compute_m_plus_model(two_var_lp_path, lift_mode='kk')
        expected = np.zeros(model.packed_dim)
        expected[mat_idx(2, 1)] = -0.5
        A = model.A.toarray()
        hits = [r for r in range(model.num_rows)
                if np.array_equal(A[r], expected)
                and model.row_rhs[r] == 0.0 and model.row_senses[r] == '<']
        assert len(hits) == 1

    def test_pair_filter(self, triangle_lp_path):
        """Keeping only LP x bound pairs reproduces the row count of the
        'ls' model without bound lifting (classes 1-2 products)."""
        model = _compute_m_plus_model(
            triangle_lp_path, lift_mode='kk',
            pair_filter=lambda H, K, L: (K < 3) & (L >= 3))
        ls = _compute_m_plus_model(triangle_lp_path, lift_bounds=False)
        assert model.num_rows == ls.num_rows

    def test_skip_func_rejected_in_kk(self, triangle_lp_path):
        with pytest.raises(ValueError, match='skip_func'):
            _compute_m_plus_model(triangle_lp_path, lift_mode='kk',
                                  skip_func=lovasz_schrijver_filter)

    def test_pair_filter_rejected_in_ls(self, triangle_lp_path):
        with pytest.raises(ValueError, match='pair_filter'):
            _compute_m_plus_model(triangle_lp_path,
                                  pair_filter=lambda H, K, L: K >= 0)

    def test_invalid_mode_rejected(self, triangle_lp_path):
        with pytest.raises(ValueError, match='lift_mode'):
            _compute_m_plus_model(triangle_lp_path, lift_mode='bogus')


# ---------------------------------------------------------------------------
# to_lp export (M-operator LP)
# ---------------------------------------------------------------------------

class TestToLP:

    @pytest.fixture
    def triangle_model(self, triangle_lp_path):
        return _compute_m_plus_model(triangle_lp_path)

    def test_shapes_and_keys(self, triangle_model):
        d = triangle_model.to_lp()
        assert d['num_cols'] == triangle_model.packed_dim
        assert d['num_rows'] == triangle_model.num_rows
        assert d['A'].shape == (d['num_rows'], d['num_cols'])
        assert d['obj'].shape == (d['num_cols'],)
        assert len(d['col_names']) == d['num_cols']
        assert d['obj_sense'] == triangle_model.obj_sense
        assert d['obj_offset'] == triangle_model.obj_offset

    def test_off_diagonal_coefficients_doubled(self, triangle_model):
        """A row's LP coefficients are the canonical packed entries with the
        implicit x2 made explicit on off-diagonal positions."""
        model = triangle_model
        d = model.to_lp()
        s = model._lp_scale()
        pi, pj = model._packed_ij()
        assert ((s == 2.0) == (pi != pj)).all()
        np.testing.assert_array_equal(
            d['A'].toarray(), model.A.toarray() * s)
        # objective: C is diagonal here, so no doubling applies
        C_diag = np.asarray(model._C_full().todense()).diagonal()
        diag_idx = np.array([mat_idx(i, i) for i in range(model.dim)])
        np.testing.assert_array_equal(d['obj'][diag_idx], C_diag)

    def test_lp_inner_product_matches_canonical(self, triangle_model):
        """obj . svec-free packed y == <C, Y> for a random symmetric Y."""
        model = triangle_model
        d = model.to_lp()
        rng = np.random.default_rng(3)
        Y_half = rng.standard_normal((model.dim, model.dim))
        Y = Y_half + Y_half.T
        pi, pj = model._packed_ij()
        y_packed = Y[pi, pj]   # plain packed entries, no scaling
        C_full = np.asarray(model._C_full().todense())
        C_sym = C_full + C_full.T - np.diag(C_full.diagonal())
        np.testing.assert_allclose(d['obj'] @ y_packed,
                                   np.sum(C_sym * Y), atol=1e-10)

    def test_bounds(self, triangle_model):
        # M+ models carry L = 0, mapped onto the variable lower bounds
        d = triangle_model.to_lp()
        assert (d['lb'] == 0.0).all()
        assert (d['ub'] == np.inf).all()
        triangle_model.L = None
        d_free = triangle_model.to_lp()
        assert (d_free['lb'] == -np.inf).all()

    def test_gurobi_m_operator_bound(self, triangle_lp_path):
        """The M(K,K)-operator LP bound on the triangle FRAC must be 1.0."""
        gp = pytest.importorskip('gurobipy')
        model = _compute_m_plus_model(triangle_lp_path, lift_mode='kk',
                                      include_squares=True)
        d = model.to_lp()
        lp = gp.Model('m_operator')
        lp.Params.OutputFlag = 0
        y = lp.addMVar(d['num_cols'], lb=d['lb'], ub=d['ub'])
        lp.setObjective(d['obj'] @ y + d['obj_offset'],
                        gp.GRB.MINIMIZE if d['obj_sense'] > 0
                        else gp.GRB.MAXIMIZE)
        for sense in ('<', '>', '='):
            mask = d['senses'] == sense
            if mask.any():
                lp.addMConstr(d['A'][mask], y, sense, d['rhs'][mask])
        lp.optimize()
        assert lp.Status == gp.GRB.OPTIMAL
        assert abs(-lp.ObjVal - 1.0) < 1e-6
