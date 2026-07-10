# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# SDP model generation via the Lovász-Schrijver M+ lift-and-project operator
# and Theta/Theta+ SDP formulations.

import os, time, itertools
import numpy as np
import re
from scipy import sparse
from scipy.sparse import lil_matrix

from pyModules.SDPModel import SDPModel
from tqdm import tqdm

"""
    The main functions in this file build SDPModel objects representing the SDP
    programs where matrices are stored in lower-triangular format.

        m_plus_lifting(): Lovász-Schrijver Lift-and-Project formulation.

             Theta_SDP(): Theta function SDP formulation.

        Theta_plus_SDP(): Theta+ function SDP formulation.

    To serialize to a specific solver format use the SDPModel export methods:

        model.to_sdpnal()          -- returns dict suitable for SDPNAL+
        model.to_adal()            -- returns dict suitable for the internal ADMM solver
        model.write_sdpnal_mat()   -- saves a SDPNAL+ .mat file

    Indexing conventions
    --------------------
    In the (n+1) x (n+1) SDP matrix variable Y, index 0 is the homogenization
    row/column and graph node i (0-based, see Graphs.py) lives at index i+1,
    i.e. x_i = Y[i+1, i+1].

        m_plus_lifting() parses an LP file whose variables x[1]..x[n] are
        already 1-based, so it passes them to mat_idx(i, var) with no shift —
        the graph itself never enters the M+ builder.

        Theta_SDP() / Theta_plus_SDP() build from a 0-based networkx graph,
        hence the mat_idx(i + 1, j + 1) shifts in their constraint loops.
"""


# ---------------------------------------------------------------------------
# LP file parsing helpers
# ---------------------------------------------------------------------------

def read_index_from_lp(path_to_file):
    with open(path_to_file) as f:
        exp_id = re.compile(r"^\s*.*:")

        i = 0
        for line in f:
            i += 1
            if '\\' in line or not line:
                continue
            if 'bounds' in line.lower():
                yield i
                break

            if re.search(exp_id, line):
                yield i


def read_expr_from_lp(path_to_file):
    exp_id = re.compile(r"^\s*[a-zA-Z]+[0-9]*:\s*")
    exp_idx = read_index_from_lp(path_to_file)
    try:
        with open(path_to_file) as f:
            i = 0
            line = ''
            start = None
            while True:
                if not start:
                    start = next(exp_idx)
                ends = next(exp_idx)
                s = ''
                while i < start:
                    i += 1
                    line = f.readline().strip().lower()
                while start <= i < ends:
                    match = re.search(exp_id, line)
                    if match:
                        line = line.replace(match.group(), '')
                    if 'subject to' not in line:
                        s += line
                    i += 1
                    line = f.readline().strip().lower()

                yield s
                start = ends

    except StopIteration:
        pass


def read_constr_from_lp(path_to_file):
    """Parse an LP file and yield (coefficient_dict, rhs) pairs.

    The first yielded item is the objective (a coefficient_dict with no rhs).
    Subsequent items are constraints: ({var: coeff}, rhs_value).

    Only '<=' rows are supported (the restricted dialect written by
    LinearFormulations); the lifting assumes every row is an upper bound,
    so any other sense raises ValueError instead of being silently flipped.
    """
    constr = re.compile(r"[\+\-]?\s?[\d]*\s?x\[?\d+\]?")
    rhs = re.compile(r"([\<\>\=]+)\s*([\d]+)")

    lines = read_expr_from_lp(path_to_file)
    for l in lines:
        dic = {}
        rhs_val = re.search(rhs, l)
        if rhs_val:
            sense, b_str = rhs_val.groups()
            if '<' not in sense:
                raise ValueError(
                    "Unsupported constraint sense %r in %s: only '<=' rows "
                    "are handled" % (sense, path_to_file))
            b = float(b_str)

        match = re.findall(constr, l)
        for m in match:
            c = m.replace(' ', '').split('x')
            try:
                val = float(c[0])
            except ValueError:
                if not c[0] or '+' == c[0]:
                    val = 1.
                else:
                    val = -1.
            var = int(c[1].replace('[', '').replace(']', ''))
            dic[var] = val

        yield (dic, b) if rhs_val else dic


# ---------------------------------------------------------------------------
# Sparse matrix helpers
# ---------------------------------------------------------------------------

def sp_unique(sp_matrix, axis=0):
    """Return a sparse matrix with the unique rows (axis=0) or columns (axis=1)."""
    if axis == 1:
        sp_matrix = sp_matrix.T

    old_format = sp_matrix.getformat()
    dt = np.dtype(sp_matrix)
    ncols = sp_matrix.shape[1]

    if old_format != 'lil':
        sp_matrix = sp_matrix.tolil()

    _, ind = np.unique(sp_matrix.data + sp_matrix.rows, return_index=True)
    rows = sp_matrix.rows[ind]
    data = sp_matrix.data[ind]
    nrows_uniq = data.shape[0]

    sp_matrix = sparse.lil_matrix((nrows_uniq, ncols), dtype=dt)
    sp_matrix.data = data
    sp_matrix.rows = rows

    ret = sp_matrix.asformat(old_format)
    if axis == 1:
        ret = ret.T
    return ret, ind


def mat_idx(i, j):
    """Map 2D symmetric matrix indices to 1D packed (lower-triangular) index."""
    if i >= j:
        return ((i + 1) * i) // 2 + j
    else:
        return ((j + 1) * j) // 2 + i


def push(_i, _j, _d, _Bt, step, dim):
    coo = sparse.coo_matrix((_d, (_i, _j)), shape=(dim, step))
    del _i, _j, _d
    _Bt.append(coo)
    return [], [], [], 0


# ---------------------------------------------------------------------------
# Constraint class mapping helpers
# ---------------------------------------------------------------------------

def lovasz_schrijver_filter(constr, i):
    """Return True if variable i appears in constr (used to skip edge lifting)."""
    return i in constr


def map_edge_constraints(filename):
    constrs = list(read_constr_from_lp(filename))
    n = len(constrs[0])
    map = {}
    idx = 1
    for constr, b in constrs[1:]:
        for j in range(1, n + 1):
            if lovasz_schrijver_filter(constr, j):
                continue
            map[idx] = 'lift_edge1:1'
            idx += 1
            map[idx] = 'lift_edge1:2'
            idx += 1

    for _ in itertools.combinations(range(1, n + 1), 2):
        map[idx] = 'bound2'
        idx += 1
        map[idx] = 'bound4'
        idx += 1
    return map


def map_clique_constraints(filename):
    constrs = list(read_constr_from_lp(filename))
    n = len(constrs[0])
    map = {}
    idx = 1
    for constr, b in constrs[1:]:
        for j in range(1, n + 1):
            if j in constr:
                map[idx] = 'lift_clq1:1'
            else:
                map[idx] = 'lift_clq2:1'
            idx += 1

            if j in constr:
                map[idx] = 'lift_clq1:2'
            else:
                map[idx] = 'lift_clq2:2'
            idx += 1

    for _ in itertools.combinations(range(1, n + 1), 2):
        map[idx] = 'bound2'
        idx += 1
        map[idx] = 'bound4'
        idx += 1
    return map


def map_nod_constraints(filename):
    constrs = list(read_constr_from_lp(filename))
    n = len(constrs[0])
    map = {}
    idx = 1
    for constr, b in constrs[1:]:
        nod = 0
        for i in constr:
            if constr[i] != 1.:
                nod = i
        for j in range(1, n + 1):
            if j == nod:
                map[idx] = 'lift_nod1:1'
            elif j in constr:
                map[idx] = 'lift_nod2:1'
            else:
                map[idx] = 'lift_nod3:1'
            idx += 1

            if j == nod:
                map[idx] = 'lift_nod1:2'
            elif j in constr:
                map[idx] = 'lift_nod2:2'
            else:
                map[idx] = 'lift_nod3:2'
            idx += 1

    for _ in itertools.combinations(range(1, n + 1), 2):
        map[idx] = 'bound2'
        idx += 1
        map[idx] = 'bound4'
        idx += 1
    return map


# ---------------------------------------------------------------------------
# M+ Lift-and-Project
# ---------------------------------------------------------------------------
# Every constraint of K is kept in nonnegative homogenized form
# g = (beta, -alpha) in R^{1+n}, meaning beta - alpha.x >= 0 on y = (1, x).
# The product of two such forms is valid on the moment matrix
# Y = [1 x^T; x X]:
#
#     g_k^T Y g_l >= 0,   i.e.   <sym(g_k g_l^T), Y> >= 0,
#
# emitted in the model's '<' sense as <-sym(g_k g_l^T), Y> <= 0, then
# mechanically rewritten under the homogenization equalities
# (use of X_{0i} = X_{ii} folds packed column (i,0) into (i,i); X_{00} = 1
# moves the (0,0) coefficient to the rhs).  Which pairs (k, l) are formed
# is controlled by lift_mode: 'ls' reproduces the classic Lovász-Schrijver
# operator (LP row x bound, plus bound x bound — cut classes 1-4), 'kk' is
# the full M_+(K,K) with all distinct pairs of the extended system.


def _homogenized_system(lp_rows, n):
    """Build the homogenized constraint system of K as a CSR matrix H.

    Row ordering (pair filters may rely on it):

        rows 0..m-1      : LP rows  a.x <= b      -> (b, -a)
        rows m..m+n-1    : bounds   x_i >= 0      -> e_i        (i = 1..n)
        rows m+n..m+2n-1 : bounds   1 - x_i >= 0  -> (1, -e_i)  (i = 1..n)
    """
    m = len(lp_rows)
    ri = []; ci = []; d = []
    for r, (constr, b) in enumerate(lp_rows):
        ri.append(r); ci.append(0); d.append(b)
        for var, coeff in constr.items():
            ri.append(r); ci.append(var); d.append(-coeff)
    for i in range(1, n + 1):
        ri.append(m + i - 1); ci.append(i); d.append(1.)
    for i in range(1, n + 1):
        ri.append(m + n + i - 1); ci.append(0); d.append(1.)
        ri.append(m + n + i - 1); ci.append(i); d.append(-1.)
    return sparse.coo_matrix(
        (d, (ri, ci)), shape=(m + 2 * n, n + 1)).tocsr()


def _ls_pairs(H, m, n, lift_bounds, skip_func, lp_dicts):
    """Pair index arrays reproducing the classic LS operator row order.

    Per LP row r, per variable i ascending: (r, x_i) then (r, 1-x_i)
    — classes 1 and 2; then per i < j: (1-x_i, x_j) and (1-x_i, 1-x_j)
    — classes 3 and 4.  Returns (K_idx, L_idx, cut_classes).

    Cut classes, for an LP row a.x <= b (after the homogenization
    identities X_0i = X_ii, X_00 = 1):

      1  (b - a.x) x_i >= 0        ->  sum_j a_j X_ij <= b X_ii
                                       (the row conditional on x_i = 1)
      2  (b - a.x)(1 - x_i) >= 0   ->  sum_j a_j (X_jj - X_ij) <= b (1 - X_ii)
                                       (the row conditional on x_i = 0)
      3  (1 - x_i) x_j >= 0        ->  X_ij <= X_jj
                                       (i < j only: one orientation per pair,
                                       against the larger-index diagonal;
                                       NOT X_ij >= 0 — the x_i x_j products
                                       are never formed in 'ls' mode)
      4  (1 - x_i)(1 - x_j) >= 0   ->  X_ii + X_jj - X_ij <= 1  (McCormick)

    Used by MATLAB cut selection (kelley_cutting_plane.m).
    """
    K_parts = []; L_parts = []; cls_parts = []
    all_vars = np.arange(1, n + 1)

    for r in range(m):
        if skip_func is lovasz_schrijver_filter:
            # vectorized equivalent: skip i iff column i is in row r's support
            support = H.indices[H.indptr[r]:H.indptr[r + 1]]
            vi = np.setdiff1d(all_vars, support)
        elif skip_func:
            vi = np.array([i for i in all_vars
                           if not skip_func(lp_dicts[r], i)], dtype=np.int64)
        else:
            vi = all_vars
        if len(vi) == 0:
            continue
        L = np.empty(2 * len(vi), dtype=np.int64)
        L[0::2] = m + vi - 1          # x_i >= 0 multiplier      (class 1)
        L[1::2] = m + n + vi - 1      # 1 - x_i >= 0 multiplier  (class 2)
        K_parts.append(np.full(2 * len(vi), r, dtype=np.int64))
        L_parts.append(L)
        cls_parts.append(np.tile([1., 2.], len(vi)))

    if lift_bounds:
        oi, oj = np.triu_indices(n, k=1)   # 0-based: pairs i < j with i=oi+1, j=oj+1
        npairs = len(oi)
        K = np.repeat(m + n + oi, 2)       # (1 - x_i) multiplier
        L = np.empty(2 * npairs, dtype=np.int64)
        L[0::2] = m + oj                   # x_j        -> class 3: X_ij <= X_jj
        L[1::2] = m + n + oj               # (1 - x_j)  -> class 4: McCormick
        K_parts.append(K)
        L_parts.append(L)
        cls_parts.append(np.tile([3., 4.], npairs))

    if not K_parts:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty.copy(), np.empty(0)
    return (np.concatenate(K_parts), np.concatenate(L_parts),
            np.concatenate(cls_parts))


def _kk_pairs(num_forms, m, lift_bounds, include_squares):
    """All unordered pairs over the extended system for M_+(K,K)."""
    K, L = np.triu_indices(num_forms, k=0 if include_squares else 1)
    if not lift_bounds:
        keep = ~((K >= m) & (L >= m))
        K, L = K[keep], L[keep]
    return K, L


def _pair_products(H, K_idx, L_idx, packed_dim, step=100000):
    """Lifted product rows for the given pairs, in '<' sense.

    Row p is <-sym(g_{K_idx[p]} g_{L_idx[p]}^T), Y> <= u_p after the
    diag-identity rewrites.  Computed as a row-wise Kronecker product on
    the CSR internals of H, in chunks of `step` pairs to bound peak memory.
    Returns (B, u) with B a (len(K_idx) x packed_dim) CSR matrix.
    """
    blocks = []; u_parts = []
    chunks = range(0, len(K_idx), step)
    if len(chunks) > 1:
        chunks = tqdm(chunks, desc='M+ Lifting', unit='chunk')
    for s in chunks:
        B, u = _pair_products_chunk(
            H, K_idx[s:s + step], L_idx[s:s + step], packed_dim)
        blocks.append(B)
        u_parts.append(u)
    if not blocks:
        return sparse.csr_matrix((0, packed_dim)), np.empty(0)
    return sparse.vstack(blocks, format='csr'), np.concatenate(u_parts)


def _pair_products_chunk(H, K_idx, L_idx, packed_dim):
    indptr, indices, data = H.indptr, H.indices, H.data
    kstart = indptr[K_idx]
    lstart = indptr[L_idx]
    kcnt = indptr[K_idx + 1] - kstart
    lcnt = indptr[L_idx + 1] - lstart

    # Cross product of the two nonzero supports for every pair, ragged-style:
    # triplet t of pair p combines k-entry (t // lcnt[p]) with l-entry (t % lcnt[p])
    cnt = kcnt * lcnt
    num_pairs = len(K_idx)
    pair_of = np.repeat(np.arange(num_pairs), cnt)
    offsets = np.concatenate(([0], np.cumsum(cnt)))
    t = np.arange(offsets[-1]) - offsets[pair_of]
    nl = lcnt[pair_of]
    ki = kstart[pair_of] + t // nl
    li = lstart[pair_of] + t % nl

    p = indices[ki].astype(np.int64)   # row index in Y
    q = indices[li].astype(np.int64)   # col index in Y
    d = -(data[ki] * data[li])         # negate: emit in '<' sense
    d = np.where(p == q, d, 0.5 * d)   # canonical implicit x2 off-diagonal

    mx = np.maximum(p, q)
    mn = np.minimum(p, q)

    # X_00 = 1: the (0,0) coefficient moves to the rhs
    u = np.zeros(num_pairs)
    on_00 = mx == 0
    np.subtract.at(u, pair_of[on_00], d[on_00])

    # X_0i = X_ii: packed column (i,0) folds into (i,i); the off-diagonal
    # canonical weight 0.5 becomes 1 on the diagonal, hence the x2
    on_0col = (mn == 0) & (mx > 0)
    packed = (mx * (mx + 1)) // 2 + mn
    packed[on_0col] += mx[on_0col]     # (i,0) -> (i,i): mat_idx shifts by i
    d[on_0col] *= 2.0

    keep = ~on_00
    B = sparse.coo_matrix(
        (d[keep], (pair_of[keep], packed[keep])),
        shape=(num_pairs, packed_dim)).tocsr()
    return B, u


def _homogenization_equalities(n, packed_dim):
    """Equality rows X_00 = 1 and X_0i = X_ii (i = 1..n), canonical form."""
    i = np.arange(1, n + 1)
    ri = np.concatenate(([0], i, i))
    ci = np.concatenate(([0],
                         i * (i + 1) // 2,       # mat_idx(i, 0)
                         i * (i + 1) // 2 + i))  # mat_idx(i, i)
    d = np.concatenate(([1.], np.full(n, 0.5), np.full(n, -1.)))
    A_eq = sparse.coo_matrix(
        (d, (ri, ci)), shape=(n + 1, packed_dim)).tocsr()
    b = np.concatenate(([1.], np.zeros(n)))
    return A_eq, b


def _compute_m_plus_model(path_to_file, step=100000, lift_bounds=True,
                          skip_func=None, lift_mode='ls',
                          include_squares=False, pair_filter=None):
    """Build the M+ SDP model as an SDPModel.

    Parameters
    ----------
    step            : pair-chunk size bounding peak memory during construction
    lift_bounds     : include the bound x bound products
    skip_func       : ('ls' only) function(constr_dict, i) -> True to skip
                      lifting LP row `constr` by both x_i and 1 - x_i
    lift_mode       : 'ls' — classic Lovász-Schrijver pairs, cut classes 1-4
                      (LP row x bound and bound x bound); 'kk' — full
                      M_+(K,K), all distinct pairs of {LP rows u bounds},
                      O((m+2n)^2) rows, cut_classes=None
    include_squares : ('kk' only) also emit the k = l products; redundant
                      under Y >= 0 (PSD) but part of M(K,K) proper when the
                      PSD constraint is dropped (see SDPModel.to_lp)
    pair_filter     : ('kk' only) function(H, K_idx, L_idx) -> boolean keep
                      mask over candidate pairs, applied before any numeric
                      work; row ordering of H is documented in
                      _homogenized_system

    All constraint matrices are stored in lower-triangular form with no
    solver-specific scaling; call model.to_sdpnal() or model.to_adal() to
    obtain a solver-ready representation.
    """
    constrs = list(read_constr_from_lp(path_to_file))

    obj = constrs[0]
    lp_rows = constrs[1:]
    C_dense = -np.diag([0.] + list(obj.values()))
    n = C_dense.shape[0] - 1
    dim = n + 1
    packed_dim = (dim * (dim + 1)) // 2

    # C is diagonal — lower-tri = diagonal itself
    C = sparse.csr_matrix(C_dense)

    m = len(lp_rows)
    H = _homogenized_system(lp_rows, n)

    if lift_mode == 'ls':
        if pair_filter is not None:
            raise ValueError("pair_filter is only supported in 'kk' mode; "
                             "use skip_func in 'ls' mode")
        K_idx, L_idx, cut_classes = _ls_pairs(
            H, m, n, lift_bounds, skip_func, [c for c, _ in lp_rows])
    elif lift_mode == 'kk':
        if skip_func is not None:
            raise ValueError("skip_func is only supported in 'ls' mode; "
                             "use pair_filter in 'kk' mode")
        K_idx, L_idx = _kk_pairs(H.shape[0], m, lift_bounds, include_squares)
        if pair_filter is not None:
            keep = np.asarray(pair_filter(H, K_idx, L_idx), dtype=bool)
            K_idx, L_idx = K_idx[keep], L_idx[keep]
        cut_classes = None
    else:
        raise ValueError("lift_mode must be 'ls' or 'kk', got %r" % lift_mode)

    B, u = _pair_products(H, K_idx, L_idx, packed_dim, step=step)

    if lift_mode == 'kk':
        # drop products that became identically zero under the diag-identity
        # rewrites (e.g. x_i (1 - x_i) >= 0 reduces to 0 >= 0)
        B.eliminate_zeros()
        keep = np.diff(B.indptr) > 0
        B, u = B[keep], u[keep]

    A_eq, b = _homogenization_equalities(n, packed_dim)
    A = sparse.vstack([A_eq, B], format='csr')

    m_eq, m_ineq = A_eq.shape[0], B.shape[0]

    # L = 0: X >= 0 comes from the x_i x_j multiplications of the operator
    return SDPModel(
        dim=dim, n=n, num_rows=m_eq + m_ineq,
        C=C, A=A,
        row_rhs=np.concatenate([b, u]),
        row_senses=np.array(['='] * m_eq + ['<'] * m_ineq, dtype='U1'),
        L=0.,
        cut_classes=cut_classes,
    )


def m_plus_lifting(filename, model_out_dir='', work_dir='', step=100000,
                   lift_bounds=True, skip_func=None, do_split=False,
                   lift_mode='ls', include_squares=False, pair_filter=None):
    """Apply the Lovász-Schrijver M+ operator to an LP and save the SDP model.

    Parameters
    ----------
    filename       : LP file name (without directory)
    model_out_dir  : directory where the .mat output is written
    work_dir       : directory containing the LP file
    step           : pair-chunk size for building the sparse constraint matrix
    lift_bounds    : whether to include bound x bound lifting constraints
                     (McCormik inequalities)
    skip_func      : optional function(constr, i) -> bool to skip constraints
                     ('ls' mode only)
    do_split       : split Bt into 4 files when True (for large models)
    lift_mode      : 'ls' (classic operator, default) or 'kk' (full M_+(K,K));
                     see _compute_m_plus_model
    include_squares: 'kk' mode only — also emit the squared-constraint products
    pair_filter    : 'kk' mode only — boolean mask function over candidate pairs

    Returns
    -------
    model : SDPModel in canonical lower-triangular format
    """
    start = time.time()
    print('File: ', filename)
    path_to_file = os.path.join(work_dir, filename)
    out_filename = '.'.join(filename.split('.')[:-1])

    model = _compute_m_plus_model(path_to_file, step=step,
                                  lift_bounds=lift_bounds, skip_func=skip_func,
                                  lift_mode=lift_mode,
                                  include_squares=include_squares,
                                  pair_filter=pair_filter)

    out_path = os.path.join(model_out_dir, out_filename)
    print('Saving MAT at: %s.mat' % out_path)
    model.write_sdpnal_mat(out_path, do_split=do_split)

    elapsed = time.time() - start
    eq_rows   = int((model.row_senses == '=').sum())
    ineq_rows = int((model.row_senses != '=').sum())
    print('Finished! Total time: %.2f' % elapsed)
    print('Order of the Matrix Variable: %d' % model.dim)
    print('Number of Equalities: %d' % eq_rows)
    print('Number of Inequalities: %d' % ineq_rows)
    print('-' * 12)

    return model


# ---------------------------------------------------------------------------
# Theta SDP
# ---------------------------------------------------------------------------

def Theta_SDP(G, step=10000):
    """Build and return the Theta SDP for graph G as a canonical SDPModel.

    Constraints (all equalities):
      - X_{i+1,j+1} = 0         for each edge (i, j) in G
      - X_{0,0} = 1              (trace normalisation)
      - X_{0,i+1} = X_{i+1,i+1} for each node i
    Objective: min -sum_i X_{i+1,i+1}

    Use model.to_adal() to obtain the ADMM-compatible format, or
    model.write_sdpnal_mat(path) to save a SDPNAL+ .mat file.
    """
    n = len(G.nodes())
    dim = n + 1
    packed_dim = (dim * (dim + 1)) // 2

    C_lil = lil_matrix((dim, dim))
    for i in range(n):
        C_lil[i + 1, i + 1] = -1.
    C = C_lil.tocsr()

    _At = []
    _b  = []
    _i = []; _j = []; _d = []
    p = 0

    for i, j in G.edges():
        # graph nodes are 0-based; +1 skips the homogenization index 0
        ii, jj = max(i, j), min(i, j)
        _i.append(mat_idx(ii + 1, jj + 1)); _j.append(p); _d.append(0.5)
        _b.append(0.)
        p += 1
        if p == step:
            _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    _i.append(mat_idx(0, 0)); _j.append(p); _d.append(1.)
    _b.append(1.)
    p += 1
    if p == step:
        _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    for i in G.nodes():
        _i += [mat_idx(i + 1, 0), mat_idx(i + 1, i + 1)]
        _j += [p, p]
        _d += [0.5, -1.]
        _b.append(0.)
        p += 1
        if p == step:
            _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    if _i:
        _i, _j, _d, p = push(_i, _j, _d, _At, p % step, packed_dim)

    num_rows = len(_b)
    A = sparse.hstack(_At, format='csc').T.tocsr()
    row_rhs    = np.array(_b)
    row_senses = np.full(num_rows, '=', dtype='U1')

    return SDPModel(dim=dim, n=n, num_rows=num_rows,
                    C=C, A=A, row_rhs=row_rhs, row_senses=row_senses)


# ---------------------------------------------------------------------------
# Theta+ SDP
# ---------------------------------------------------------------------------

def Theta_plus_SDP(G, step=10000):
    """Build and return the Theta+ SDP for graph G as a canonical SDPModel.

    Extends Theta with non-edge nonnegativity inequalities:
      - <A_ij, X> >= 0  (i.e. X_{i+1,j+1} >= 0)  for each non-edge (i, j)
    followed by the same equality constraints as Theta_SDP.

    Use model.write_sdpnal_mat(path) to save a SDPNAL+ .mat file, or
    model.to_adal() for the ADMM-compatible format.
    """
    n = len(G.nodes())
    dim = n + 1
    packed_dim = (dim * (dim + 1)) // 2

    edge_set = set(G.edges()) | {(j, i) for i, j in G.edges()}
    complement_edges = [
        (i, j) for i, j in itertools.combinations(G.nodes(), 2)
        if (i, j) not in edge_set
    ]

    C_lil = lil_matrix((dim, dim))
    for i in range(n):
        C_lil[i + 1, i + 1] = -1.
    C = C_lil.tocsr()

    # --- Inequality constraints: non-edge non-positivity ---
    _Bt = []
    _u  = []
    _ib = []; _jb = []; _db = []
    l = 0

    for i, j in complement_edges:
        # graph nodes are 0-based; +1 skips the homogenization index 0
        ii, jj = max(i, j), min(i, j)
        _ib.append(mat_idx(ii + 1, jj + 1)); _jb.append(l); _db.append(0.5)
        _u.append(0.)
        l += 1
        if l % step == 0:
            _ib, _jb, _db, l = push(_ib, _jb, _db, _Bt, step, packed_dim)

    if _ib:
        _ib, _jb, _db, l = push(_ib, _jb, _db, _Bt, l % step, packed_dim)

    # --- Equality constraints: same as Theta_SDP ---
    _At = []
    _b  = []
    _i = []; _j = []; _d = []
    p = 0

    for i, j in G.edges():
        ii, jj = max(i, j), min(i, j)
        _i.append(mat_idx(ii + 1, jj + 1)); _j.append(p); _d.append(0.5)
        _b.append(0.)
        p += 1
        if p == step:
            _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    _i.append(mat_idx(0, 0)); _j.append(p); _d.append(1.)
    _b.append(1.)
    p += 1
    if p == step:
        _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    for i in G.nodes():
        _i += [mat_idx(i + 1, 0), mat_idx(i + 1, i + 1)]
        _j += [p, p]
        _d += [0.5, -1.]
        _b.append(0.)
        p += 1
        if p == step:
            _i, _j, _d, p = push(_i, _j, _d, _At, step, packed_dim)

    if _i:
        _i, _j, _d, p = push(_i, _j, _d, _At, p % step, packed_dim)

    # Combine: inequalities first, then equalities
    A_ineq = sparse.hstack(_Bt, format='csc').T.tocsr() if _Bt else sparse.csr_matrix((0, packed_dim))
    A_eq   = sparse.hstack(_At, format='csc').T.tocsr()
    A      = sparse.vstack([A_ineq, A_eq], format='csr')

    m_ineq   = len(_u)
    m_eq     = len(_b)
    num_rows = m_ineq + m_eq

    row_rhs    = np.concatenate([np.array(_u), np.array(_b)])
    row_senses = np.array(['>'] * m_ineq + ['='] * m_eq, dtype='U1')

    return SDPModel(dim=dim, n=n, num_rows=num_rows,
                    C=C, A=A, row_rhs=row_rhs, row_senses=row_senses,
                    L=0.)


if __name__ == "__main__":
    pass
