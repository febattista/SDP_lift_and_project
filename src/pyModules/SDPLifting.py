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
    programs in a solver-agnostic canonical lower-triangular format.

        m_plus_lifting(): Lovász-Schrijver Lift-and-Project formulation.

             Theta_SDP(): Theta function SDP formulation.

        Theta_plus_SDP(): Theta+ function SDP formulation.

    To serialize to a specific solver format use the SDPModel export methods:

        model.to_sdpnal()          -- returns dict suitable for SDPNAL+
        model.to_adal()            -- returns dict suitable for the internal ADMM solver
        model.write_sdpnal_mat()   -- saves a SDPNAL+ .mat file
"""


# ---------------------------------------------------------------------------
# LP file parsing helpers
# ---------------------------------------------------------------------------

def read_index_from_lp(path_to_file):
    try:
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

    except IOError:
        print('Unable to open %s. Try again.' % path_to_file)


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

    except IOError:
        print('Unable to open %s. Try again.' % path_to_file)
    except StopIteration:
        pass


def read_constr_from_lp(path_to_file):
    """Parse an LP file and yield (coefficient_dict, rhs) pairs.

    The first yielded item is the objective (a coefficient_dict with no rhs).
    Subsequent items are constraints: ({var: coeff}, rhs_value).
    """
    constr = re.compile(r"[\+\-]?\s?[\d]*\s?x\[?\d+\]?")
    rhs = re.compile(r"[\<\>\=]+\s*([\d]+)")

    lines = read_expr_from_lp(path_to_file)
    for l in lines:
        dic = {}
        rhs_val = re.search(rhs, l)
        if rhs_val:
            b = float(rhs_val.groups()[0])

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

def _compute_m_plus_model(path_to_file, step=100000, lift_bounds=True, skip_func=None):
    """Build the M+ SDP model as a canonical SDPModel.

    All constraint matrices are stored in lower-triangular form with no
    solver-specific scaling; call model.to_sdpnal() or model.to_adal() to
    obtain a solver-ready representation.
    """
    constrs = list(read_constr_from_lp(path_to_file))

    obj = constrs[0]
    C_dense = -np.diag([0.] + list(obj.values()))
    n = C_dense.shape[0] - 1
    dim = n + 1
    packed_dim = (dim * (dim + 1)) // 2

    # C is diagonal — lower-tri = diagonal itself
    C = sparse.csr_matrix(C_dense)

    p = 0  # equality column counter
    l = 0  # inequality column counter

    _i_ineq = []; _j_ineq = []; _d_ineq = []
    _i_eq   = []; _j_eq   = []; _d_eq   = []

    _Bt = []   # inequality COO blocks (packed_dim x step each)
    _u  = []   # inequality RHS values

    _At = []   # equality COO blocks
    _b  = []   # equality RHS values

    cuts_classes = []

    for constr, b in tqdm(constrs[1:], desc='M+ Lifting', unit='constr'):
        for i in range(1, n + 1):

            if skip_func and skip_func(constr, i):
                continue

            # (x_i)(a_j x - b) <= 0  lifts to  <A_cut, X> >= 0
            # Canonical A entries: A_{i,var} = constr[var]/2 for var != i (off-diag),
            #                      A_{i,i} accumulated with constr[i] (if in constr) + (-b)
            for var in constr:
                idx = mat_idx(i, var)
                _i_ineq.append(idx)
                _j_ineq.append(l)
                _d_ineq.append(constr[var] if i == var else 0.5 * constr[var])

            idx = mat_idx(i, i)
            _i_ineq.append(idx)
            _j_ineq.append(l)
            _d_ineq.append(-b)
            _u.append(0.)
            cuts_classes.append(1.)
            l += 1
            if l % step == 0:
                _i_ineq, _j_ineq, _d_ineq, l = push(
                    _i_ineq, _j_ineq, _d_ineq, _Bt, step, packed_dim)

            # (1 - x_i)(a_j x - b) <= 0  lifts to  <A_cut, X> >= rhs
            for var in constr:
                idx = mat_idx(i, var)
                _i_ineq.append(idx)
                _j_ineq.append(l)
                _d_ineq.append(-constr[var] if i == var else -0.5 * constr[var])
                idx = mat_idx(var, var)
                _i_ineq.append(idx)
                _j_ineq.append(l)
                _d_ineq.append(constr[var])

            idx = mat_idx(i, i)
            _i_ineq.append(idx)
            _j_ineq.append(l)
            _d_ineq.append(b)
            _u.append(b)
            cuts_classes.append(2.)
            l += 1
            if l % step == 0:
                _i_ineq, _j_ineq, _d_ineq, l = push(
                    _i_ineq, _j_ineq, _d_ineq, _Bt, step, packed_dim)

    if lift_bounds:
        for i, j in itertools.combinations(range(1, n + 1), 2):
            # X_{ij} - X_{ii} <= 0  lifts to  <A, X> >= 0
            # Canonical: A_{i,i} = -1, A_{i,j} = 0.5 (off-diagonal, i > j assumed)
            ii, jj = max(i, j), min(i, j)
            idx = mat_idx(ii, ii)
            _i_ineq.append(idx); _j_ineq.append(l); _d_ineq.append(-1.)
            idx = mat_idx(ii, jj)
            _i_ineq.append(idx); _j_ineq.append(l); _d_ineq.append(0.5)
            _u.append(0.)
            cuts_classes.append(3.)
            l += 1
            if l % step == 0:
                _i_ineq, _j_ineq, _d_ineq, l = push(
                    _i_ineq, _j_ineq, _d_ineq, _Bt, step, packed_dim)

            # X_{ii} + X_{jj} - X_{ij} >= 1
            # Canonical: A_{i,i} = 1, A_{j,j} = 1, A_{i,j} = -0.5
            idx = mat_idx(ii, ii)
            _i_ineq.append(idx); _j_ineq.append(l); _d_ineq.append(1.)
            idx = mat_idx(jj, jj)
            _i_ineq.append(idx); _j_ineq.append(l); _d_ineq.append(1.)
            idx = mat_idx(ii, jj)
            _i_ineq.append(idx); _j_ineq.append(l); _d_ineq.append(-0.5)
            _u.append(1.)
            cuts_classes.append(4.)
            l += 1
            if l % step == 0:
                _i_ineq, _j_ineq, _d_ineq, l = push(
                    _i_ineq, _j_ineq, _d_ineq, _Bt, step, packed_dim)

    if _i_ineq:
        _i_ineq, _j_ineq, _d_ineq, l = push(
            _i_ineq, _j_ineq, _d_ineq, _Bt, l % step, packed_dim)

    # Equality constraints ---------------------------------------------------
    # X_{00} = 1
    _i_eq.append(0); _j_eq.append(p); _d_eq.append(1.)
    _b.append(1.)
    p += 1

    # X_{0,i} = X_{i,i} for i = 1..n
    for i in range(1, n + 1):
        idx = mat_idx(i, 0)   # off-diagonal: canonical A_{i,0} = 0.5
        _i_eq.append(idx); _j_eq.append(p); _d_eq.append(0.5)
        idx = mat_idx(i, i)
        _i_eq.append(idx); _j_eq.append(p); _d_eq.append(-1.)
        _b.append(0.)
        p += 1
        if p % step == 0:
            _i_eq, _j_eq, _d_eq, p = push(
                _i_eq, _j_eq, _d_eq, _At, step, packed_dim)

    if _i_eq:
        _i_eq, _j_eq, _d_eq, p = push(
            _i_eq, _j_eq, _d_eq, _At, p % step, packed_dim)

    # Assemble SDPModel.A = (num_rows x packed_dim), equalities first
    A_eq   = sparse.hstack(_At, format='csc').T.tocsr()   # (n+1) x packed_dim
    A_ineq = sparse.hstack(_Bt, format='csc').T.tocsr()   # m_ineq x packed_dim
    A      = sparse.vstack([A_eq, A_ineq], format='csr')

    m_eq   = A_eq.shape[0]
    m_ineq = A_ineq.shape[0]
    num_rows = m_eq + m_ineq

    row_rhs    = np.concatenate([np.array(_b), np.array(_u)])
    row_senses = np.array(['='] * m_eq + ['>'] * m_ineq, dtype='U1')

    return SDPModel(
        dim=dim, n=n, num_rows=num_rows,
        C=C, A=A,
        row_rhs=row_rhs, row_senses=row_senses,
        cut_classes=np.array(cuts_classes),
    )


def m_plus_lifting(filename, model_out_dir='', work_dir='', step=100000,
                   lift_bounds=True, skip_func=None, do_split=False):
    """Apply the Lovász-Schrijver M+ operator to an LP and save the SDP model.

    Parameters
    ----------
    filename     : LP file name (without directory)
    model_out_dir: directory where the .mat output is written
    work_dir     : directory containing the LP file
    step         : batch size for building the sparse constraint matrix
    lift_bounds  : whether to include 0-1 box lifting constraints
    skip_func    : optional function(constr, i) -> bool to skip constraints
    do_split     : split Bt into 4 files when True (for large models)

    Returns
    -------
    model : SDPModel in canonical lower-triangular format
    """
    start = time.time()
    print('File: ', filename)
    path_to_file = os.path.join(work_dir, filename)
    out_filename = '.'.join(filename.split('.')[:-1])

    model = _compute_m_plus_model(path_to_file, step=step,
                                   lift_bounds=lift_bounds, skip_func=skip_func)

    out_path = os.path.join(model_out_dir, out_filename)
    print('Saving MAT at: %s.mat' % out_path)
    model.write_sdpnal_mat(out_path, do_split=do_split)

    elapsed = time.time() - start
    eq_rows   = int((model.row_senses == '=').sum())
    ineq_rows = int((model.row_senses == '>').sum())
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

    Extends Theta with non-edge non-positivity inequalities:
      - <A_ij, X> >= 0  (i.e. X_{i+1,j+1} <= 0)  for each non-edge (i, j)
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
        ii, jj = max(i, j), min(i, j)
        _ib.append(mat_idx(ii + 1, jj + 1)); _jb.append(l); _db.append(-0.5)
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
                    C=C, A=A, row_rhs=row_rhs, row_senses=row_senses)


if __name__ == "__main__":
    pass
