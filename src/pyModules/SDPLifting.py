# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# SDP model generation via the Lovász-Schrijver M+ lift-and-project operator
# and Theta/Theta+ SDP formulations.

import os, time, itertools
import numpy as np
import re
from scipy import sparse
from scipy.io import savemat
from scipy.sparse import lil_matrix, hstack
from tqdm import tqdm

"""
    The main functions in this file provide the SDP models.
    The outputs are .mat files saved on the specified path.

        m_plus_lifting(): Lovász-Schrijver Lift-and-Project formulation.

             Theta_SDP(): Theta function SDP formulation.

        Theta_plus_SDP(): Theta+ function SDP formulation.

    Format-specific serialization is separated from the computation:

        write_sdpnal_mat(): Serialize a model dict to SDPNAL+ .mat format.

    This separation allows targeting a different SDP solver by implementing
    an alternative write function without touching the lifting logic.
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
    """Build the M+ SDP model matrices from an LP file.

    Returns a solver-agnostic dict with keys:
        At, b      : equality constraint matrix and RHS
        Bt, u      : inequality constraint matrix and RHS
        C          : objective matrix (negated node weights)
        cut_classes: class label (1-4) for each inequality
        n          : number of graph nodes
    """
    constrs = list(read_constr_from_lp(path_to_file))

    obj = constrs[0]
    C = -np.diag([0.] + list(obj.values()))
    n = C.shape[0] - 1
    dim = ((n + 1) * (n + 2)) // 2

    _2 = np.round(np.sqrt(2) / 2, 10)

    p = 0  # Equality counter
    l = 0  # Inequality counter

    _i = []
    _j = []
    _d = []

    _At = []
    _b = []

    _Bt = []
    _u = []

    cuts_classes = []

    for constr, b in tqdm(constrs[1:], desc='M+ Lifting', unit='constr'):
        for i in range(1, n + 1):

            if skip_func and skip_func(constr, i):
                continue

            # (x_i)(a_j x - b) <= 0
            for var in constr:
                idx = mat_idx(i, var)
                _i.append(idx)
                _j.append(l)
                _d.append(constr[var] if i == var else _2 * constr[var])

            idx = mat_idx(i, i)
            _i.append(idx)
            _j.append(l)
            _d.append(-b)
            _u.append(0.)
            cuts_classes.append(1.)
            l += 1
            if l % step == 0:
                _i, _j, _d, l = push(_i, _j, _d, _Bt, step, dim)

            # (1 - x_i)(a_j x - b) <= 0
            for var in constr:
                idx = mat_idx(i, var)
                _i.append(idx)
                _j.append(l)
                _d.append(-constr[var] if i == var else -_2 * constr[var])
                idx = mat_idx(var, var)
                _i.append(idx)
                _j.append(l)
                _d.append(constr[var])

            idx = mat_idx(i, i)
            _i.append(idx)
            _j.append(l)
            _d.append(b)
            _u.append(b)
            cuts_classes.append(2.)
            l += 1
            if l % step == 0:
                _i, _j, _d, l = push(_i, _j, _d, _Bt, step, dim)

    if lift_bounds:
        for i, j in itertools.combinations(range(1, n + 1), 2):
            idx = mat_idx(i, i)
            _i.append(idx)
            _j.append(l)
            _d.append(-1.)
            idx = mat_idx(i, j)
            _i.append(idx)
            _j.append(l)
            _d.append(_2)
            _u.append(0.)
            cuts_classes.append(3.)
            l += 1
            if l % step == 0:
                _i, _j, _d, l = push(_i, _j, _d, _Bt, step, dim)

            idx = mat_idx(i, i)
            _i.append(idx)
            _j.append(l)
            _d.append(1.)
            idx = mat_idx(j, j)
            _i.append(idx)
            _j.append(l)
            _d.append(1.)
            idx = mat_idx(i, j)
            _i.append(idx)
            _j.append(l)
            _d.append(-_2)
            _u.append(1.)
            cuts_classes.append(4.)
            l += 1
            if l % step == 0:
                _i, _j, _d, l = push(_i, _j, _d, _Bt, step, dim)

    if _i:
        _i, _j, _d, l = push(_i, _j, _d, _Bt, l % step, dim)

    _i.append(0)
    _j.append(p)
    _d.append(1.)
    _b.append(1.)
    p += 1

    for i in range(1, n + 1):
        idx = mat_idx(i, 0)
        _i.append(idx)
        _j.append(p)
        _d.append(_2)

        idx = mat_idx(i, i)
        _i.append(idx)
        _j.append(p)
        _d.append(-1.)

        _b.append(0.)
        p += 1
        if p % step == 0:
            _i, _j, _d, p = push(_i, _j, _d, _At, step, dim)

    if _i:
        _i, _j, _d, p = push(_i, _j, _d, _At, p % step, dim)

    At = sparse.hstack(_At, format='csc')
    Bt = sparse.hstack(_Bt, format='csc')
    u = np.array(_u).reshape(-1, 1)
    b = np.array(_b).reshape(-1, 1)

    return {'At': At, 'Bt': Bt, 'b': b, 'u': u, 'C': C,
            'cut_classes': np.array(cuts_classes), 'n': n}


def write_sdpnal_mat(model, out_path, do_split=False):
    """Serialize a model dict to SDPNAL+ .mat format.

    Parameters
    ----------
    model    : dict returned by _compute_m_plus_model
    out_path : output file path without extension (e.g. '.../graphname_edge')
    do_split : if True, split Bt across 4 files (for very large models)
    """
    At = model['At']
    Bt = model['Bt']
    b  = model['b']
    u  = model['u']
    C  = model['C']
    cut_classes = model['cut_classes']
    n  = model['n']

    if not do_split:
        savemat(out_path + '.mat',
                {'At': At, 'b': b, 'Bt': Bt, 'u': u, 'L': 0.,
                 'C': C, 's': float(n + 1), 'cut_classes': cut_classes})
    else:
        savemat(out_path + '_1.mat',
                {'At': At, 'b': b, 'u': u, 'L': 0.,
                 'C': C, 's': float(n + 1), 'cut_classes': cut_classes})
        splt = int(np.ceil(Bt.shape[1] / 4))
        savemat(out_path + '_2.mat', {'Bt_1': Bt[:, :splt]})
        savemat(out_path + '_3.mat', {'Bt_2': Bt[:, splt:2*splt]})
        savemat(out_path + '_4.mat', {'Bt_3': Bt[:, 2*splt:3*splt]})
        savemat(out_path + '_5.mat', {'Bt_4': Bt[:, 3*splt:]})


def m_plus_lifting(filename, model_out_dir='', work_dir='', step=100000,
                   lift_bounds=True, skip_func=None, do_split=False):
    """Apply the Lovász-Schrijver M+ operator to an LP and save the SDP model.

    Parameters
    ----------
    filename     : LP file name (without directory)
    model_out_dir: directory where the .mat output is written
    work_dir     : directory containing the LP file
    step         : batch size for building the sparse Bt matrix
    lift_bounds  : whether to include 0-1 box lifting constraints
    skip_func    : optional function(constr, i) -> bool to skip constraints
    do_split     : split Bt into 4 pieces when True (for large models)

    Returns
    -------
    model : solver-agnostic dict (At, Bt, b, u, C, cut_classes, n)
    """
    start = time.time()
    print('File: ', filename)
    path_to_file = os.path.join(work_dir, filename)
    out_filename = '.'.join(filename.split('.')[:-1])

    model = _compute_m_plus_model(path_to_file, step=step,
                                   lift_bounds=lift_bounds, skip_func=skip_func)

    out_path = os.path.join(model_out_dir, out_filename)
    print('Saving MAT at: %s.mat' % out_path)
    write_sdpnal_mat(model, out_path, do_split=do_split)

    elapsed = time.time() - start
    print('Finished! Total time: %.2f' % elapsed)
    print('Order of the Matrix Variable: %d' % (model['n'] + 1))
    print('Number of Equalities: %d' % model['At'].shape[1])
    print('Number of Inequalities: %d' % model['Bt'].shape[1])
    print('-' * 12)

    return model


# ---------------------------------------------------------------------------
# Theta SDP
# ---------------------------------------------------------------------------

def _compute_theta_adal(G, step=10000):
    """Build Theta SDP matrices in ADAL (full-matrix) format."""
    n = len(G.nodes())
    dim = n + 1

    C = lil_matrix((dim, dim))
    for i in range(n):
        C[i + 1, i + 1] = -1.

    A = lil_matrix((dim**2, 0))
    b = np.zeros((0, 1))
    _A = lil_matrix((dim**2, step))
    _b = np.zeros((step, 1))
    p = 0

    def _flush():
        nonlocal A, b, _A, _b, p
        A = hstack([A, _A])
        b = np.vstack([b, _b])
        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

    for i, j in G.edges():
        i1 = (i + 1)*dim + (j + 1)
        i2 = (j + 1)*dim + (i + 1)
        _A[i1, p] = .5
        _A[i2, p] = .5
        p += 1
        if p == step:
            _flush()

    _A[0, p] = 1.
    _b[p] = 1.
    p += 1
    if p == step:
        _flush()

    for i in G.nodes():
        i1 = (i + 1)
        i2 = (i + 1)*dim
        _A[i1, p] = .5
        _A[i2, p] = .5
        _A[(i + 1)*dim + (i + 1), p] = -1.
        p += 1
        if p == step:
            _flush()

    if p > 0:
        A = hstack([A, _A.tocsc()[:, :p]])
        b = np.vstack([b, _b[:p]])

    return {'A': A.T, 'b': b, 'C': C, 'mleq': 0, 'n': n}


def _compute_theta_sdpnal(G, step=10000):
    """Build Theta SDP matrices in SDPNAL+ (packed symmetric) format."""
    n = len(G.nodes())
    dim = ((n + 1) * (n + 2)) // 2
    _2 = np.sqrt(2)

    C = lil_matrix((n + 1, n + 1))
    for i in range(n):
        C[i + 1, i + 1] = -1.

    At = lil_matrix((dim, 0))
    b = np.zeros((0, 1))
    _At = lil_matrix((dim, step))
    _b = np.zeros((step, 1))
    p = 0

    def _flush():
        nonlocal At, b, _At, _b, p
        At = hstack([At, _At])
        b = np.vstack([b, _b])
        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

    for i, j in G.edges():
        ii, jj = sorted([i, j], reverse=True)
        idx = ((ii + 2) * (ii + 1)) // 2 + jj + 1
        _At[idx, p] = _2 * 0.5
        p += 1
        if p == step:
            _flush()

    _At[0, p] = 1.
    _b[p] = 1.
    p += 1
    if p == step:
        _flush()

    for i in G.nodes():
        idx = ((i + 2) * (i + 1)) // 2
        _At[idx, p] = _2 * 0.5
        idx2 = ((i + 2) * (i + 1)) // 2 + (i + 1)
        _At[idx2, p] = -1.
        p += 1
        if p == step:
            _flush()

    if p > 0:
        At = hstack([At, _At.tocsc()[:, :p]])
        b = np.vstack([b, _b[:p]])

    return {'At': At, 'b': b, 'C': C, 'n': n}


def Theta_SDP(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    """Compute and save the Theta SDP for graph G.

    Parameters
    ----------
    model_out : 'adal' or 'sdpnal' — selects the output matrix format
    """
    assert model_out.lower() in {'adal', 'sdpnal'}, 'Supported models: adal, sdpnal'
    if debug:
        print('Theta_SDP')

    start_time = time.time()

    if model_out.lower() == 'adal':
        model = _compute_theta_adal(G, step=step)
        savemat(os.path.join(model_out_dir, filename) + '.mat',
                {'A': model['A'], 'b': model['b'], 'C': model['C'], 'mleq': model['mleq']})
    else:
        model = _compute_theta_sdpnal(G, step=step)
        savemat(os.path.join(model_out_dir, filename) + '.mat',
                {'At': model['At'], 'b': model['b'], 'C': model['C'], 's': float(model['n'] + 1)})

    if debug:
        elapsed = time.time() - start_time
        print('Finished! Time elapsed: %.2f' % elapsed)
        print('Dimension of matrix variable: %d' % (model['n'] + 1))
        print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


# ---------------------------------------------------------------------------
# Theta+ SDP
# ---------------------------------------------------------------------------

def _compute_theta_plus_adal(G, step=10000):
    """Build Theta+ SDP matrices in ADAL format (adds non-edge positivity constraints)."""
    n = len(G.nodes())
    dim = n + 1
    G_complement_edges = [e for e in itertools.combinations(G.nodes(), 2) if e not in G.edges()]
    mleq = 0

    C = lil_matrix((dim, dim))
    for i in range(n):
        C[i + 1, i + 1] = -1.

    A = lil_matrix((dim**2, 0))
    b = np.zeros((0, 1))
    _A = lil_matrix((dim**2, step))
    _b = np.zeros((step, 1))
    p = 0

    def _flush():
        nonlocal A, b, _A, _b, p
        A = hstack([A, _A])
        b = np.vstack([b, _b])
        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

    for i, j in G_complement_edges:
        i1 = (i + 1)*dim + (j + 1)
        i2 = (j + 1)*dim + (i + 1)
        _A[i1, p] = -.5
        _A[i2, p] = -.5
        mleq += 1
        p += 1
        if p == step:
            _flush()

    for i, j in G.edges():
        i1 = (i + 1)*dim + (j + 1)
        i2 = (j + 1)*dim + (i + 1)
        _A[i1, p] = .5
        _A[i2, p] = .5
        p += 1
        if p == step:
            _flush()

    _A[0, p] = 1.
    _b[p] = 1.
    p += 1
    if p == step:
        _flush()

    for i in G.nodes():
        i1 = (i + 1)
        i2 = (i + 1)*dim
        _A[i1, p] = .5
        _A[i2, p] = .5
        _A[(i + 1)*dim + (i + 1), p] = -1.
        p += 1
        if p == step:
            _flush()

    if p > 0:
        A = hstack([A, _A.tocsc()[:, :p]])
        b = np.vstack([b, _b[:p]])

    return {'A': A.T, 'b': b, 'C': C, 'mleq': mleq, 'n': n}


def _compute_theta_plus_sdpnal(G, step=10000):
    """Build Theta+ SDP matrices in SDPNAL+ format."""
    n = len(G.nodes())
    dim = ((n + 1) * (n + 2)) // 2
    _2 = np.sqrt(2)

    C = lil_matrix((n + 1, n + 1))
    for i in range(n):
        C[i + 1, i + 1] = -1.

    At = lil_matrix((dim, 0))
    b = np.zeros((0, 1))
    _At = lil_matrix((dim, step))
    _b = np.zeros((step, 1))
    p = 0

    def _flush():
        nonlocal At, b, _At, _b, p
        At = hstack([At, _At])
        b = np.vstack([b, _b])
        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

    for i, j in G.edges():
        ii, jj = sorted([i, j], reverse=True)
        idx = ((ii + 2) * (ii + 1)) // 2 + jj + 1
        _At[idx, p] = _2 * 0.5
        p += 1
        if p == step:
            _flush()

    _At[0, p] = 1.
    _b[p] = 1.
    p += 1
    if p == step:
        _flush()

    for i in G.nodes():
        idx = ((i + 2) * (i + 1)) // 2
        _At[idx, p] = _2 * 0.5
        idx2 = ((i + 2) * (i + 1)) // 2 + (i + 1)
        _At[idx2, p] = -1.
        p += 1
        if p == step:
            _flush()

    if p > 0:
        At = hstack([At, _At.tocsc()[:, :p]])
        b = np.vstack([b, _b[:p]])

    return {'At': At, 'b': b, 'C': C, 'n': n}


def Theta_plus_SDP(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    """Compute and save the Theta+ SDP for graph G.

    Parameters
    ----------
    model_out : 'adal' or 'sdpnal' — selects the output matrix format
    """
    assert model_out.lower() in {'adal', 'sdpnal'}, 'Supported models: adal, sdpnal'
    if debug:
        print('Theta_plus_SDP')

    start_time = time.time()

    if model_out.lower() == 'adal':
        model = _compute_theta_plus_adal(G, step=step)
        savemat(os.path.join(model_out_dir, filename) + '.mat',
                {'A': model['A'], 'b': model['b'], 'C': model['C'], 'mleq': model['mleq']})
    else:
        model = _compute_theta_plus_sdpnal(G, step=step)
        savemat(os.path.join(model_out_dir, filename) + '.mat',
                {'At': model['At'], 'b': model['b'], 'C': model['C'],
                 'L': 0., 's': float(model['n'] + 1)})

    if debug:
        elapsed = time.time() - start_time
        print('Finished! Time elapsed: %.2f' % elapsed)
        print('Dimension of matrix variable: %d' % (model['n'] + 1))
        print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


if __name__ == "__main__":
    pass
