import os, time, itertools
import numpy as np
import re
from scipy import sparse
from scipy.io import savemat
from scipy.sparse import lil_matrix, hstack

"""
    The main functions in this file provide the SDP models.
    The outputs are .mat files saved on the specified path.   

        m_plus_lifting(): Lovasz-Schrijver Lift-and-Project 
                          formulation.

             Theta_SDP(): Theta function SDP formulation.
    
    Auxiliary functions parse constraints from the .lp files.
"""


def read_index_from_lp(path_to_file):
    try:
        with open(path_to_file) as f:
            exp_id = re.compile(r"^\s*.*:")

            i = 0
            for line in f:
                i += 1
                if '\\' in line or not line:
                    # Skip comments
                    continue
                if 'bounds' in line.lower():
                    # End at bounds section
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
            j = 0
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
        # print('LP file has been read.')
        pass


def read_constr_from_lp(path_to_file):
    constr = re.compile(r"[\+\-]?\s?[\d]*\s?x\[?\d+\]?")
    rhs = re.compile(r"[\<\>\=]+\s*([\d]+)")

    lines = read_expr_from_lp(path_to_file)
    j = 0
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


def sp_unique(sp_matrix, axis=0):
    ''' Returns a sparse matrix with the unique rows (axis=0)
    or columns (axis=1) of an input sparse matrix sp_matrix'''
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

    sp_matrix = sparse.lil_matrix((nrows_uniq, ncols), dtype=dt)  #  or sp_matrix.resize(nrows_uniq, ncols)
    sp_matrix.data = data
    sp_matrix.rows = rows

    ret = sp_matrix.asformat(old_format)
    if axis == 1:
        ret = ret.T        
    return ret, ind


def mat_idx(i, j):
    return int(((i + 1) * (i)) / 2 + (j)) if i >= j else int(((j + 1) * (j)) / 2 + (i))


def push(_i, _j, _d, _Bt, step, dim):
    coo = sparse.coo_matrix((_d, (_i, _j)), shape=(dim, step))
    del _i, _j, _d
    _Bt.append(coo)
    return [], [], [], 0


def lovasz_schrijver_filter(constr, i):
    return i in constr


def map_edge_constraints(filename):
    constrs = list(read_constr_from_lp(filename))
    n = len(constrs[0])
    map = {}
    idx = 1
    nod = 0
    for constr, b in constrs[1:]:
        for j in range(1, n + 1):

            if lovasz_schrijver_filter(constr, j):
                continue

            # k != i,j
            map[idx] = 'lift_edge1:1'
            idx += 1

            # k != i,j
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
    nod = 0
    for constr, b in constrs[1:]:
        for j in range(1, n + 1):
            if j in constr:
                # j \in C
                map[idx] = 'lift_clq1:1'
            else:
                # j \notin C
                map[idx] = 'lift_clq2:1'

            idx += 1

            if j in constr:
                # j \in C
                map[idx] = 'lift_clq1:2'
            else:
                # j \notin C
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
    nod = 0
    for constr, b in constrs[1:]:
        for i in constr:
            if constr[i] != 1.:
                nod = i
        for j in range(1, n + 1):
            if j == nod:
                # i \in V
                map[idx] = 'lift_nod1:1'
            elif j in constr:
                # j \in \Gamma(i)
                map[idx] = 'lift_nod2:1'
            else:
                # j \notin \Gamma(i)
                map[idx] = 'lift_nod3:1'
            idx += 1

            if j == nod:
                # i \in V
                map[idx] = 'lift_nod1:2'
            elif j in constr:
                # j \in \Gamma(i)
                map[idx] = 'lift_nod2:2'
            else:
                # j \notin \Gamma(i)
                map[idx] = 'lift_nod3:2'

            idx += 1

    for _ in itertools.combinations(range(1, n + 1), 2):
        map[idx] = 'bound2'
        idx += 1
        map[idx] = 'bound4'
        idx += 1
        
    return map


def m_plus_lifting(filename, model_out_dir='', work_dir='', step=100000, lift_bounds=True, skip_func=None, do_split=False):
    start = time.time()
    print('File: ', filename)
    path_to_file = os.path.join(work_dir, filename)
    out_filename = '.'.join(filename.split('.')[:-1])
    constrs = list(read_constr_from_lp(path_to_file))

    obj = constrs[0]
    C = -np.diag([0.] + list(obj.values()))
    del obj
    n = C.shape[0] - 1
    dim = int(((n + 1) * (n + 2)) / 2)

    _2 = np.round(np.sqrt(2)/2, 10)

    p = 0  # Equalities
    l = 0  # Inequalities

    _i = []
    _j = []
    _d = []

    _At = []  # Equalities
    _b = []

    _Bt = []  # Inequalities
    _u = []
    
    cuts_classes = []

    print('M+ Lifting...')
    for constr, b in constrs[1:]:
        for i in range(1, n + 1):

            if skip_func and skip_func(constr, i):
                continue

            # (x_i)(a_jx - b) <= 0
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

            # (1 - x_i)(a_jx - b) <= 0
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

    u = np.matrix(_u).T
    b = np.matrix(_b).T

    print('Saving MAT at: %s' % os.path.join(model_out_dir, out_filename) + '.mat')
    # Exporting the model to output file
    if not do_split:
        savemat(os.path.join(model_out_dir, out_filename) + '.mat', \
                {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'L': 0., 'C' : C, 's' : float(n + 1), 'cut_classes' : np.array(cuts_classes)})
    else:
        savemat(os.path.join(model_out_dir, out_filename) + '_1.mat', \
                {'At' : At, 'b' : b, 'u' : u, 'L': 0., 'C' : C, 's' : float(n + 1), 'cut_classes' : np.array(cuts_classes)})
        splt = int(np.ceil(Bt.shape[1]/4))
        # print((Bt[:, :splt].shape[1] + Bt[:, splt:2*splt].shape[1] + Bt[:, 2*splt:3*splt].shape[1] + Bt[:, 3*splt:].shape[1]))
        # print(Bt.shape[1])
        savemat(os.path.join(model_out_dir, out_filename) + '_2.mat', {'Bt_1' : Bt[:, :splt]})
        savemat(os.path.join(model_out_dir, out_filename) + '_3.mat', {'Bt_2' : Bt[:, splt:2*splt]})
        savemat(os.path.join(model_out_dir, out_filename) + '_4.mat', {'Bt_3' : Bt[:, 2*splt:3*splt]})
        savemat(os.path.join(model_out_dir, out_filename) + '_5.mat', {'Bt_4' : Bt[:, 3*splt:]})
    end = time.time() - start
    print('Finished! Total time: %.2f' % end)
    print('Order of the Matrix Variable: %d' % (n + 1))
    print('Number of Equalities: %d' % At.shape[1])
    print('Number of Inequalities: %d' % Bt.shape[1])
    print('-'*12)


def Theta_SDP(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_SDP')
        print('Reading the model...')

    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation 
        n = len(G.nodes())
        dim = n + 1
        m = len(G.edges())
        
        l = 0
        mleq = 0

        # Check in advance if the problem size will exceed the limit
        if limit and l > limit:
            if debug:
                print("Problem is too big: %d rows" % l)
            return

        # Objective function (minimization problem)
        C = lil_matrix((dim, dim))
        for i in range(n):
            C[i + 1, i + 1] = -1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        # x_ij = 0 , for all (i, j) in E
        for i, j in G.edges():
            i1 = (i + 1)*dim + (j + 1) 
            i2 = (j + 1)*dim + (i + 1)
            _A[i1, p] = .5
            _A[i2, p] = .5
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # Add constraint: X_00 = 1
        _A[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        if p == step:
            A = hstack([A, _A])
            b = np.vstack([b, _b])
            _A = lil_matrix((dim**2, step))
            _b = np.zeros((step, 1))
            p = 0

        # Add constraint: X_ii = X_0i = X_i0 for i in V
        for i in G.nodes():
            i1 = (i + 1)
            i2 = (i + 1)*dim
            _A[i1, p] = .5
            _A[i2, p] = .5
            i1 = (i + 1)*dim + (i + 1)
            _A[i1, p] = -1.
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0
        if p > 0:
            A = hstack([A,  _A.tocsc()[:, :p]]);
            b = np.vstack([b, _b[:p]])
        A = A.transpose()

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % dim)
            print('Number of constraints: %d' % l)
        
        savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')

    elif model_out.lower() == 'sdpnal': # SDPNAL model creation

        n = len(G.nodes())
        m = len(G.edges())

        dim = int(((n + 1) * (n + 2))/2)
        
        l = 0

        _2 = np.sqrt(2)

        # Objective function (minimization problem)
        C = lil_matrix((n + 1, n + 1))
        for i in range(n):
            C[i + 1, i + 1] = -1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        # x_ij == 0
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 2)*(ii + 1))/2 + jj + 1 
            _At[idx, p] = _2*0.5
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # Add constraint: X_00 = 1
        _At[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        if p == step:
            At = hstack([At, _At])
            b = np.vstack([b, _b])
            _At = lil_matrix((dim, step))
            _b = np.zeros((step, 1))
            p = 0

        # Add constraint: X_ii = X_0i = X_i0 for i in V
        for i in G.nodes():
            idx = ((i + 2)*(i + 1))/2 
            _At[idx, p] = _2*0.5
            
            idx = ((i + 2)*(i + 1))/2 + (i + 1)
            _At[idx, p] = -1.
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0
        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : C, 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_plus_SDP(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_plus_SDP2')
        print('Reading the model...')

    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation
        n = len(G.nodes())
        dim = n + 1
        m = len(G.edges())
        G_complement_edges = [i for i in itertools.combinations(G.nodes(), 2) if i not in G.edges()]
        
        l = 0
        mleq = 0

        # Check in advance if the problem size will exceed the limit
        if limit and l > limit:
            if debug:
                print("Problem is too big: %d rows" % l)
            return

        # Objective function (minimization problem)
        C = lil_matrix((dim, dim))
        for i in range(n):
            C[i + 1, i + 1] = -1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0


        # x_ij > 0 , for all (i, j) not in E
        for i, j in G_complement_edges:
            # proceeding column wise
            i1 = (i + 1)*dim + (j + 1) 
            i2 = (j + 1)*dim + (i + 1)
            _A[i1, p] = -.5
            _A[i2, p] = -.5
            mleq += 1
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # x_ij = 0 , for all (i, j) in E
        for i, j in G.edges():
            i1 = (i + 1)*dim + (j + 1) 
            i2 = (j + 1)*dim + (i + 1)
            _A[i1, p] = .5
            _A[i2, p] = .5
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # Add constraint: X_00 = 1
        _A[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        if p == step:
            A = hstack([A, _A])
            b = np.vstack([b, _b])
            _A = lil_matrix((dim**2, step))
            _b = np.zeros((step, 1))
            p = 0

        # Add constraint: X_ii = X_0i = X_i0 for i in V
        for i in G.nodes():
            i1 = (i + 1)
            i2 = (i + 1)*dim
            _A[i1, p] = .5
            _A[i2, p] = .5
            i1 = (i + 1)*dim + (i + 1)
            _A[i1, p] = -1.
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0
        if p > 0:
            A = hstack([A,  _A.tocsc()[:, :p]]);
            b = np.vstack([b, _b[:p]])
        A = A.transpose()

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % dim)
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')

    elif model_out.lower() == 'sdpnal': # SDPNAL model creation

        n = len(G.nodes())
        m = len(G.edges())

        dim = int(((n + 1) * (n + 2))/2)
        
        l = 0

        _2 = np.sqrt(2)

        # Objective function (minimization problem)
        C = lil_matrix((n + 1, n + 1))
        for i in range(n):
            C[i + 1, i + 1] = -1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        # x_ij == 0
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 2)*(ii + 1))/2 + jj + 1 
            idx = int(idx)
            _At[idx, p] = _2*0.5
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # Add constraint: X_00 = 1
        _At[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        if p == step:
            At = hstack([At, _At])
            b = np.vstack([b, _b])
            _At = lil_matrix((dim, step))
            _b = np.zeros((step, 1))
            p = 0

        # Add constraint: X_ii = X_0i = X_i0 for i in V
        for i in G.nodes():
            idx = ((i + 2)*(i + 1))/2 
            idx = int(idx)
            _At[idx, p] = _2*0.5
            
            idx = ((i + 2)*(i + 1))/2 + (i + 1)
            idx = int(idx)
            _At[idx, p] = -1.
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0
        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : C, 'L' : 0., 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


if __name__ == "__main__":
    pass
    # filepath = "/Users/feb223/projects/SDP_MSSP_GCP/StableSets/DIMACS/lp/brock200_1_nod_alpha.lp"
    # maps = map_nod_constraints(filepath)
    # h = loadmat("/Users/feb223/projects/SDP_MSSP_GCP/StableSets/DIMACS/violated_cuts/brock200_1_nod_alpha_viol_test.mat")
    # # m_plus_lifting('G_300_0.1_5_edge.lp', skip_func=lovasz_schrijver_filer, lift_bounds=False)
    # # maps = map_nod_constraints('/home/federico/adal/StableSets/Random-np/300/G_300_0.9_1_nod_theta.lp')
    # # h = loadmat('/home/federico/adal/StableSets/Random-np/violated_cuts/G_300_0.9_1_nod_theta_viol.mat')
    # # print(h['added_cuts_idx'])
    # print(len(maps))
    # print(count_values(maps, h['added_cuts_idx']))

    # filepath = "/Users/feb223/projects/SDP_MSSP_GCP/StableSets/DIMACS/lp/brock200_1_cov.lp"
    # maps = map_clique_constraints(filepath)
    # h = loadmat("/Users/feb223/projects/SDP_MSSP_GCP/StableSets/DIMACS/violated_cuts/brock200_1_cov_viol_test.mat")
    # print(len(maps))
    # print(count_values(maps, h['added_cuts_idx']))

# maps = map_nod_constraints('/home/federico/adal/StableSets/Random-np/300/G_300_0.9_1_nod_theta.lp')
# h = loadmat('/home/federico/adal/StableSets/Random-np/violated_cuts/G_300_0.9_1_nod_theta_viol_test.mat')
# # print(h['added_cuts_idx'])
# print(len(maps))
# print(count_values(maps, h['added_cuts_idx']))

