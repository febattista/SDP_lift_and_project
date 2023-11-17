import os
import re
import numpy as np
from scipy import sparse
from scipy.io import savemat, loadmat
import time
import itertools


def mat_idx(i, j):
    return int(((i + 1) * (i)) / 2 + (j)) if i >= j else int(((j + 1) * (j)) / 2 + (i))


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
        print('LP file has been read.')


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


# %%
def push(_i, _j, _d, _Bt, step, dim):
    coo = sparse.coo_matrix((_d, (_i, _j)), shape=(dim, step))
    del _i, _j, _d
    _Bt.append(coo)
    return [], [], [], 0


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
                print('lift_nod1:1 ', idx)
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


def lovasz_schrijver_filer(constr, i):
    return i in constr


# %%
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
        print((Bt[:, :splt].shape[1] + Bt[:, splt:2*splt].shape[1] + Bt[:, 2*splt:3*splt].shape[1] + Bt[:, 3*splt:].shape[1]))
        print(Bt.shape[1])
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


def count_values(input_dict, idx):
    output_dict = {}
    for i in idx:
        value = input_dict[i[0]]
        if value not in output_dict:
            output_dict[value] = 1
        else:
            output_dict[value] += 1
    return output_dict


# m_plus_lifting('G_300_0.1_5_edge.lp', skip_func=lovasz_schrijver_filer, lift_bounds=False)
# maps = map_nod_constraints('/home/federico/adal/StableSets/Random-np/300/G_300_0.9_1_nod_theta.lp')
# h = loadmat('/home/federico/adal/StableSets/Random-np/violated_cuts/G_300_0.9_1_nod_theta_viol.mat')
# # print(h['added_cuts_idx'])
# print(len(maps))
# print(count_values(maps, h['added_cuts_idx']))

# maps = map_nod_constraints('/home/federico/adal/StableSets/Random-np/300/G_300_0.9_1_nod_theta.lp')
# h = loadmat('/home/federico/adal/StableSets/Random-np/violated_cuts/G_300_0.9_1_nod_theta_viol_test.mat')
# # print(h['added_cuts_idx'])
# print(len(maps))
# print(count_values(maps, h['added_cuts_idx']))

