#coding:utf-8
from graph import *
import time
import pandas as pd
#import gurobipy as gb
import numpy as np
import networkx as nx
#import mosek
import random
import math
from scipy.sparse import lil_matrix, vstack, hstack, save_npz, csc_matrix, lil_matrix
from scipy.io import savemat
import sys
import itertools
import subprocess
import os

def mat_idx(i, j):
    return int(((i + 1)*(i))/2 + (j) if i >= j else ((j + 1)*(j))/2 + (i))

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

    sp_matrix = lil_matrix((nrows_uniq, ncols), dtype=dt)  #  or sp_matrix.resize(nrows_uniq, ncols)
    sp_matrix.data = data
    sp_matrix.rows = rows

    ret = sp_matrix.asformat(old_format)
    if axis == 1:
        ret = ret.T        
    return ret, ind

def iter_sample_fast(iterable, samplesize):
    results = []
    # Fill in the first samplesize elements:
    try:
        for _ in range(samplesize):
            results.append(next(iterable))
    except StopIteration:
        raise ValueError("Sample larger than population.")
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterable, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items
    return results

def GCP(G):
    gcp_G = gb.Model()
    x = gcp_G.addVars(G.nodes(), range(len(G.nodes())), vtype=gb.GRB.BINARY, name='x')
    w = gcp_G.addVars(range(len(G.nodes())), vtype=gb.GRB.BINARY, name='w')
    gcp_G.setObjective(w.sum(), gb.GRB.MINIMIZE)
    gcp_G.addConstrs(x.sum(i, '*') == 1 for i in G.nodes())
    for i, k in G.edges():
        for j in range(len(G.nodes())):
            gcp_G.addConstr(x[i, j] + x[k, j] <= w[j])
    gcp_G.update()

    return gcp_G


def REP(G):
    G_c = nx.complement(G)
    mapping = {}
    i = 0
    for u in G.nodes():
        for v in itertools.chain([u], G_c.neighbors(u)):
            mapping[u, v] = i
            i += 1

    rep_G = gb.Model()
    x = rep_G.addVars(mapping.keys(), vtype=gb.GRB.BINARY, name='x')
    rep_G.update()

    c = gb.LinExpr((1., x[u,u]) for u in G.nodes())
    rep_G.setObjective(c, gb.GRB.MINIMIZE)

    rep_G.addConstrs((x.sum('*', u) >= 1 for u in G.nodes()), name='repres')
    rep_G.update()

    for u in G.nodes():
        for v,w in G.subgraph(itertools.chain(G_c.neighbors(u), [u])).edges():
            rep_G.addConstr(x[u, v] + x[u, w] <= x[u, u], name='neigh')

    rep_G.update()

    return rep_G


def Theta(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta')
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
        C = np.zeros((dim, dim))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        # z_ii - t = -1, for all i in V
        for i in G.nodes():
            i1 = (i)*dim + (i)
            i2 = (n)*dim + (n)
            _A[i1, p] = 1.
            _A[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0


        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
            _A[i1, p] = .5
            _A[i2, p] = .5
            _b[p] = -1.
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
        C = np.zeros((n + 1, n + 1))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        # z_ii - t = -1
        for i in G.nodes():
            i1 = ((i + 1)*(i))/2 + (i)
            i2 = ((n + 1)*(n))/2 + (n)
            _At[i1, p] = 1.
            _At[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            idx = int(idx)
            _At[idx, p] = _2*0.5
            _b[p] = -1.
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : C, 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_plus(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_plus')
        print('Reading the model...')

    start_time = time.time()

    G_complement_edges = [i for i in itertools.combinations(G.nodes(), 2) if i not in G.edges()]

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
        C = np.zeros((dim, dim))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        L = np.zeros((dim, dim))
        L.fill(-np.inf)
        for i, j in G_complement_edges:
        	L[i, j] = -1.
        	L[j, i] = -1.
        # z_ij >= -1 , for all (i, j) not in E
        # for i, j in G_complement_edges:
        #     i1 = (i)*dim + (j) 
        #     i2 = (j)*dim + (i)
        #     _A[i1, p] = -.5
        #     _A[i2, p] = -.5
        #     _b[p] = 1.
        #     mleq += 1
        #     l += 1
        #     p += 1
        #     if p == step:
        #         A = hstack([A, _A])
        #         b = np.vstack([b, _b])
        #         _A = lil_matrix((dim**2, step))
        #         _b = np.zeros((step, 1))
        #         p = 0

        # z_ii - t = -1, for all i in V
        for i in G.nodes():
            i1 = (i)*dim + (i)
            i2 = (n)*dim + (n)
            _A[i1, p] = 1.
            _A[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
            _A[i1, p] = .5
            _A[i2, p] = .5
            _b[p] = -1.
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq, 'L' : L})
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
        C = np.zeros((n + 1, n + 1))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        Bt = lil_matrix((dim, 0))
        u = np.zeros((0, 1))

        _Bt = lil_matrix((dim, step))
        _u = np.zeros((step, 1))
        p1 = 0

        # z_ii - t = -1
        for i in G.nodes():
            i1 = ((i + 1)*(i))/2 + (i)
            i2 = ((n + 1)*(n))/2 + (n)
            _At[i1, p] = 1.
            _At[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            idx = int(idx)
            _At[idx, p] = _2*0.5
            _b[p] = -1.
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij >= -1, for all (i, j) not in E
        L = np.zeros((n + 1, n + 1))
        L.fill(-np.inf)
        for i, j in G_complement_edges:
            L[i, j] = -1.
            L[j, i] = -1.

        # for i, j in G_complement_edges:
        #     ii, jj = sorted([i, j], reverse=True)
        #     idx = ((ii + 1)*(ii))/2 + jj
        #     idx = int(idx)
        #     _Bt[idx, p1] = -_2*0.5
        #     _u[p1] = 1.
        #     p1 += 1
        #     l += 1
        #     if p1 == step:
        #         Bt = hstack([Bt, _Bt])
        #         u = np.vstack([u, _u])
        #         _Bt = lil_matrix((dim, step))
        #         _u = np.zeros((step, 1))
        #         p1 = 0
        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        if p1 > 0:
            Bt = hstack([Bt,  _Bt.tocsc()[:, :p1]]);
            u = np.vstack([u, _u[:p1]])


        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'L': L, 'C' : C, 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_plus_triangle(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000, sample=None, seed=123456):
    # Dukanovic and Rendl
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('DukanovicRendl')
        print('Reading the model...')

    random.seed(seed)

    start_time = time.time()

    G_complement_edges = [set(i) for i in itertools.combinations(G.nodes(), 2) if i not in G.edges()]

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
        C = np.zeros((dim, dim))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        # z_ij + z_ik - z_jk - t <= -1, for ij, jk not in E [as modeled in DukanovicRendl2005 (7)]
        # for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
        #     if set((i, j)) in G_complement_edges and set((j, k)) in G_complement_edges:
        #         # z_ij
        #         i1 = (i)*dim + (j) 
        #         i2 = (j)*dim + (i)
        #         _A[i1, p] = .5
        #         _A[i2, p] = .5
        #         # z_jk
        #         i1 = (j)*dim + (k) 
        #         i2 = (k)*dim + (j)
        #         _A[i1, p] = -.5
        #         _A[i2, p] = -.5
        #         # z_ik
        #         i1 = (k)*dim + (i) 
        #         i2 = (i)*dim + (k)
        #         _A[i1, p] = .5
        #         _A[i2, p] = .5
        #         # t
        #         i1 = (n)*dim + (n)
        #         _A[i1, p] = -1.
        #         _b[p] = -1.
        #         mleq += 1
        #         l += 1
        #         p += 1
        #         if p == step:
        #             A = hstack([A, _A])
        #             b = np.vstack([b, _b])
        #             _A = lil_matrix((dim**2, step))
        #             _b = np.zeros((step, 1))
        #             p = 0

        # z_ij + z_ik - z_jk - t <= -1, for i, j, k in V [as modeled in DAM (17)]
        try:
            if sample:
                pool = iter_sample_fast(itertools.combinations(G.nodes(), 3), sample)
            else:
                pool = itertools.combinations(G.nodes(), 3)
        except Exception:
            pool = itertools.combinations(G.nodes(), 3)

        for (i, j, k) in pool:
            # z_ij
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
            _A[i1, p] = .5
            _A[i2, p] = .5
            # z_jk
            i1 = (j)*dim + (k) 
            i2 = (k)*dim + (j)
            _A[i1, p] = -.5
            _A[i2, p] = -.5
            # z_ik
            i1 = (k)*dim + (i) 
            i2 = (i)*dim + (k)
            _A[i1, p] = .5
            _A[i2, p] = .5
            # t
            i1 = (n)*dim + (n)
            _A[i1, p] = -1.
            _b[p] = -1.
            mleq += 1
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij >= -1 , for all (i, j) not in E
        L = np.zeros((dim, dim))
        L.fill(-np.inf)
        for i, j in G_complement_edges:
            L[i, j] = -1.
            L[j, i] = -1.

        # for i, j in G_complement_edges:
        #     i1 = (i)*dim + (j) 
        #     i2 = (j)*dim + (i)
        #     _A[i1, p] = -.5
        #     _A[i2, p] = -.5
        #     _b[p] = 1.
        #     mleq += 1
        #     l += 1
        #     p += 1
        #     if p == step:
        #         A = hstack([A, _A])
        #         b = np.vstack([b, _b])
        #         _A = lil_matrix((dim**2, step))
        #         _b = np.zeros((step, 1))
        #         p = 0

        # z_ii - t = -1, for all i in V
        for i in G.nodes():
            i1 = (i)*dim + (i)
            i2 = (n)*dim + (n)
            _A[i1, p] = 1.
            _A[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                A = hstack([A, _A])
                b = np.vstack([b, _b])
                _A = lil_matrix((dim**2, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
            _A[i1, p] = .5
            _A[i2, p] = .5
            _b[p] = -1.
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'L' : L, 'mleq' : mleq})
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
        C = np.zeros((n + 1, n + 1))
        C[n, n] = 1.

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        Bt = lil_matrix((dim, 0))
        u = np.zeros((0, 1))

        _Bt = lil_matrix((dim, step))
        _u = np.zeros((step, 1))
        p1 = 0

        # z_ii - t = -1
        for i in G.nodes():
            i1 = ((i + 1)*(i))/2 + (i)
            i2 = ((n + 1)*(n))/2 + (n)
            _At[i1, p] = 1.
            _At[i2, p] = -1.
            _b[p] = -1.
            l += 1
            p += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij = -1, for all (i, j) in E
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            idx = int(idx)
            _At[idx, p] = _2*0.5
            _b[p] = -1.
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # z_ij >= -1, for all (i, j) not in E
        L = np.zeros((n + 1, n + 1))
        L.fill(-np.inf)
        for i, j in G_complement_edges:
            L[i, j] = -1.
            L[j, i] = -1.

        # for i, j in G_complement_edges:
        #     ii, jj = sorted([i, j], reverse=True)
        #     idx = ((ii + 1)*(ii))/2 + jj
        #     _Bt[idx, p1] = -_2*0.5
        #     _u[p1] = 1.
        #     p1 += 1
        #     l += 1
        #     if p1 == step:
        #         Bt = hstack([Bt, _Bt])
        #         u = np.vstack([u, _u])
        #         _Bt = lil_matrix((dim, step))
        #         _u = np.zeros((step, 1))
        #         p1 = 0

        # z_ij + z_ik - z_jk - t <= -1 , for all (i, j, k) stable
        # for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
        #     if set((i, j)) in G_complement_edges and set((j, k)) in G_complement_edges:
        #         # z_ik
        #         ii, kk = sorted([i, k], reverse=True)
        #         idx = ((ii + 1)*(ii))/2 + kk
        #         _Bt[idx, p1] = _2*.5
 
        #         # z_jk
        #         jj, kk = sorted([j, k], reverse=True)
        #         idx = ((jj + 1)*(jj))/2 + kk
        #         _Bt[idx, p1] = -_2*.5
 
        #         # z_ij
        #         ii, jj = sorted([i, j], reverse=True)
        #         idx = ((ii + 1)*(ii))/2 + jj
        #         _Bt[idx, p1] = _2*.5
 
        #         # t
        #         idx = ((n + 1)*(n))/2 + n
        #         _Bt[idx, p1] = -1.
 
        #         _u[p1] = -1.
 
        #         p1 += 1
        #         l += 1
        #         if p1 == step:
        #             Bt = hstack([Bt, _Bt])
        #             u = np.vstack([u, _u])
        #             _Bt = lil_matrix((dim, step))
        #             _u = np.zeros((step, 1))
        #             p1 = 0

        # z_ij + z_ik - z_jk - t <= -1 , for all (i, j, k) stable
        try:
            if sample:
                pool = iter_sample_fast(itertools.combinations(G.nodes(), 3), sample)
            else:
                pool = itertools.combinations(G.nodes(), 3)
        except Exception:
            pool = itertools.combinations(G.nodes(), 3)

        for (i, j, k) in pool:
            # z_ik
            ii, kk = sorted([i, k], reverse=True)
            idx = ((ii + 1)*(ii))/2 + kk
            idx = int(idx)
            _Bt[idx, p1] = _2*.5

            # z_jk
            jj, kk = sorted([j, k], reverse=True)
            idx = ((jj + 1)*(jj))/2 + kk
            idx = int(idx)
            _Bt[idx, p1] = -_2*.5

            # z_ij
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            idx = int(idx)
            _Bt[idx, p1] = _2*.5

            # t
            idx = ((n + 1)*(n))/2 + n
            idx = int(idx)
            _Bt[idx, p1] = -1.

            _u[p1] = -1.

            p1 += 1
            l += 1
            if p1 == step:
                Bt = hstack([Bt, _Bt])
                u = np.vstack([u, _u])
                _Bt = lil_matrix((dim, step))
                _u = np.zeros((step, 1))
                p1 = 0

        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        if p1 > 0:
            Bt = hstack([Bt,  _Bt.tocsc()[:, :p1]]);
            u = np.vstack([u, _u[:p1]])


        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'L' : L, 'C' : C, 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def M_plus_REP_GCP(G, filename, model_out_dir='', step=10000, debug=False):
    start_time = time.time()
    
    n = len(G.nodes())
    m = len(G.edges())
    
    G_c = nx.complement(G)
    
    # Constraints counter
    _l = 0

    # Sqrt of 2, needed for off-diagonal matrix elements
    _2 = np.sqrt(2)
    
    mapping = {}
    i = 1
    for u in G.nodes():
        for v in itertools.chain([u], G_c.neighbors(u)):
            mapping[u, v] = i
            i += 1
    
    # Number of vars
    _n = i - 1
    
    # Dimension of Y: number of vars (_n) + 1
    dim = int((_n + 1) * (_n + 2) / 2)
    
    # Objective function (minimization problem)
    C = np.zeros((_n + 1, _n + 1))
    for i in range(n):
        C[mapping[i, i], mapping[i, i]] = 1.

    # At will be the output, _At is used just as temporary var
    # each step _At will be concatenated to At
    At = lil_matrix((dim, 0))
    b = np.zeros((0, 1))

    _At = lil_matrix((dim, step))
    _b = np.zeros((step, 1))
    _p = 0

    Bt = lil_matrix((dim, 0))
    u = np.zeros((0, 1))

    _Bt = lil_matrix((dim, step))
    _u = np.zeros((step, 1))
    _p1 = 0
    
    # (1.1) - (1.3)
    for v in G.nodes():
        for i in G.nodes():
            for j in itertools.chain([i], G_c.neighbors(i)):
                # x_ij
                _id1 = mapping[i, j]
                i1 = mat_idx(_id1, _id1)
                _Bt[i1, _p1] -= 1.
                for w in itertools.chain([v], G_c.neighbors(v)):
                    # x_ijx_wv
                    i2 = mat_idx(mapping[w, v], _id1)
                    _Bt[i2, _p1] += _2*.5 if _id1 != mapping[w, v] else 1.
                
                _p1 += 1
                _l += 1
                if _p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    _p1 = 0
                    
    # (1.2) - (1.4)
    for v in G.nodes():
        for i in G.nodes():
            for j in itertools.chain([i], G_c.neighbors(i)):
                # x_ij
                _id1 = mapping[i, j]
                i1 = mat_idx(_id1, _id1)
                _Bt[i1, _p1] += 1.
                for w in itertools.chain([v], G_c.neighbors(v)):
                    # x_wv
                    _id2 = mapping[w, v]
                    i2 = mat_idx(_id2, _id2)
                    _Bt[i2, _p1] += 1.
                    # x_ijx_uv
                    i3 = mat_idx(_id1, _id2)
                    _Bt[i3, _p1] -= _2*.5 if _id1 != _id2 else 1.
                
                _u[_p1] = 1.
                
                _p1 += 1
                _l += 1
                if _p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    _p1 = 0
                    
    # (2.1)
    for uu in G.nodes():
        for v,w in G.subgraph(itertools.chain(G_c.neighbors(uu), [uu])).edges():
            for i in G.nodes():
                for j in itertools.chain([i], G_c.neighbors(i)):
                    # x_ijx_uu
                    _id1 = mapping[i, j]
                    _id2 = mapping[uu, uu]
                    i1 = mat_idx(_id1, _id2)
                    _Bt[i1, _p1] += _2*.5 if _id1 != _id2 else 1.
                    
                    # x_ijx_uw
                    _id2 = mapping[uu, w]
                    i2 = mat_idx(_id1, _id2)
                    _Bt[i2, _p1] -= _2*.5 if _id1 != _id2 else 1.
                    
                    # x_ijx_uv
                    _id2 = mapping[uu, v]
                    i3 = mat_idx(_id1, _id2)
                    _Bt[i3, _p1] -= _2*.5 if _id1 != _id2 else 1.
                    
                    _p1 += 1
                    _l += 1
                    if _p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        _p1 = 0
    
    # (2.2)
    #for uu in G.nodes():
    #    for v,w in G.subgraph(itertools.chain(G_c.neighbors(uu), [uu])).edges():
    #        for i in G.nodes():
    #            for j in itertools.chain([i], G_c.neighbors(i)):
    #                # x_ijx_uu
    #                _id1 = mapping[i, j]
    #                _id2 = mapping[uu, uu]
    #                i1 = mat_idx(_id1, _id2)
    #                _Bt[i1, _p1] -= _2*.5 if _id1 != _id2 else 1.
    #                
    #                # x_ijx_uw
    #                _id2 = mapping[uu, w]
    #                i2 = mat_idx(_id1, _id2)
    #                _Bt[i2, _p1] += _2*.5 if _id1 != _id2 else 1.
    #                
    #                # x_ijx_uv
    #                _id2 = mapping[uu, v]
    #                i3 = mat_idx(_id1, _id2)
    #                _Bt[i3, _p1] += _2*.5 if _id1 != _id2 else 1.
    #                
    #                # x_uu
    #                _id1 = mapping[uu, uu]
    #                i4 = mat_idx(_id1, _id1)
    #                _Bt[i4, _p1] += 1.
    #                
    #                # x_uw
    #                _id1 = mapping[uu, w]
    #                i5 = mat_idx(_id1, _id1)
    #                _Bt[i5, _p1] -= 1.
    #                
    #                # x_uv
    #                _id1 = mapping[uu, v]
    #                i6 = mat_idx(_id1, _id1)
    #                _Bt[i6, _p1] -= 1.
    #                
    #                _p1 += 1
    #                _l += 1
    #                if _p1 == step:
    #                    Bt = hstack([Bt, _Bt])
    #                    u = np.vstack([u, _u])
    #                    _Bt = lil_matrix((dim, step))
    #                    _u = np.zeros((step, 1))
    #                    _p1 = 0
                    
    # Add non-negativity on Y
    L = 0.

    # Add constraint: Y_00 = 1
    _At[0, _p] = 1.
    _b[_p] = 1.
    _p += 1
    _l += 1
    if _p == step:
        At = hstack([At, _At])
        b = np.vstack([b, _b])
        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        _p = 0

    # Add constraint: Y_ii = Y_0i = Y_i0 
    for i in range(_n):
        idx = mat_idx(i + 1, 0)
        
        _At[idx, _p] = 1
        idx = mat_idx(i + 1, i + 1)
        _At[idx, _p] = -1
        _p += 1
        _l += 1
        if _p == step:
            At = hstack([At, _At])
            b = np.vstack([b, _b])
            _At = lil_matrix((dim, step))
            _b = np.zeros((step, 1))
            _p = 0
                    
    if _p > 0:
        At = hstack([At,  _At.tocsc()[:, :_p]])
        b = np.vstack([b, _b[:_p]])
    if _p1 > 0:
        Bt = hstack([Bt,  _Bt.tocsc()[:, :_p1]]);
        u = np.vstack([u, _u[:_p1]])

    elapsed_time = time.time() - start_time
    
    size_f = (At.shape[1], Bt.shape[1])
    
    At, id_At = sp_unique(At, axis=1)
    Bt, id_Bt = sp_unique(Bt, axis=1)
    
    size_u = (At.shape[1], Bt.shape[1])
    
    At.eliminate_zeros()
    Bt.eliminate_zeros()

    if debug:
        print('Finished! Time elapsed: %.2f' % elapsed_time)
        print('Dimension of matrix variable: %d' % (_n + 1))
        print('Number of constraints: %d' % (size_f[0] + size_f[1]))
        print('Number of UNIQUE constraints: %d' % (size_u[0] + size_u[1]))
    savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b[id_At], 'Bt' : Bt, 'u' : u[id_Bt], 'C' : C, 'L' : L, 's' : float(_n + 1)})
    # Exporting the model to output file
    if debug:
        print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def M_plus_GCP(G, filename, ub=0, cliques=False, model_out_dir='', step=10000, debug=False):
    start_time = time.time()
    n = len(G.nodes())
    m = len(G.edges())

    if ub > 0:
        _n = ub + n*ub
        n = ub
    else:
        _n = n*n + n

    # Constraints counter
    _l = 0

    # Sqrt of 2, needed for off-diagonal matrix elements
    _2 = np.sqrt(2) 

    # Dimension of Y: number of vars (_n) + 1
    dim = int((_n + 1) * (_n + 2) / 2)
    
    idxs = list(set(list(range(dim))) - set([mat_idx(i, i) for i in range(_n)]))

    # Objective function (minimization problem)
    C = np.zeros((_n + 1, _n + 1))
    for i in range(n):
        C[i + 1, i + 1] = 1.

    # At will be the output, _At is used just as temporary var
    # each step _At will be concatenated to At
    At = lil_matrix((dim, 0))
    b = np.zeros((0, 1))

    _At = lil_matrix((dim, step))
    _b = np.zeros((step, 1))
    _p = 0

    Bt = lil_matrix((dim, 0))
    u = np.zeros((0, 1))

    _Bt = lil_matrix((dim, step))
    _u = np.zeros((step, 1))
    _p1 = 0

    # Indexes of variables that are fixed to 0 by (1.5), will be removed from constraints at the end of the formulation
    fix_var = []

    for i in G.nodes():
        # 1.6
        for j in range(n):
            _id1 = n + n*i + j + 1
            # x_ij
            i1 = mat_idx(_id1, _id1)
            _At[i1, _p] += 1.
 
        # Add 1.6
        _b[_p] = 1.
        _l += 1
        _p += 1
        if _p == step:
            At = hstack([At, _At])
            b = np.vstack([b, _b])
            _At = lil_matrix((dim, step))
            _b = np.zeros((step, 1))
            _p = 0

        for k in range(n):
            # 1.1 - 1.2 - 1.5
            _a = lil_matrix((dim, 3))
            _id1 = k + 1
            # w_k
            i1 = mat_idx(_id1, _id1)
            _a[i1, 0] += -1.
            _a[i1, 1] += 1.
            for j in range(n):
                _id2 = n + n*i + j + 1
                # w_k-x_ij
                i2 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                _a[i2, 0] += _2*.5
                _a[i2, 1] += -_2*.5
                # x_ij
                i3 = mat_idx(_id2, _id2)
                _a[i3, 1] += 1.
                if j != k:
                    _id3 = n + n*i + k + 1
                    # x_ik-x_ij
                    i4 = mat_idx(_id2, _id3) if _id2 >= _id3 else mat_idx(_id3, _id2)
                    _a[i4, 2] += _2*.5
                    fix_var.append(i4)

            # Add 1.1
            _At[:, _p] = _a[:, 0]
            _l += 1
            _p += 1
            if _p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                _p = 0

            # Add 1.2
            _At[:, _p] = _a[:, 1]
            _b[_p] = 1.
            _l += 1
            _p += 1
            if _p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                _p = 0

            # Add 1.5
            _At[:, _p] = _a[:, 2]
            _l += 1
            _p += 1
            if _p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                _p = 0

            for l in G.nodes():
                if l != i:
                    # 1.3 - 1.4
                    _a = lil_matrix((dim, 2))
                    _id1 = n + n*l + k + 1
                    # x_lk
                    i1 = mat_idx(_id1, _id1)
                    _a[i1, 0] += -1.
                    _a[i1, 1] += 1.
                    for j in range(n):
                        _id2 = n + n*i + j + 1
                        # x_ij
                        i2 = mat_idx(_id2, _id2)
                        _a[i2, 1] += 1.
                        # x_lk-x_ij
                        i3 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                        _a[i3, 0] += _2*.5
                        _a[i3, 1] += -_2*.5

                    # Add 1.3
                    _At[:, _p] = _a[:, 0]
                    _l += 1
                    _p += 1
                    if _p == step:
                        At = hstack([At, _At])
                        b = np.vstack([b, _b])
                        _At = lil_matrix((dim, step))
                        _b = np.zeros((step, 1))
                        _p = 0

                    # Add 1.4
                    _At[:, _p] = _a[:, 1]
                    _b[_p] = 1.
                    _l += 1
                    _p += 1
                    if _p == step:
                        At = hstack([At, _At])
                        b = np.vstack([b, _b])
                        _At = lil_matrix((dim, step))
                        _b = np.zeros((step, 1))
                        _p = 0

    if debug:
        print('------------------------')
        print('End Equality Constraints:', _l)

    for i, k in G.edges():
        for j in range(n):
            for l in G.nodes():
                for p in range(n):
                    # 2.1 - 2.2
                    _a = lil_matrix((dim, 2))
                    _id1, _id2 = n + n*l + p + 1, n + n*i + j + 1
                    # x_lp-x_ij
                    i1 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                    _a[i1, 0] += 1.  if _id1 == _id2 else _2*.5
                    _a[i1, 1] += -1. if _id1 == _id2 else -_2*.5
                    _id1, _id2 = n + n*l + p + 1, n + n*k + j + 1
                    # x_lp-x_kj
                    i2 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                    _a[i2, 0] += 1.  if _id1 == _id2 else _2*.5
                    _a[i2, 1] += -1. if _id1 == _id2 else -_2*.5
                    _id1, _id2 = n + n*l + p + 1, j + 1
                    # w_j-x_lp
                    i3 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                    _a[i3, 0] += -_2*.5
                    _a[i3, 1] += _2*.5
                    _id1 = n + i*n + j + 1
                    # x_ij
                    i4 = mat_idx(_id1, _id1)
                    _a[i4, 1] += 1.
                    _id1 = n + k*n + j + 1
                    # x_kj
                    i5 = mat_idx(_id1, _id1)
                    _a[i5, 1] += 1.
                    _id1 = j + 1
                    # w_j
                    i6 = mat_idx(_id1, _id1)
                    _a[i6, 1] += -1.

                    # Add 2.1
                    _Bt[:, _p1] = _a[:, 0]
                    _p1 += 1
                    _l += 1
                    if _p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        _p1 = 0

                    # Add 2.2
                    _Bt[:, _p1] = _a[:, 1]
                    _p1 += 1
                    _l += 1
                    if _p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        _p1 = 0

            for p in range(n):
                # 2.3 - 2.4
                _a = lil_matrix((dim, 2))
                _id1, _id2 = n + i*n + j + 1, p + 1
                # w_p-x_ij
                i1 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                _a[i1, 0] += 1.
                _a[i1, 1] += -1.
                _id1, _id2 = n + k*n + j + 1, p + 1
                # w_p-x_kj
                i2 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                _a[i2, 0] += 1.
                _a[i2, 1] += -1.
                _id1, _id2 = j + 1, p + 1
                # w_p-w_j
                i3 = mat_idx(_id1, _id2) if _id1 >= _id2 else mat_idx(_id2, _id1)
                _a[i3, 0] += -1. if _id1 == _id2 else -_2*.5
                _a[i3, 1] += 1.  if _id1 == _id2 else _2*.5
                _id1 = n + n*i + j + 1
                # x_ij
                i1 = mat_idx(_id1, _id1)
                _a[i1, 1] += 1.
                _id2 = n + n*k + j + 1
                # x_kj
                i2 = mat_idx(_id2, _id2)
                _a[i2, 1] += 1.
                _id3 = j + 1
                # w_j
                i3 = mat_idx(_id3, _id3)
                _a[i3, 1] += -1.

                # Add 2.3
                _Bt[:, _p1] = _a[:, 0]
                _p1 += 1
                _l += 1
                if _p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    _p1 = 0

                # Add 2.4
                _Bt[:, _p1] = _a[:, 1]
                _p1 += 1
                _l += 1
                if _p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    _p1 = 0

    if debug:
        print('------------------------')
        print('End Inequality Constraints:', _l)

    # Add non-negativity on Y
    L = 0.

    # Add constraint: Y_00 = 1
    _At[0, _p] = 1.
    _b[_p] = 1.
    _p += 1
    _l += 1
    if _p == step:
        At = hstack([At, _At])
        b = np.vstack([b, _b])
        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        _p = 0

    # Add constraint: Y_ii = Y_0i = Y_i0 
    for i in range(_n):
        idx = mat_idx(i + 1, 0)
        _At[idx, _p] = 1
        idx = mat_idx(i + 1, i + 1)
        _At[idx, _p] = -1
        _p += 1
        _l += 1
        if _p == step:
            At = hstack([At, _At])
            b = np.vstack([b, _b])
            _At = lil_matrix((dim, step))
            _b = np.zeros((step, 1))
            _p = 0

    clique_count = 0

    if cliques:
        for c in nx.find_cliques(G):
            # Add the clique inequality for each color
            for j in range(n):
                # w_j
                i1 = mat_idx(j + 1, j + 1)
                _Bt[i1, _p1] = -1.
                for i in c:
                    # x_ij
                    i2 = mat_idx(n + i*n + j + 1, n + i*n + j + 1)
                    _Bt[i2, _p1] = 1.

                _p1 += 1
                _l += 1
                clique_count += 1
                if _p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    _p1 = 0

                    
    if _p > 0:
        At = hstack([At,  _At.tocsc()[:, :_p]])
        b = np.vstack([b, _b[:_p]])
    if _p1 > 0:
        Bt = hstack([Bt,  _Bt.tocsc()[:, :_p1]]);
        u = np.vstack([u, _u[:_p1]])

    elapsed_time = time.time() - start_time

    start_time = time.time()

    size_f = (At.shape[1], Bt.shape[1])
    
    print(size_f[0] + size_f[1])

    At, id_At = sp_unique(At, axis=1)
    Bt, id_Bt = sp_unique(Bt, axis=1)

    size_u = (At.shape[1], Bt.shape[1])

    At = At.tocsr()
    Bt = Bt.tocsr()

    At[fix_var, :] = 0.
    Bt[fix_var, :] = 0.
    
    At.eliminate_zeros()
    Bt.eliminate_zeros()

    elapsed_time1 = time.time() - start_time

    if debug:
        print('Finished! Time elapsed: %.2f' % elapsed_time)
        print('Variable fixing Time elapsed: %.2f' % elapsed_time1)
        print('Variable fixed: %d' % len(fix_var))
        print('Dimension of matrix variable: %d' % (_n + 1))
        print('Number of constraints: %d' % (size_f[0] + size_f[1]))
        print('Number of UNIQUE constraints: %d' % (size_u[0] + size_u[1]))
        print('Number of CLIQUES inequalities: %d' % clique_count)
    savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b[id_At], 'Bt' : Bt, 'u' : u[id_Bt], 'C' : C, 'L' : L, 's' : float(_n + 1)})
    # Exporting the model to output file
    if debug:
        print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')

