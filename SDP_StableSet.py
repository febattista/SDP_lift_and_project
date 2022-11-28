#coding:utf-8
from graph import *
import time
import pandas as pd
import gurobipy as gb
import numpy as np
import networkx as nx
import random
import math
from scipy.sparse import lil_matrix, vstack, hstack, save_npz, csc_matrix, lil_matrix
from scipy.io import savemat
import sys
import itertools
import subprocess
import os

intpnt_co_tolerance = 1.0e-7

def read_model_from_file(path):
    try:
        model = gb.read(path)
    except Exception as error:
        print(error)
    return model

def FRAC(G, debug=False):
    stab_G = gb.Model()
    x = stab_G.addVars(G.nodes(), vtype=gb.GRB.BINARY, name='x')
    stab_G.setObjective(x.sum(), gb.GRB.MAXIMIZE)
    stab_G.addConstrs(x[u] + x[v] <= 1 for u, v in G.edges())
    stab_G.update()
    model = stab_G

    return model


def NSTAB(G, model_name, dir_path=''):
    path = os.path.join(dir_path, model_name + '_nod.lp')
    nod = gb.Model()
    x = nod.addVars(G.nodes(), vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    nod.setObjective(x.sum(), gb.GRB.MAXIMIZE)

    for i in G.nodes():
        neighbors = list(G.neighbors(i))
        if neighbors:
            # Compute alpha of G[N(i)] by computing the max clique on the complement
            alpha = nx.graph_clique_number(nx.complement(G.subgraph(neighbors)))
            nod.addConstr(x.sum(neighbors) + alpha*x[i] <= alpha, name='nod' + str(i))

    nod.update()
    nod.write(path)
    return nod

def Theta_SDP1(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_SDP1')
        print('Reading the model...')

    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation 
        n = len(G.nodes())
        dim = n 
        m = len(G.edges())
        
        l = 0
        mleq = 0

        # Check in advance if the problem size will exceed the limit
        if limit and l > limit:
            if debug:
                print("Problem is too big: %d rows" % l)
            return

        # Objective function (minimization problem)
        J = -np.ones((dim, dim))

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        # x_ij = 0 , for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
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

        # tr(X) = 1
        for i in G.nodes():
        	i1 = (i)*dim + (i)
        	_A[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1

        if p > 0:
            A = hstack([A,  _A.tocsc()[:, :p]]);
            b = np.vstack([b, _b[:p]])
        A = A.transpose()

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % dim)
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : J, 'mleq' : mleq})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


    elif model_out.lower() == 'sdpnal': # SDPNAL model creation

        n = len(G.nodes())
        m = len(G.edges())

        dim = int(((n) * (n + 1))/2)
        
        l = 0

        _2 = np.sqrt(2)

        # Objective function (minimization problem)
        J = -np.ones((n, n))

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
            idx = ((ii + 1)*(ii))/2 + jj
            _At[idx, p] = _2*0.5
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # tr(X) = 1
        for i in G.nodes():
        	i1 = ((i + 1)*(i))/2 + (i)
        	_At[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1
        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : J, 's' : float(n)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_plus_SDP1(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_SDP1')
        print('Reading the model...')

    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation 
        n = len(G.nodes())
        dim = n 
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
        J = -np.ones((dim, dim))

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0

        # Z_ij > 0 , for all (i, j) not in E
        for i, j in G_complement_edges:
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
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

        # Z_ij = 0 , for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
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

        # tr(Z) = 1
        for i in G.nodes():
        	i1 = (i)*dim + (i)
        	_A[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1

        if p > 0:
            A = hstack([A,  _A.tocsc()[:, :p]]);
            b = np.vstack([b, _b[:p]])
        A = A.transpose()



        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % dim)
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : J, 'mleq' : mleq, 'L': np.zeros((n, n))})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


    elif model_out.lower() == 'sdpnal': # SDPNAL model creation

        n = len(G.nodes())
        m = len(G.edges())

        dim = int(((n) * (n + 1))/2)
        
        l = 0

        _2 = np.sqrt(2)

        # Objective function (minimization problem)
        J = -np.ones((n, n))

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        At = lil_matrix((dim, 0))
        b = np.zeros((0, 1))

        _At = lil_matrix((dim, step))
        _b = np.zeros((step, 1))
        p = 0

        # Z_ij == 0
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            _At[idx, p] = _2*0.5
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # tr(Z) = 1
        for i in G.nodes():
        	i1 = ((i + 1)*(i))/2 + (i)
        	_At[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1
        
        if p > 0:
            At = hstack([At,  _At.tocsc()[:, :p]])
            b = np.vstack([b, _b[:p]])

        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % (n + 1))
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : J, 'L' : 0., 's' : float(n)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def DukanovicRendl(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000, dnn=False):
    # Dukanovic and Rendl
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('DukanovicRendl')
        print('Reading the model...')

    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation 
        n = len(G.nodes())
        dim = n 
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
        J = -np.ones((dim, dim))

        # A will be the output, _A is used just as temporary var
        # each step _A will be concatenated to A
        A = lil_matrix((dim**2, 0))
        b = np.zeros((0, 1))

        _A = lil_matrix((dim**2, step))
        _b = np.zeros((step, 1))
        p = 0
        
        # Z_ik + Z_jk <= Z_ij + Z_kk, for i, j, k stable
        for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
            if len(G.subgraph((i, j, k)).edges()) == 0:
            	# Z_ik
                i1 = (i)*dim + (k) 
                i2 = (k)*dim + (i)
                _A[i1, p] = .5
                _A[i2, p] = .5
                # Z_jk
                i1 = (j)*dim + (k) 
                i2 = (k)*dim + (j)
                _A[i1, p] = .5
                _A[i2, p] = .5
                # Z_ij
                i1 = (j)*dim + (i) 
                i2 = (i)*dim + (j)
                _A[i1, p] = -.5
                _A[i2, p] = -.5
                # Z_kk
                i1 = (k)*dim + (k)
                _A[i1, p] = -1.
                mleq += 1
                l += 1
                p += 1
                if p == step:
                    A = hstack([A, _A])
                    b = np.vstack([b, _b])
                    _A = lil_matrix((dim**2, step))
                    _b = np.zeros((step, 1))
                    p = 0

        # Z_ik + Z_jk <=  Z_kk, for i,j in E, k \= i, j
        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # x_kk
                    i1 = (k)*dim + (k) 
                    _A[i1, p] = -1.
                    # x_jk
                    i1 = (j)*dim + (k)
                    i2 = (k)*dim + (j)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    # x_ik
                    i1 = (i)*dim + (k)
                    i2 = (k)*dim + (i)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    l += 1
                    mleq += 1
                    p += 1
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0

        # Z_ij > 0 , for all (i, j) not in E
        if not dnn:
            for i, j in G_complement_edges:
                i1 = (i)*dim + (j) 
                i2 = (j)*dim + (i)
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

        # Z_ij = 0 , for all (i, j) in E
        for i, j in G.edges():
            i1 = (i)*dim + (j) 
            i2 = (j)*dim + (i)
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

        # tr(Z) = 1
        for i in G.nodes():
        	i1 = (i)*dim + (i)
        	_A[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1

        if p > 0:
            A = hstack([A,  _A.tocsc()[:, :p]]);
            b = np.vstack([b, _b[:p]])
        A = A.transpose()



        elapsed_time = time.time() - start_time
        if debug:
            print('Finished! Time elapsed: %.2f' % elapsed_time)
            print('Dimension of matrix variable: %d' % dim)
            print('Number of constraints: %d' % l)
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : J, 'mleq' : mleq, 'L': np.zeros((n, n))})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


    elif model_out.lower() == 'sdpnal': # SDPNAL model creation

        n = len(G.nodes())
        m = len(G.edges())

        dim = int(((n) * (n + 1))/2)
        
        l = 0

        _2 = np.sqrt(2)

        # Objective function (minimization problem)
        J = -np.ones((n, n))

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

        # Z_ik + Z_jk < Z_ij + Z_kk , for all (i, j, k) stable
        for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
            if len(G.subgraph((i, j, k)).edges()) == 0:
                # Z_kk
                idx = ((k + 1)*(k))/2 + k
                _Bt[idx, p1] = -1.

                # Z_ik
                ii, kk = sorted([i, k], reverse=True)
                idx = ((ii + 1)*(ii))/2 + kk
                _Bt[idx, p1] = _2*.5

                # Z_jk
                jj, kk = sorted([j, k], reverse=True)
                idx = ((jj + 1)*(jj))/2 + kk
                _Bt[idx, p1] = _2*.5

                # Z_ij
                ii, jj = sorted([i, j], reverse=True)
                idx = ((ii + 1)*(ii))/2 + jj
                _Bt[idx, p1] = -_2*.5

                p1 += 1
                l += 1
                if p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    p1 = 0

        # Z_ik + Z_jk < Z_kk, for all (i, j) \in E, k \= i, j
        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # Z_kk
                    idx = ((k + 1)*(k))/2 + k
                    _Bt[idx, p1] = -1.

                    # Z_ik
                    ii, kk = sorted([i, k], reverse=True)
                    idx = ((ii + 1)*(ii))/2 + kk
                    _Bt[idx, p1] = _2*.5

                    # Z_jk
                    jj, kk = sorted([j, k], reverse=True)
                    idx = ((jj + 1)*(jj))/2 + kk
                    _Bt[idx, p1] = _2*.5

                    p1 += 1
                    l += 1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

        # Z_ij == 0
        for i, j in G.edges():
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 1)*(ii))/2 + jj
            _At[idx, p] = _2*0.5
            p += 1
            l += 1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0

        # tr(Z) = 1
        for i in G.nodes():
        	i1 = ((i + 1)*(i))/2 + (i)
        	_At[i1, p] = 1.

        _b[p] = 1.
        l += 1
        p += 1
        
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat',  {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'C' : J, 'L' : 0., 's' : float(n)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_SDP2(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('Theta_SDP2')
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : C, 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def Theta_plus_SDP2(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq, 'L': np.zeros((dim, dim))})
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'C' : C, 'L' : 0., 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def GruberRendl(G, filename, limit=None, export_matrices=False, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    # Gruber and Rendl
    if debug:
        print('GruberRendl')
        print('Reading the model...')

    start_time = time.time()
    n = len(G.nodes())
    dim = n + 1
    m = len(G.edges())
    G_complement_edges = [i for i in itertools.combinations(G.nodes(), 2) if i not in G.edges()]
    m_complement = len(G_complement_edges)

    if model_out.lower() == 'adal':
        l = 0
        mleq = 0

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

        # x_ik + x_jk <  x_kk , for all (i, j) in E, k \not= i, j
        # x_ii + x_jj + x_kk < 1 + x_ik + x_jk , for all (i, j) in E, k \not= i, j
        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1) 
                    _A[i1, p] = -1.
                    # x_jk
                    i1 = (j + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (j + 1)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    # x_ik
                    i1 = (i + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (i + 1)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    l += 1
                    mleq += 1
                    p += 1
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0
                    
                    # x_ii
                    i1 = (i + 1)*dim + (i + 1) 
                    _A[i1, p] = 1.
                    # x_jj
                    i1 = (j + 1)*dim + (j + 1) 
                    _A[i1, p] = 1.
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1) 
                    _A[i1, p] = 1.
                    # x_jk
                    i1 = (j + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (j + 1)
                    _A[i1, p] = -.5
                    _A[i2, p] = -.5
                    # x_ik
                    i1 = (i + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (i + 1)
                    _A[i1, p] = -.5
                    _A[i2, p] = -.5
                    _b[p] = 1.
                    l += 1
                    mleq += 1
                    p += 1
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0
        
        # x_ik + x_jk < x_ij + x_kk , for all (i, j, k) stable
        # x_ii + x_jj + x_kk < 1 + x_ij + x_ik + x_jk , for all (i, j, k) stable
        for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
            if len(G.subgraph((i, j, k)).edges()) == 0:
                # (5)
                # x_ik
                i1 = (i + 1)*dim + (k + 1) 
                i2 = (k + 1)*dim + (i + 1)
                _A[i1, p] = .5
                _A[i2, p] = .5
                # x_jk
                i1 = (j + 1)*dim + (k + 1) 
                i2 = (k + 1)*dim + (j + 1)
                _A[i1, p] = .5
                _A[i2, p] = .5
                # x_ij
                i1 = (j + 1)*dim + (i + 1) 
                i2 = (i + 1)*dim + (j + 1)
                _A[i1, p] = -.5
                _A[i2, p] = -.5
                # x_kk
                i1 = (k + 1)*dim + (k + 1)
                _A[i1, p] = -1.
                mleq += 1
                l += 1
                p += 1
                if p == step:
                    A = hstack([A, _A])
                    b = np.vstack([b, _b])
                    _A = lil_matrix((dim**2, step))
                    _b = np.zeros((step, 1))
                    p = 0

                # x_ii
                i1 = (i + 1)*dim + (i + 1)
                _A[i1, p] = 1.
                # x_jj
                i1 = (j + 1)*dim + (j + 1)
                _A[i1, p] = 1.
                # x_kk
                i1 = (k + 1)*dim + (k + 1)
                _A[i1, p] = 1.
                # x_ij
                i1 = (j + 1)*dim + (i + 1) 
                i2 = (i + 1)*dim + (j + 1)
                _A[i1, p] = -.5
                _A[i2, p] = -.5
                # x_ik
                i1 = (k + 1)*dim + (i + 1) 
                i2 = (i + 1)*dim + (k + 1)
                _A[i1, p] = -.5
                _A[i2, p] = -.5
                # x_jk
                i1 = (j + 1)*dim + (k + 1) 
                i2 = (k + 1)*dim + (j + 1)
                _A[i1, p] = -.5
                _A[i2, p] = -.5
                _b[p] = 1.
                mleq += 1
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq, 'L': np.zeros((dim, dim))})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')

    elif model_out.lower() == 'sdpnal':
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

        Bt = lil_matrix((dim, 0))
        u = np.zeros((0, 1))

        _Bt = lil_matrix((dim, step))
        _u = np.zeros((step, 1))
        p1 = 0

        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] = -1.

                    # x_ik
                    ii, kk = sorted([i, k], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                    _Bt[idx, p1] = _2*.5

                    # x_jk
                    jj, kk = sorted([j, k], reverse=True)
                    idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                    _Bt[idx, p1] = _2*.5

                    p1 += 1
                    l += 1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

                    # x_ii
                    idx = ((i + 2)*(i + 1))/2 + i + 1
                    _Bt[idx, p1] = 1.

                    # x_jj
                    idx = ((j + 2)*(j + 1))/2 + j + 1
                    _Bt[idx, p1] = 1.

                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] = 1.

                    # x_ik
                    ii, kk = sorted([i, k], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                    _Bt[idx, p1] = -_2*.5

                    # x_jk
                    jj, kk = sorted([j, k], reverse=True)
                    idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                    _Bt[idx, p1] = -_2*.5

                    _u[p1] = 1.

                    p1 += 1
                    l += 1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

        # x_ik + x_jk < x_ij + x_kk , for all (i, j, k) stable
        # x_ii + x_jj + x_kk < 1 + x_ij + x_ik + x_jk , for all (i, j, k) stable
        for (i, j, k) in [i for i in  itertools.combinations(G.nodes(), 3)]:
            if len(G.subgraph((i, j, k)).edges()) == 0:
                # x_kk
                idx = ((k + 2)*(k + 1))/2 + k + 1
                _Bt[idx, p1] = -1.

                # x_ik
                ii, kk = sorted([i, k], reverse=True)
                idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                _Bt[idx, p1] = _2*.5

                # x_jk
                jj, kk = sorted([j, k], reverse=True)
                idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                _Bt[idx, p1] = _2*.5

                # x_ij
                ii, jj = sorted([i, j], reverse=True)
                idx = ((ii + 2)*(ii + 1))/2 + jj + 1
                _Bt[idx, p1] = -_2*.5

                p1 += 1
                l += 1
                if p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    p1 = 0

                # x_ii
                idx = ((i + 2)*(i + 1))/2 + i + 1
                _Bt[idx, p1] = 1.

                # x_jj
                idx = ((j + 2)*(j + 1))/2 + j + 1
                _Bt[idx, p1] = 1.

                # x_kk
                idx = ((k + 2)*(k + 1))/2 + k + 1
                _Bt[idx, p1] = 1.

                # x_ik
                ii, kk = sorted([i, k], reverse=True)
                idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                _Bt[idx, p1] = -_2*.5

                # x_jk
                jj, kk = sorted([j, k], reverse=True)
                idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                _Bt[idx, p1] = -_2*.5

                # x_ij
                ii, jj = sorted([i, j], reverse=True)
                idx = ((ii + 2)*(ii + 1))/2 + jj + 1
                _Bt[idx, p1] = -_2*.5

                _u[p1] = 1.

                p1 += 1
                l += 1
                if p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    p1 = 0

        # Add constraint: X_00 = 1
        _At[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        # Add constraint: x_ij = 0 for i, j in E
        for i, j in G.edges():
            # x_ij
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 2)*(ii + 1))/2 + jj + 1
            _At[idx, p] = _2*.5
            
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
            l+=1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0
        
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'C' : C, 'L' : 0., 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def LovaszSchrijver(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', step=10000):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('LovaszSchrijver')
        print('Reading the model...')
    
    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation
        n = len(G.nodes())
        dim = n + 1
        m = len(G.edges())
        G_complement_edges = [i for i in itertools.combinations(G.nodes(), 2) if i not in G.edges()]
        m_complement = len(G_complement_edges)
        
        l = 0
        mleq = 0

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

        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1) 
                    _A[i1, p] = -1.
                    # x_jk
                    i1 = (j + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (j + 1)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    # x_ik
                    i1 = (i + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (i + 1)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                    l += 1
                    mleq += 1
                    p += 1
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0
                               
                    # x_ii
                    i1 = (i + 1)*dim + (i + 1) 
                    _A[i1, p] = 1.
                    # x_jj
                    i1 = (j + 1)*dim + (j + 1) 
                    _A[i1, p] = 1.
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1) 
                    _A[i1, p] = 1.
                    # x_jk
                    i1 = (j + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (j + 1)
                    _A[i1, p] = -.5
                    _A[i2, p] = -.5
                    # x_ik
                    i1 = (i + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (i + 1)
                    _A[i1, p] = -.5
                    _A[i2, p] = -.5
                    _b[p] = 1.
                    l += 1
                    mleq += 1
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'A' : A, 'b' : b, 'C' : C, 'mleq' : mleq,  'L': np.zeros((dim, dim))})
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

        Bt = lil_matrix((dim, 0))
        u = np.zeros((0, 1))

        _Bt = lil_matrix((dim, step))
        _u = np.zeros((step, 1))
        p1 = 0

        for i, j in G.edges():
            for k in G.nodes():
                if k != i and k != j:
                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] = -1.

                    # x_ik
                    ii, kk = sorted([i, k], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                    _Bt[idx, p1] = _2*.5

                    # x_jk
                    jj, kk = sorted([j, k], reverse=True)
                    idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                    _Bt[idx, p1] = _2*.5

                    p1 += 1
                    l += 1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

                    # x_ii
                    idx = ((i + 2)*(i + 1))/2 + i + 1
                    _Bt[idx, p1] = 1.

                    # x_jj
                    idx = ((j + 2)*(j + 1))/2 + j + 1
                    _Bt[idx, p1] = 1.

                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] = 1.

                    # x_ik
                    ii, kk = sorted([i, k], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                    _Bt[idx, p1] = -_2*.5

                    # x_jk
                    jj, kk = sorted([j, k], reverse=True)
                    idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                    _Bt[idx, p1] = -_2*.5

                    _u[p1] = 1.

                    p1 += 1
                    l += 1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

        # Add constraint: X_00 = 1
        _At[0, p] = 1.
        _b[p] = 1.
        p += 1
        l += 1

        # Add constraint: x_ij = 0 for i, j in E
        for i, j in G.edges():
            # x_ij
            ii, jj = sorted([i, j], reverse=True)
            idx = ((ii + 2)*(ii + 1))/2 + jj + 1
            _At[idx, p] = _2*.5
            
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
            l+=1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0
        
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'C' : C, 'L' : 0., 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')


def M_plus_NOD(G, filename, limit=None, debug=False, model_out_dir='', model_out='adal', to_delete=[], step=10000, alpha=None):
    assert type(model_out) == type('') and model_out.lower() in {'adal', 'sdpnal'}, 'Supported models : adal, sdpnal'
    if debug:
        print('M_plus_NOD')
        print('Reading the model...')
    start_time = time.time()

    if model_out.lower() == 'adal': # ADAL model creation
        dim = len(G.nodes()) + 1
        n = len(G.nodes())
        m = len(G.edges())
        
        l = 0
        mleq = 0

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

        # If not given, compute alphas for each node
        if alpha is None:
            alpha = {}
            for i in G.nodes():
                neighbors = list(G.neighbors(i))
                if neighbors:
                    # Compute alpha of G[N(i)] by computing the max clique on the complement
                    alpha[i] = nx.graph_clique_number(nx.complement(G.subgraph(neighbors)))
        
        # (1.1)
        if 1.1 not in to_delete:
            for i in G.nodes():
                for j in G.neighbors(i):
                    # x_ij
                    i1 = (i + 1)*dim + (j + 1) 
                    i2 = (j + 1)*dim + (i + 1)
                    _A[i1, p] = .5
                    _A[i2, p] = .5
                l += 1
                mleq += 1
                p += 1
                if p == step:
                    A = hstack([A, _A])
                    b = np.vstack([b, _b])
                    _A = lil_matrix((dim**2, step))
                    _b = np.zeros((step, 1))
                    p = 0

        # (2.1)
        if 2.1 not in to_delete:
            for i in G.nodes():
                for k in G.neighbors(i):
                    neighbors_but_k = [n for n in G.neighbors(i) if n != k]
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1)
                    _A[i1, p] = float(1 - alpha[i])
                    for j in neighbors_but_k:
                        if (j, k) not in G.edges():
                            # x_jk
                            i1 = (j + 1)*dim + (k + 1)
                            i2 = (k + 1)*dim + (j + 1)
                            _A[i1, p] = .5
                            _A[i2, p] = .5
                    l += 1
                    mleq += 1
                    p += 1       
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0                 

        # (3.1)
        if 3.1 not in to_delete:           
            for i in G.nodes():
                not_neighbors = [n for n in G.nodes() if n != i and n not in G.neighbors(i)]
                for k in not_neighbors:
                    # x_kk
                    i1 = (k + 1)*dim + (k + 1)
                    _A[i1, p] = float(-alpha[i])
                    # x_ik
                    i1 = (i + 1)*dim + (k + 1)
                    i2 = (k + 1)*dim + (i + 1)
                    _A[i1, p] = alpha[i]/2
                    _A[i2, p] = alpha[i]/2
                    for j in G.neighbors(i):
                        if (j, k) not in G.edges():    
                            # x_jk
                            i1 = (j + 1)*dim + (k + 1)
                            i2 = (k + 1)*dim + (j + 1)
                            _A[i1, p] = .5
                            _A[i2, p] = .5
                    l += 1
                    mleq += 1
                    p += 1       
                    if p == step:
                        A = hstack([A, _A])
                        b = np.vstack([b, _b])
                        _A = lil_matrix((dim**2, step))
                        _b = np.zeros((step, 1))
                        p = 0

        for i,j in itertools.combinations(G.nodes(), 2):
            # (4.1)
            if 4.1 not in to_delete:
                # x_ij
                i1 = (i + 1)*dim + (j + 1)
                i2 = (j + 1)*dim + (i + 1)
                _A[i1, p] = -.5 
                _A[i2, p] = -.5
                l += 1
                mleq += 1
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


    elif model_out.lower() == 'sdpnal':   # SDPNAL+ model creation
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

        Bt = lil_matrix((dim, 0))
        u = np.zeros((0, 1))

        _Bt = lil_matrix((dim, step))
        _u = np.zeros((step, 1))
        p1 = 0

        # If not given, compute alphas for each node
        if alpha is None:
            alpha = {}
            for i in G.nodes():
                neighbors = list(G.neighbors(i))
                if neighbors:
                    # Compute alpha of G[N(i)] by computing the max clique on the complement
                    alpha[i] = nx.graph_clique_number(nx.complement(G.subgraph(neighbors)))

        # (1.1)
        if 1.1 not in to_delete:
            for i in G.nodes():
                for j in G.neighbors(i):
                    # x_ij
                    ii, jj = sorted([i, j], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + jj + 1
                    _Bt[idx, p1] = _2*0.5
                p1 += 1
                l+=1
                if p1 == step:
                    Bt = hstack([Bt, _Bt])
                    u = np.vstack([u, _u])
                    _Bt = lil_matrix((dim, step))
                    _u = np.zeros((step, 1))
                    p1 = 0

        # (2.1)
        if 2.1 not in to_delete:
            for i in G.nodes():
                for k in G.neighbors(i):
                    neighbors_but_k = [n for n in G.neighbors(i) if n != k]
                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] = float(1 - alpha[i])
                    for j in neighbors_but_k:
                        if (j, k) not in G.edges():
                            # x_jk
                            jj, kk = sorted([k, j], reverse=True)
                            idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                            _Bt[idx, p1] = _2*.5
                    p1 += 1
                    l+=1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

        # (3.1)
        if 3.1 not in to_delete:           
            for i in G.nodes():
                not_neighbors = [n for n in G.nodes() if n != i and n not in G.neighbors(i)]
                for k in not_neighbors:
                    # x_kk
                    idx = ((k + 2)*(k + 1))/2 + k + 1
                    _Bt[idx, p1] =  float(-alpha[i])
                    # x_ik
                    ii, kk = sorted([k, i], reverse=True)
                    idx = ((ii + 2)*(ii + 1))/2 + kk + 1
                    _Bt[idx, p1] = _2*(alpha[i]/2)

                    for j in G.neighbors(i):
                        if (j, k) not in G.edges():    
                            # x_jk
                            jj, kk = sorted([k, j], reverse=True)
                            idx = ((jj + 2)*(jj + 1))/2 + kk + 1
                            _Bt[idx, p1] = _2*0.5
                    p1 += 1
                    l+=1
                    if p1 == step:
                        Bt = hstack([Bt, _Bt])
                        u = np.vstack([u, _u])
                        _Bt = lil_matrix((dim, step))
                        _u = np.zeros((step, 1))
                        p1 = 0

        
        # Add constraint: X_00 = 1
        _At[0, p] = 1.
        _b[p] = 1.
        p += 1
        l+=1

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
            l +=1
            if p == step:
                At = hstack([At, _At])
                b = np.vstack([b, _b])
                _At = lil_matrix((dim, step))
                _b = np.zeros((step, 1))
                p = 0
        
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
        
        io.savemat(os.path.join(model_out_dir, filename) + '.mat', {'At' : At, 'b' : b, 'Bt' : Bt, 'u' : u, 'C' : C, 'L' : 0., 's' : float(n + 1)})
        # Exporting the model to output file
        if debug:
            print('Saving the model at: ' + os.path.join(model_out_dir, filename) + '.mat')
    
