import os
_NUM_THREADS = "20"
os.environ["OMP_NUM_THREADS"] = _NUM_THREADS # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = _NUM_THREADS # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = _NUM_THREADS # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = _NUM_THREADS # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = _NUM_THREADS # export NUMEXPR_NUM_THREADS=6

import time
import numpy as np
import scipy.linalg as la
from sksparse.cholmod import cholesky
import scipy.sparse.linalg as sp_linalg
# from Node import Node


class Problem:

    def __init__(self, A=None, b=None, C=None, c=0.):
        self.A, self.b, self.C, self.c = A, b, C, c
        self.n = C.shape[0]

    def getProblem(self):
        return self.A, self.b, self.C, self.c

    def getSubproblem(self, node):
        res = []
        # Get the idxs of rows and cols to be deleted
        for k in node.fixed_vars:
            i, j = get_sub_idx(k, self.n)
            res = res + [i, j]
        res = np.array(res).reshape((1, -1))
        if not res.any():
            return self
        mask = np.ones(self.A.shape[1], dtype=bool)
        mask[res] = False
        # Select the submatrix A
        A_sub = self.A[:, mask]
        # Some rows (i.e. constraints) may become empty
        mask_empty_rows = A_sub.getnnz(1)>0
        A_sub = A_sub[mask_empty_rows]
        b_sub = self.b[mask_empty_rows]
        mask_C = np.ones(self.C.shape[0], dtype=bool)
        mask_C[node.fixed_vars] = False
        C_sub = self.C[mask_C, :]
        C_sub = C_sub[:, mask_C]

        return Problem(A=A_sub, b=b_sub, C=C_sub, c=node.c)

    def write(self, filename, node, G):
        S = G.subgraph(set(G.nodes()) - set([i - 1 for i in node.fixed_vars]))
        # Inverse permutation from current model indexes to original graph nodes
        ppam = {i: j for i, j in enumerate(S.nodes())}
        ppam[-1] = -1
        with open(filename, 'w') as f:
            obj = ''
            c_nnz = self.C.nonzero()
            for i in range(len(c_nnz[0])):
                obj += '%.1fx%s' % (self.C[c_nnz[0][i], c_nnz[0][i]], '(%s, %s)' % (ppam[c_nnz[0][i] - 1], ppam[c_nnz[0][i] - 1]))
            f.write('min ' + obj + '\n')
            f.write('s.t.\n')
            n = int(np.sqrt(self.A.shape[1]))
            for k in range(self.A.shape[0]):
                constr = {}
                curr_a = self.A[k, :].reshape((n, n)).todense()
                nnz = curr_a.nonzero()
                for l in range(len(nnz[0])):
                    i, j = sorted([nnz[0][l], nnz[1][l]])
                    if (i, j) in constr:
                        constr[(i, j)] += curr_a[i, j]
                    else:
                        constr[(i, j)] = curr_a[i, j]
                s = ''
                for i, j in constr:
                    s += '%.1fx%s' % (constr[(i, j)], (ppam[i - 1], ppam[j - 1]))
                s += ' = %.1f' % self.b[k]
                f.write(s + '\n')
            text = ''
            for n in node.node_pick:
                text += 'x_{} = {}, '.format(n[0], n[1])
            f.write(text + '\n')


class Solution:

    def __init__(self, primal=np.Inf, dual=-np.Inf, X=None, Z=None, S=None, safe_dual=-np.Inf, sigma=1.):
        self.primal, self.dual = primal, dual
        self.X, self.Z, self.S, self.sigma = X, Z, S, sigma
        self.safe_dual = safe_dual

    def getSubsolution(self, fixed_vars):
        # Get the corresponding X,Z,S which will be used as warmstart
        # from the optimal solution of the father node
        mask = np.ones(self.X.shape[0], dtype=bool)
        mask[fixed_vars] = False
        X_sub = self.X[mask, :]
        X_sub = X_sub[:, mask]
        Z_sub = self.Z[mask, :]
        Z_sub = Z_sub[:, mask]
        S_sub = None
        if self.S is not None:
            S_sub = self.S[mask, :]
            S_sub = S_sub[:, mask]

        return Solution(self.primal, self.dual,
                        X_sub, Z_sub, S=S_sub,
                        safe_dual=self.safe_dual, sigma=self.sigma)


class ADMMsolver:
    @staticmethod
    def solve(P, S=None, target=np.Inf, options={}, norm_bound=None, c=0, use_patience=False):
        if not options:
            options['tolerance'] = 1e-4
            options['max_iter'] = 1000000
            options['timelimit'] = 1800
            options['print_it'] = 1000
            options['debug'] = True
        if S:
            # A warmstart is given
            X, y, Z, S, safebounds, dual, primal, sigma, factor, norm_bound, pruned, secs = ADMM_3b(P.A, P.b, P.C, sigma=S.sigma, Y=S.X, Z=S.Z, S=S.S, options=options, target=target, norm_bound=norm_bound, c=c)
        else:
            # Otherwise the default call
            X, y, Z, S, safebounds, dual, primal, sigma, factor, norm_bound, pruned, secs = ADMM_3b(P.A, P.b, P.C, options=options, target=target, norm_bound=norm_bound, c=c, use_patience=use_patience)

        opt_sol = Solution(primal, dual, X, Z, S=S, safe_dual=np.max(safebounds), sigma=sigma)
        return opt_sol, norm_bound, pruned, secs


def get_sub_idx(i, dim):
    i = i + 1
    return np.arange((i - 1) * dim, i * dim), np.arange(dim) * dim + (i - 1)


def projf(X, L, U):
    if X.shape == ():
        return L if X < L else (U if X > U else X)
    # Works only with Numpy Data types
    idx_l = X < L
    idx_u = X > U
    X[idx_l] = L[idx_l]
    X[idx_u] = U[idx_u]
    return X


def safebound_error(dual, Aty, C, X, S, mu=1.1, max_lamb_x=None):
    n = C.shape[0]
    # If no upper bound on the maximum eigenvalue of the optimal X is given,
    # then it is estimated by the current X, scaled by a factor mu
    max_lamb_x = mu * np.max(la.eigh(X)[0]) if not max_lamb_x else max_lamb_x

    # Compute a feasible Znew (which in general will not be Positive Semidefinite)
    Znew = C - Aty - S
    lb0 = dual
    pert = 0.
    lam, _ = la.eigh(Znew)
    I = np.where(lam < 0)[0]

    # print(lam)
    if I.shape[0] > 0:
        mineig = np.min(lam[I])
        pert += max_lamb_x * np.sum(lam[I])
        #print('numneg_ev in Znew: %3.0d, mineigZnew : %3.2e, perturbation: %3.2e' % (I.shape[0], mineig, pert))
    lb = lb0 + pert
    #print('dual obj: ', dual, 'safebound: ', lb)
    return lb, max_lamb_x


def safebound(dual, K, U, debug=False, iter=0):
    bound = dual - U * K
    for i in range(iter):
        U = np.min([np.abs(np.floor(bound)) + 1, U])
        bound = dual - U * K
    if debug:
        print('Perturbation: %13.10f' % (-U*K))
        print('Norm Bound: %d' % U, flush=True)
        print('Safe dual bound: %13.4f' % bound, flush=True)
    return bound


def can_be_pruned(safe_bound, target):
    return safe_bound >= target


def ADMM_3b(A, b, C, sigma=1., Y=None, Z=None, S=None, norm_bound=None, options={}, target=np.Inf, c=0, use_patience=False):

    ############################################
    # TESTING
    ############################################
    initial_patience = 500
    patience = initial_patience
    tailOff = 1e-2
    old_bound, curr_bound = None, None

    ############################################
    # Initialization
    tstart = time.time()
    # Read option from dict
    tol = 1e-6 if 'tolerance' not in options else options['tolerance']
    max_iter = 10000 if 'max_iter' not in options else options['max_iter']
    timelimit = 3600 if 'timelimit' not in options else options['timelimit']
    print_it = 10 if 'print_it' not in options else options['print_it']
    debug = True if 'debug' not in options else options['debug']
    num_iter = 1
    done = False

    # initialization of sigma box
    t_min = np.float64(1e-4)
    t_max = np.float64(1e+7)

    m, n2 = A.shape
    n = int(np.sqrt(n2))
    At, AAT = A.T, A @ A.T

    factor = cholesky(AAT)
    secs = time.time() - tstart
    if debug:
        print('Cholesky factorization completed after: %12.5f' % secs)

    if Y is None:
        Y = np.zeros((n, n))
    if Z is None:
        Z = np.zeros((n, n))
    if S is None:
        S = np.zeros((n, n))

    nonnegative = True

    normb, normC = la.norm(b), sp_linalg.norm(C)
    if (not norm_bound) or (norm_bound == np.Inf):
        norm_bound = int(n)
    bound = -np.Inf
    err_bound = -np.Inf

    all_err_bound = []
    all_norm_bounds = []
    all_norms = []
    all_duals = []

    g = b - A * Y.reshape((-1, 1))

    if debug:
        print(' it     secs   safebound       dual        primal       lg(rrd)    lg(rrp)   sigma   ||X||-bound')

    while not done:
        # weight for sigma update
        w = np.power(2, -(num_iter - 1) / 100)

        # given Y, Z, S and sigma, solve for y
        rhs = A * (C.reshape((-1, 1)) - Z.reshape((-1, 1)) - S.reshape((-1, 1))) + g / sigma
        y = factor(rhs)
        Aty = (At * y).reshape((n, n))
        # now form W1 = W1(y, Y, sigma, Z)

        if nonnegative:
            W1 = Aty - C + Y / sigma + Z
            S = -W1
            # S = (S + S.T) / 2
            S[S < 0] = 0

        M = Aty - C + Y / sigma + S
        # M = (M + M.T) / 2
        lam, ev = la.eigh(M)
        I = np.where(lam > 0)[0]
        j = len(I)
        if j < n / 2:
            evp = np.zeros((n, j))
            for r in range(j):
                ic = I[r]
                evp[:, r] = ev[:, ic] * np.sqrt(lam[ic])
            if j == 0:
                evp = np.zeros((n, 1))
            Mp = evp @ evp.T
            Mn = M - Mp
        else:
            I = np.where(lam < 0)[0]
            j = len(I)  # should be <= n/2
            evn = np.zeros((n, j))
            for r in range(j):
                ic = I[r]
                evn[:, r] = ev[:, ic] * np.sqrt(-lam[ic])
            if j == 0:
                evn = np.zeros((n, 1))
            Mn = -evn @ (evn.T)
            Mp = M - Mn

        # Ensure these matrices are symmetric
        # Mn = (Mn + Mn.T) / 2
        # Mp = (Mp + Mp.T) / 2

        Z = -Mn
        X = sigma * Mp

        # Y update
        Y = X

        g = b - A * Y.reshape((-1, 1))
        G = C - Aty - Z - S

        err_d = la.norm(G)
        dual = (b.T @ y).item()
        err_p = la.norm(g)
        primal = np.sum(C.reshape((1, -1)) * Y.reshape((-1, 1)))
        rel_err_p, rel_err_d = err_p / (1 + normb), err_d / (1 + normC)

        secs = time.time() - tstart
        num_iter = num_iter + 1

        # Safe bound by weak duality, in general weaker than error bound but cheap computationally,
        # computed at every iteration of ADMM
        # It is needed to know an upper bound on the norm(X, 'fro'), X optimal solution of the (P)
        # For Stable Set if the graph is connected, we might assume norm(X, 'fro') <= floor(n/2) + 1
        new_bound = safebound(dual, err_d, norm_bound, iter=5, debug=False)
        #all_norm_bounds.append(new_bound)
        bound = np.max([new_bound, err_bound])
        norm_bound = np.floor(-bound) + 1 if np.floor(-bound) < norm_bound else norm_bound

        # Computing safe bounds
        # Error bound, stronger but expensive (eigh decomposition), computed every 100 iterations of ADMM
        # It is needed to know an upper bound on the max(eigh(X)), X optimal solution of the (P)
        if num_iter % 1000 == 0:
            new_err_bound, _ = safebound_error(dual, Aty, C, X, S, max_lamb_x=norm_bound)
            #all_err_bound.append(new_err_bound)
            err_bound = np.max([new_err_bound, err_bound])
            norm_bound = np.floor(np.abs(err_bound)) + 1 if np.floor(np.abs(err_bound)) < norm_bound else norm_bound

        # Printing
        if (num_iter % print_it) == 0:
            if debug:
                print('%3.0d %8.2f %13.5e %13.5e %13.5e %8.3f %8.3f  %9.6f  %d' %
                  (num_iter, secs, bound, dual, primal, np.log10(rel_err_d), np.log10(rel_err_p), sigma, norm_bound))


        ############################################
        # TESTING
        ############################################
        if use_patience:
            curr_bound = np.max([bound, err_bound])
            # patience = patience - 1 if old_bound and (curr_bound - old_bound < tailOff) and np.max([rel_err_d, rel_err_p]) < tailOff else initial_patience
            # if old_bound:
            #     print("%10.5f, %10.5f" % (curr_bound, old_bound))
            #     print("%13.5e" % (curr_bound - old_bound))
            if old_bound and (np.max([rel_err_d, rel_err_p]) < tol*10 or np.min([rel_err_d, rel_err_p]) > tol*1e-1) :
                improvement = curr_bound - old_bound
                if improvement > 0 and 1e-7 <= improvement < tailOff:
                    patience = patience - 1 
                else:
                    patience = initial_patience
            old_bound = curr_bound
        ############################################

        # Stopping criteria
        if np.max([rel_err_d, rel_err_p]) < tol or num_iter > max_iter or secs > timelimit or can_be_pruned(err_bound, target) or patience <= 0:
            if debug:
                print('%3.0d %8.2f %13.5e %13.5e %13.5e %8.3f %8.3f  %9.6f  %d' %
                  (num_iter, secs, bound, dual, primal, np.log10(rel_err_d), np.log10(rel_err_p), sigma, norm_bound))
            new_bound = safebound(dual, err_d, norm_bound, iter=5)
            #all_norm_bounds.append(new_bound)
            bound = np.max([new_bound, err_bound])
            norm_bound = np.floor(-bound) + 1 if np.floor(-bound) < norm_bound else norm_bound

            if err_bound == -np.Inf:
                new_err_bound, _ = safebound_error(dual, Aty, C, X, S, max_lamb_x=norm_bound)
                #all_err_bound.append(new_err_bound)
                err_bound = np.max([new_err_bound, err_bound])
                norm_bound = np.floor(np.abs(err_bound)) + 1 if np.floor(np.abs(err_bound)) < norm_bound else norm_bound

            if debug:
                print('total time: %10.3f' % secs)
            if num_iter > max_iter:
                if debug:
                    print('max outer iterations reached.')
            if secs > timelimit:
                if debug:
                    print('Time limit exceeded')

            if patience <= 0 and debug:
                print('No more patience!')

            pruned = True if can_be_pruned(err_bound, target) else False

            done = 1
        #all_norms.append(norm_bound)
        #all_duals.append(dual)
        ratio = la.norm(X.reshape((-1, 1))) / la.norm(Z.reshape((-1, 1)))
        sigma = (1 - w) * sigma + w * projf(ratio, t_min, t_max)
    # print('Dual obj: ', dual + c)
    # print('Safe bound: ', bound + c)
    # print('Norm bound: ', norm_bound - c)
    # print('-'*10, flush=True)
    return Y, y, Z, S, [bound, err_bound], dual, primal, sigma, factor, norm_bound, pruned, secs
