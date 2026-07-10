# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from dataclasses import dataclass, field
from scipy import sparse
from scipy.io import savemat
from scipy.sparse import csr_matrix


@dataclass
class SDPModel:
    """
    Standard representation of a semidefinite program (SDP).

        min/max  <C, X>
        s.t.     <A_i, X>  ['<' | '>' | '=']  row_rhs_i,  i = 1..num_rows
                 X >= 0 (PSD)
                 L <= X <= U (elementwise, optional)

    All constraint matrices are stored in canonical lower-triangular form:
    A[i, mat_idx(j,k)] holds the symmetric entry (A_i)_{jk} for j >= k. 
    The inner product is defined as:

        <A_i, X> = sum_j (A_i)_{jj} X_{jj} + 2 * sum_{j>k} (A_i)_{jk} X_{jk}

    Row senses : '<' (<=),  '>' (>=),  '=' (==)
    """
    dim:        int                # dimension of PSD variable X (e.g. n+1)
    n:          int                # number of graph nodes
    num_rows:   int                # total constraints (equalities + inequalities)

    C:          csr_matrix         # (dim x dim) sparse lower-tri objective matrix
    obj_sense:  float = 1.0        # 1.0 = minimize, -1.0 = maximize
    obj_offset: float = 0.0

    A:          csr_matrix = None  # (num_rows x packed_dim) constraint matrix
    row_rhs:    np.ndarray = None  # (num_rows,)
    row_senses: np.ndarray = None  # (num_rows,)  dtype='U1': '<', '>', '='

    # Elementwise bounds on X: None (absent), a scalar broadcast to every
    # entry, or a (dim x dim) array with -inf/+inf on unbounded entries.
    # SDPNAL+ takes them natively (L <= X <= U); other exports convert
    # (to_mosek: explicit rows, to_lp: variable bounds).
    L:          object = None
    U:          object = None

    cut_classes: np.ndarray = None  # class labels aligned with '<' rows (M+ 'ls' mode only; None in 'kk' mode)
    col_names:  list = field(default_factory=list)
    row_names:  list = field(default_factory=list)

    @property
    def packed_dim(self):
        """Number of entries in the packed lower-triangular representation."""
        return (self.dim * (self.dim + 1)) // 2

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _svec_scale(self):
        """Build the svec scaling vector: sqrt(2) for off-diagonal, 1 for diagonal."""
        dim = self.dim
        s = np.full(self.packed_dim, np.sqrt(2))
        for i in range(dim):
            s[i * (i + 1) // 2 + i] = 1.0   # mat_idx(i, i) = i*(i+1)//2 + i
        return s

    def _lp_scale(self):
        """Build the LP scaling vector: 2 for off-diagonal, 1 for diagonal.

        Converts canonical packed rows to plain linear coefficients over the
        scalar variables y_k = X_{ij}, since
        <A, X> = sum_j A_jj X_jj + 2 * sum_{j>k} A_jk X_jk.
        """
        s = np.full(self.packed_dim, 2.0)
        for i in range(self.dim):
            s[i * (i + 1) // 2 + i] = 1.0   # mat_idx(i, i) = i*(i+1)//2 + i
        return s

    def _bounds_packed(self):
        """Broadcast L/U to packed arrays (lb, ub) over the packed entries.

        None -> -inf/+inf everywhere; scalar -> constant; (dim x dim)
        array-like -> entry (i, j) with i >= j in packed order.
        """
        def expand(bound, fill):
            if bound is None:
                return np.full(self.packed_dim, fill)
            b = np.asarray(bound, dtype=float)
            if b.ndim == 0:
                return np.full(self.packed_dim, float(b))
            pi, pj = self._packed_ij()
            return b[pi, pj]

        return expand(self.L, -np.inf), expand(self.U, np.inf)

    def _C_full(self):
        """Unpack lower-tri C to a full (dim x dim) symmetric sparse matrix."""
        C_coo = self.C.tocoo()
        off = C_coo.row != C_coo.col
        rows = np.concatenate([C_coo.row, C_coo.col[off]])
        cols = np.concatenate([C_coo.col, C_coo.row[off]])
        data = np.concatenate([C_coo.data, C_coo.data[off]])
        return sparse.coo_matrix((data, (rows, cols)), shape=self.C.shape).tocsr()

    def _packed_ij(self):
        """Return arrays (packed_rows, packed_cols) mapping packed index k to (i, j) with i >= j."""
        dim = self.dim
        rows = np.empty(self.packed_dim, dtype=int)
        cols = np.empty(self.packed_dim, dtype=int)
        k = 0
        for i in range(dim):
            for j in range(i + 1):
                rows[k] = i
                cols[k] = j
                k += 1
        return rows, cols

    # ------------------------------------------------------------------
    # Solver-specific exports
    # ------------------------------------------------------------------

    def to_sdpnal(self):
        """Convert to SDPNAL+ format.

        Applies the svec sqrt(2) scaling to off-diagonal packed entries of At
        and Bt. C is returned as a full (dim x dim) symmetric matrix (no extra
        scaling: SDPNAL+ reads C as a plain matrix and computes trace(C X) directly).

        Elementwise bounds are passed through as-is (scalar or matrix):
        SDPNAL+ accepts L <= X <= U natively, so they are not converted
        into inequality rows.

        Returns
        -------
        dict with keys: At, b, Bt (or None), u (or None), C, s, L, U,
        cut_classes.
        """
        s = self._svec_scale()

        eq_mask   = self.row_senses == '='
        ineq_mask = (self.row_senses == '>') | (self.row_senses == '<')

        # Split by sense first, then scale each subset independently.
        # This avoids creating a full-size copy of A with scaling applied.
        At = self.A[eq_mask].multiply(s).T.tocsc()
        b  = self.row_rhs[eq_mask].reshape(-1, 1)
        Bt = self.A[ineq_mask].multiply(s).T.tocsc() if ineq_mask.any() else None
        u  = self.row_rhs[ineq_mask].reshape(-1, 1) if ineq_mask.any() else None

        return {
            'At': At,
            'b':  b,
            'Bt': Bt,
            'u':  u,
            'C':  self._C_full(),
            's':  float(self.dim),
            'L':  self.L,
            'U':  self.U,
            'cut_classes': self.cut_classes,
        }

    def to_adal(self):
        """Convert to ADAL format.

        Unpacks each constraint row from packed lower-triangular to a full dim^2
        flat vector. Inequalities are placed first (count = mleq), then equalities.

        Note: the internal ADMM solver (ADMMsolver.ADMM_3b) treats every row as
        an equality and ignores mleq, so only equality-only models (e.g.
        Theta_SDP) may be passed to it. mleq is exported for ADAL-format
        compatibility; inequality rows keep the direction given by row_senses.

        Elementwise bounds: ADMM_3b always enforces X >= 0, which is exactly 
        L = 0 with U absent — such bounds need no conversion. 
        Any other finite bound (L_ij != 0 or a finite U_ij) becomes an explicit 
        inequality row in '<' form (lower bounds negated), placed after the 
        model's inequality rows and counted in mleq.

        Returns
        -------
        dict with keys: A (num_rows x dim^2), b, C (dim x dim full), mleq.
        """
        dim = self.dim
        pi, pj = self._packed_ij()   # packed index k maps to (pi[k], pj[k]) with pi >= pj

        # Reorder rows: inequalities first (mleq), then equalities
        ineq_mask = (self.row_senses == '>') | (self.row_senses == '<')
        eq_mask   = self.row_senses == '='
        mleq      = int(ineq_mask.sum())

        order  = np.concatenate([np.where(ineq_mask)[0], np.where(eq_mask)[0]])
        A_ord  = self.A[order].tocoo()
        rhs    = self.row_rhs[order].reshape(-1, 1)

        # Expand each packed entry to the full dim^2 flat index (both triangles)
        ri, ci, di = [], [], []
        for r, c, v in zip(A_ord.row, A_ord.col, A_ord.data):
            i, j = pi[c], pj[c]
            ri.append(r); ci.append(i * dim + j); di.append(v)
            if i != j:
                ri.append(r); ci.append(j * dim + i); di.append(v)

        nrows  = A_ord.shape[0]
        A_full = sparse.coo_matrix(
            (di, (ri, ci)), shape=(nrows, dim * dim)
        ).tocsr()

        # Bounds differing from the built-in X >= 0 become '<' rows
        lb, ub = self._bounds_packed()
        lo = np.where(np.isfinite(lb) & (lb != 0.))[0]
        up = np.where(np.isfinite(ub))[0]
        if len(lo) or len(up):
            k    = np.concatenate([lo, up])
            sign = np.concatenate([np.full(len(lo), -1.), np.ones(len(up))])
            brhs = np.concatenate([-lb[lo], ub[up]])
            i, j = pi[k], pj[k]
            off  = i != j
            r = np.concatenate([np.arange(len(k)), np.where(off)[0]])
            c = np.concatenate([i * dim + j, (j * dim + i)[off]])
            v = np.concatenate([np.where(off, 0.5, 1.0) * sign,
                                (0.5 * sign)[off]])
            B_full = sparse.coo_matrix(
                (v, (r, c)), shape=(len(k), dim * dim)).tocsr()
            A_full = sparse.vstack(
                [A_full[:mleq], B_full, A_full[mleq:]], format='csr')
            rhs   = np.vstack([rhs[:mleq], brhs.reshape(-1, 1), rhs[mleq:]])
            mleq += len(k)

        return {
            'A':    A_full,
            'b':    rhs,
            'C':    self._C_full(),
            'mleq': mleq,
        }

    def to_mosek(self):
        """Convert to MOSEK Task API format (bar-variable SDP interface).

        MOSEK's lower-triangle inner product convention matches the canonical
        format exactly, so no sqrt(2) scaling is applied.  A lower-triangle
        entry (j, k, val) with j > k is interpreted by MOSEK as
        A[j,k] = A[k,j] = val, giving
            trace(A X) = sum_j val_jj X_jj + 2 * sum_{j>k} val_jk X_jk
        which is the canonical inner product.

        Elementwise bounds L/U: MOSEK bar-variables cannot carry entry
        bounds (the putvarbound* family only addresses scalar variables),
        so every finite bound is converted into an explicit row
        <E_ij, X> >= L_ij (sense '>') or <= U_ij (sense '<'), appended
        after the model rows — num_rows in the returned dict counts them
        and may exceed self.num_rows.

        Typical usage with the MOSEK Python API:

            import mosek
            d = model.to_mosek()
            with mosek.Task() as task:
                task.appendcons(d['num_rows'])
                task.appendbarvars([d['dim']])

                n_c = len(d['C_val'])
                task.putbarcblocktriplet(n_c, [0]*n_c,
                                         d['C_rows'], d['C_cols'], d['C_val'])

                n_a = len(d['A_val'])
                task.putbarablocktriplet(n_a, d['A_cons'], [0]*n_a,
                                         d['A_rows'], d['A_cols'], d['A_val'])

                bk = [mosek.boundkey.fx if s == '=' else
                      mosek.boundkey.lo if s == '>' else
                      mosek.boundkey.up for s in d['row_senses']]
                task.putconboundlist(range(d['num_rows']), bk, d['blc'], d['buc'])

                sense = (mosek.objsense.minimize if d['obj_sense'] > 0
                         else mosek.objsense.maximize)
                task.putobjsense(sense)
                task.optimize()

        Returns
        -------
        dict with keys:

        dim        : int, dimension of the PSD variable
        num_rows   : int, total constraints including converted bound rows
        obj_sense  : float, 1.0 = minimize, -1.0 = maximize
        C_rows     : int32 array, row indices of objective lower-tri triplets (row >= col)
        C_cols     : int32 array, col indices of objective lower-tri triplets
        C_val      : float array, values of objective lower-tri triplets
        A_cons     : int32 array, constraint indices (0-based)
        A_rows     : int32 array, PSD matrix row indices (row >= col)
        A_cols     : int32 array, PSD matrix col indices
        A_val      : float array, constraint matrix values
        row_senses : char array, '=', '>', or '<' per constraint
        blc        : float array, lower bound per constraint (rhs for '='/'>', -inf for '<')
        buc        : float array, upper bound per constraint (rhs for '='/'<', +inf for '>')
        """
        pi, pj = self._packed_ij()   # pi[k] >= pj[k] by construction

        # Objective: lower-triangular entries of C
        C_coo  = self.C.tocoo()
        C_rows = C_coo.row.astype(np.int32)
        C_cols = C_coo.col.astype(np.int32)
        C_val  = C_coo.data.copy()

        # Constraint matrices: unpack each packed column index to (row, col) in the PSD matrix
        A_coo  = self.A.tocoo()
        A_cons = A_coo.row.astype(np.int32)
        A_rows = pi[A_coo.col].astype(np.int32)   # pi[k] >= pj[k]
        A_cols = pj[A_coo.col].astype(np.int32)
        A_val  = A_coo.data.copy()

        # Constraint bounds: -inf / +inf on the unconstrained side
        blc = np.where(self.row_senses != '<', self.row_rhs, -np.inf)
        buc = np.where(self.row_senses != '>', self.row_rhs,  np.inf)

        # Elementwise bounds -> explicit rows (0.5 off-diagonal under the
        # canonical implicit x2, so each row's value is exactly X_ij)
        row_senses = self.row_senses
        num_rows   = self.num_rows
        lb, ub = self._bounds_packed()
        for vals, sense in ((lb, '>'), (ub, '<')):
            k = np.where(np.isfinite(vals))[0]
            if len(k) == 0:
                continue
            A_cons = np.concatenate(
                [A_cons, num_rows + np.arange(len(k))]).astype(np.int32)
            A_rows = np.concatenate([A_rows, pi[k]]).astype(np.int32)
            A_cols = np.concatenate([A_cols, pj[k]]).astype(np.int32)
            A_val  = np.concatenate(
                [A_val, np.where(pi[k] == pj[k], 1.0, 0.5)])
            row_senses = np.concatenate(
                [row_senses, np.full(len(k), sense, dtype='U1')])
            blc = np.concatenate(
                [blc, vals[k] if sense == '>' else np.full(len(k), -np.inf)])
            buc = np.concatenate(
                [buc, vals[k] if sense == '<' else np.full(len(k), np.inf)])
            num_rows += len(k)

        return {
            'dim':        self.dim,
            'num_rows':   num_rows,
            'obj_sense':  self.obj_sense,
            'C_rows':     C_rows,
            'C_cols':     C_cols,
            'C_val':      C_val,
            'A_cons':     A_cons,
            'A_rows':     A_rows,
            'A_cols':     A_cols,
            'A_val':      A_val,
            'row_senses': row_senses,
            'blc':        blc,
            'buc':        buc,
        }

    def to_lp(self):
        """Convert to a solver-agnostic LP dict (the M-operator relaxation).

        Drops the PSD requirement on X: the packed lower-triangular entries
        become scalar LP variables y_k = X_{ij}, k = mat_idx(i, j), and
        each canonical row turns into a plain linear constraint via the
        off-diagonal x2 scaling of _lp_scale().  Any LP solver front-end can
        be assembled from the returned dict (see scripts/solve_gurobi.py);
        no solver is imported here.

        Elementwise bounds L/U map directly onto the variable bounds lb/ub
        (the packed entries are genuine scalar variables here, so no
        inequality rows are needed).

        Note: for the LP to be M(K,K) proper, the model must be built with
        include_squares=True — the squared-constraint products are implied
        by PSD but not by the pairwise products alone.

        Returns
        -------
        dict with keys:

        num_cols   : int, number of LP variables (= packed_dim)
        num_rows   : int, number of linear constraints
        obj        : float array (num_cols,), objective coefficients
        obj_sense  : float, 1.0 = minimize, -1.0 = maximize
        obj_offset : float, constant objective offset
        A          : (num_rows x num_cols) CSR constraint matrix
        rhs        : float array (num_rows,)
        senses     : char array (num_rows,), '<', '>', or '='
        lb, ub     : float arrays (num_cols,), variable bounds
        col_names  : list of 'X_i_j' strings, packed order
        row_names  : list (may be empty)
        """
        s = self._lp_scale()

        obj = np.zeros(self.packed_dim)
        C_coo = self.C.tocoo()   # stored lower-tri: row >= col
        k = C_coo.row * (C_coo.row + 1) // 2 + C_coo.col
        obj[k] = C_coo.data * s[k]

        pi, pj = self._packed_ij()
        col_names = ['X_%d_%d' % (i, j) for i, j in zip(pi, pj)]

        lb, ub = self._bounds_packed()

        return {
            'num_cols':   self.packed_dim,
            'num_rows':   self.num_rows,
            'obj':        obj,
            'obj_sense':  self.obj_sense,
            'obj_offset': self.obj_offset,
            'A':          self.A.multiply(s).tocsr(),
            'rhs':        self.row_rhs.copy(),
            'senses':     self.row_senses.copy(),
            'lb':         lb,
            'ub':         ub,
            'col_names':  col_names,
            'row_names':  list(self.row_names),
        }

    def write_sdpnal_mat(self, out_path, do_split=False):
        """Serialize to SDPNAL+ .mat file(s).

        Elementwise bounds L/U are written as-is when set (SDPNAL+ takes
        them natively as bounds on the matrix variable) and omitted when
        None.

        Parameters
        ----------
        out_path : str — output path without extension (e.g. '.../graph_edge')
        do_split : bool — if True, split Bt across 4 files (for very large models)
        """
        d   = self.to_sdpnal()
        At  = d['At'];  b  = d['b']
        Bt  = d['Bt'];  u  = d['u']
        C   = d['C']
        cut_classes = d['cut_classes']

        base = {'At': At, 'b': b}
        if self.L is not None:
            base['L'] = self.L
        if self.U is not None:
            base['U'] = self.U
        base['C'] = C
        base['s'] = float(self.n + 1)
        if Bt is not None:
            base['u'] = u
        if cut_classes is not None:
            base['cut_classes'] = cut_classes

        if not do_split:
            if Bt is not None:
                base['Bt'] = Bt
            savemat(out_path + '.mat', base)
        else:
            savemat(out_path + '_1.mat', base)   # Bt intentionally excluded from file 1
            splt = int(np.ceil(Bt.shape[1] / 4))
            savemat(out_path + '_2.mat', {'Bt_1': Bt[:, :splt]})
            savemat(out_path + '_3.mat', {'Bt_2': Bt[:, splt:2*splt]})
            savemat(out_path + '_4.mat', {'Bt_3': Bt[:, 2*splt:3*splt]})
            savemat(out_path + '_5.mat', {'Bt_4': Bt[:, 3*splt:]})
