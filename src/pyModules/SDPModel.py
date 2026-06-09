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

    cut_classes: np.ndarray = None  # class labels aligned with '>' rows (M+ only)
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

        Returns
        -------
        dict with keys: At, b, Bt (or None), u (or None), C, s, cut_classes.
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
            'cut_classes': self.cut_classes,
        }

    def to_adal(self):
        """Convert to ADAL format.

        Unpacks each constraint row from packed lower-triangular to a full dim^2
        flat vector. Inequalities are placed first (count = mleq), then equalities.

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

        Typical usage with the MOSEK Python API:

            import mosek
            d = model.to_mosek()
            with mosek.Task() as task:
                task.appendcons(model.num_rows)
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
                task.putconboundlist(range(model.num_rows), bk, d['blc'], d['buc'])

                sense = (mosek.objsense.minimize if d['obj_sense'] > 0
                         else mosek.objsense.maximize)
                task.putobjsense(sense)
                task.optimize()

        Returns
        -------
        dict with keys:

        dim        : int, dimension of the PSD variable
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

        return {
            'dim':        self.dim,
            'obj_sense':  self.obj_sense,
            'C_rows':     C_rows,
            'C_cols':     C_cols,
            'C_val':      C_val,
            'A_cons':     A_cons,
            'A_rows':     A_rows,
            'A_cols':     A_cols,
            'A_val':      A_val,
            'row_senses': self.row_senses,
            'blc':        blc,
            'buc':        buc,
        }

    def write_sdpnal_mat(self, out_path, do_split=False):
        """Serialize to SDPNAL+ .mat file(s).

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

        base = {'At': At, 'b': b, 'L': 0., 'C': C, 's': float(self.n + 1)}
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
