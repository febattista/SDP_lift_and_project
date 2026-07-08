# Generalize M+ lifting to M_+(K,K), vectorize it, and export the M-operator LP as a solver-agnostic dict

## Context

`m_plus_lifting` (src/pyModules/SDPLifting.py:424, core in `_compute_m_plus_model` at :267) currently implements the Lovász–Schrijver M+ operator: each LP constraint of K is multiplied by the bound forms `x_i ≥ 0` and `1 − x_i ≥ 0` only. Goals of this change:

1. **Generalize** to the full `M_+(K,K)` operator, where every pair of constraints defining K (LP rows *and* variable bounds) is multiplied together.
2. **Speed up** constraint generation — the current builder appends per-nonzero into Python lists with `push`/`step` batching (SDPLifting.py:299-380), which is the hot path for large graphs.
3. **Export the M operator** (same variables/constraints, no PSD restriction) from `SDPModel` as a **solver-agnostic LP dict** — same style as the existing `to_adal`/`to_sdpnal`/`to_mosek` dict exports — that any solver front-end (Gurobi, HiGHS, …) can consume. `SDPModel` gains no solver dependency.

**Revisions vs. the previous plan (per user):**
- The export returns a plain dict describing the LP (objective, `num_cols`, `num_rows`, A, rhs, senses, bounds, …), **not** a `gurobipy.Model`. Gurobi-specific code lives only in the new `scripts/solve_gurobi.py`.
- **Cut classes are out of scope.** They depend on both the input LP and the lifting operator (currently `'ls'`); generalizing them for the new `'kk'` mode is deferred. The only requirement now is that the existing `'ls'` cut classes (labels 1–4, `map_edge/clique/nod_constraints` helpers, MATLAB pipeline) keep working unchanged. In `'kk'` mode, `cut_classes` is set to `None`.

### Bug found during planning (prerequisite fix)

The lifted inequality rows are mathematically valid only in the `≤ u` direction, but `_compute_m_plus_model` labels them `'>'` (SDPLifting.py:414). Proof: class 4 comes from `(1−x_i)(1−x_j) ≥ 0 ⇔ X_ii + X_jj − X_ij ≤ 1`; the code emits it with sense `'>'`, rhs 1, so the empty stable set's moment matrix `Y = e₀e₀ᵀ` (which gives 0) is cut off. Same inversion for classes 1–3. The SDPNAL `.mat` path is correct **only if** the MATLAB script passes `u` as an upper bound (`l ≤ B(X) ≤ u` with `l = []`); the `to_mosek` path (SDPModel.py:231-232 → solve_mosek.py:97-100) and `to_adal` enforce the wrong direction today. This must be fixed before the LP export can be correct.

**User confirmations still pending** (asked, no response — proceeding with recommendations): replace the loop builder with the vectorized one; K×K default = all distinct pairs, no squares.

**STATUS (2026-07-08): Section 0 (sense fix) is DONE.** Audit results:
- MATLAB confirmed correct: `sdpnalplus(blk,AA,C,b,L,U,BB,l,u,...)` is called with `l=[], u=cuts_u` (kelley_cutting_plane.m:105), i.e. `B(X) ≤ u`, and the violation test `AXfun(...) - u > epsilon` matches. Regenerated triangle `.mat` verified identical to pre-fix output (all arrays equal).
- `to_adal` audit: `ADMMsolver.ADMM_3b` treats **every** row as an equality and ignores `mleq` — nothing to negate. Documented in the `to_adal` docstring: only equality-only models (Theta_SDP) may be fed to the internal ADMM.
- `Theta_plus_SDP` audit: its `'>'` rows encode `X_ij ≥ 0` on non-edges — mathematically correct θ+ nonnegativity; only the inverted comment was fixed. Verified: MOSEK gives θ+(C5) = √5.
- Lifted M+ rows now labeled `'<'`; `'>'`-count sites fixed (solve_mosek.py, m_plus_lifting summary, model_building.py). Integral-feasibility invariant test added (triangle + C5, all stable sets, M+/Theta/Theta+). MOSEK post-fix: triangle M+ bound 1.0 (was 3.0), C5 M+ bound 2.0.

## Design

### 0. Sense fix

- In `_compute_m_plus_model`, label lifted inequality rows `'<'` (row values/rhs unchanged: they are already in ≤ form). Update the `cut_classes` docstring in SDPModel.py:40 ("aligned with '>' rows").
- `to_sdpnal` (SDPModel.py:87) needs no change — it splits on `('>' | '<')` already; `.mat` output stays byte-identical, so `cut_classes` semantics on the MATLAB side are untouched.
- `to_mosek` (SDPModel.py:162) and `solve_mosek.py` are sense-generic and become correct automatically.
- Audit `to_adal` (SDPModel.py:120): check what direction `pyModules/ADMMsolver.py` assumes for the first `mleq` rows; negate those rows on export if it assumes `≥`.
- Audit `Theta_plus_SDP` (SDPLifting.py:535) non-edge rows for the same inversion.
- Fix the stale expectation in `solve_mosek.py:158` and `m_plus_lifting`'s summary (SDPLifting.py:456) which count `'>'` rows — count non-`'='` rows instead.

### 1. Generalized builder: homogenized constraint forms + pair products

Represent every constraint of K in nonnegative homogenized form `g = (β, −α) ∈ R^{1+n}`, meaning `β − α·x ≥ 0` on `y = (1, x)`:
- LP row `a·x ≤ b` → `(b, −a)` (parser `read_constr_from_lp` at SDPLifting.py:92 already yields `(coeff_dict, b)`; assemble directly into a CSR matrix `G` of shape `(m_lp + 2n) × (n+1)`).
- Bound `x_i ≥ 0` → `e_i`; bound `1 − x_i ≥ 0` → `(1, −e_i)`.

For a pair `(k, l)` the product is `ĝ_kᵀ Y ĝ_l ≥ 0` with `Y = [1 xᵀ; x X]`, i.e. `⟨sym(ĝ_k ĝ_lᵀ), Y⟩ ≥ 0`. Emitted in the model's ≤ form: `⟨−sym(ĝ_k ĝ_lᵀ), Y⟩ ≤ 0`. Packed lower-tri entries under the canonical "implicit ×2 off-diagonal" convention (SDPModel.py:24): off-diag `(i>j)`: `(g_k[i]g_l[j] + g_k[j]g_l[i])/2`; diag: `g_k[i]g_l[i]`.

To reproduce the current model *exactly* in LS mode, apply the same two mechanical rewrites the current code does, as a vectorized post-transform (`use_diag_identity=True`): fold packed column `(i,0)` entries into `(i,i)` (valid under equality `X_0i = X_ii`), and move the `X_00` coefficient to the rhs (valid under `X_00 = 1`).

**Pair-selection modes** (new `lift_mode` parameter on `_compute_m_plus_model` / `m_plus_lifting`, default `'ls'` = today's behavior):
- `'ls'`: pairs = (each LP row) × (each of the 2n bound forms), plus bound×bound pairs `(1−x_i)x_j` and `(1−x_i)(1−x_j)` when `lift_bounds=True` — exactly classes 1–4, same ordering, same `cut_classes` values as today. `skip_func` keeps its `(constr, i)` signature for back-compat (`lovasz_schrijver_filter` at SDPLifting.py:174 still works, model_building.py:167 untouched). The `map_*_constraints` helpers stay as they are.
- `'kk'`: all distinct pairs `k < l` over the extended system {LP rows ∪ bounds}: LP×LP, LP×bounds, bounds×bounds. Squares `k = l` skipped by default (implied by `Y ⪰ 0`), available via `include_squares=True` — relevant when the PSD constraint is dropped (M-operator LP). Optionally skip the `x_i·x_j ≥ 0` pairs (`X_ij ≥ 0`) via a flag since the SDPNAL pipeline imposes them through `L=0`; default is to include them so the exported M LP is self-contained. **`cut_classes = None` in this mode** — designing class labels for K×K products is deferred until the cut-selection machinery is generalized (it depends on the input LP and lifting operator).
- Optional dedup of identical product rows via the existing `sp_unique` (SDPLifting.py:128), off by default.

Row count in `'kk'` mode is `O((m+2n)²)` — document this and keep `'ls'` the default.

### 2. Vectorized construction (replaces the per-nonzero Python loops)

Replace the append/`push` machinery (SDPLifting.py:163, :288-380) with sparse block algebra:

1. Build `G` (CSR) once from the parsed constraints.
2. Materialize the selected pair list as two index arrays `(K_idx, L_idx)`.
3. Compute the row-wise Kronecker product of `G[K_idx]` and `G[L_idx]` directly from CSR internals: per-pair nnz = `nnz_k · nnz_l`; output `indptr` by `cumsum`; `indices`/`data` via `np.repeat`/tile arithmetic on the CSR arrays (standard ragged-repeat technique, no Python loop over pairs). This yields COO triplets in full `(n+1)²` column space.
4. Fold full columns `i·dim + j` to packed `mat_idx(i,j)` with a vectorized `np.where(i>=j, i(i+1)//2+j, j(j+1)//2+i)`, weight `0.5` for `i≠j` entries and `1.0` for diagonal, negate, and let `coo_matrix(...).tocsr()` sum the symmetric duplicates.
5. Apply the diag-identity fold (column-0 → diagonal, `X_00` → rhs) as an index remap on the same triplet arrays.
6. Equality rows (`X_00 = 1`, `X_0i = X_ii`, SDPLifting.py:382-402) built the same way in one shot.

Total work and memory are `O(Σ_pairs nnz_k·nnz_l)` in numpy instead of interpreted Python — for stable-set LPs supports are tiny (edges: 3 nnz homogenized), so this is orders of magnitude faster and the `step`/`SDP_LIFTING_STEP` batching becomes unnecessary. Keep the `step` parameter accepted (used as an optional pair-chunk size to bound peak memory on huge instances; note the changed meaning in scripts/parameters.py:80).

The wrapper `m_plus_lifting` keeps its signature + printing + `.mat` writing unchanged; in `'ls'` mode `cut_classes` is produced exactly as today.

### 3. Solver-agnostic M-operator LP export from `SDPModel`

The M operator is exactly the current model with the packed entries of `Y` as free scalar LP variables. Add to `SDPModel` — **no solver imports**, dict export in the same style as `to_adal`/`to_sdpnal`/`to_mosek`:

- `_lp_scale()`: vector with 2 on off-diagonal packed positions, 1 on diagonal (sibling of `_svec_scale`, SDPModel.py:53) — converts canonical rows to plain linear coefficients, since `⟨A,X⟩ = Σ A_jj X_jj + 2 Σ_{j>k} A_jk X_jk`.
- `to_lp(nonneg=False)` → returns a dict describing the LP over the `packed_dim` scalar variables `y_k = X_{ij}` (k = `mat_idx(i,j)`):

  | key | value |
  |---|---|
  | `num_cols` | `packed_dim` |
  | `num_rows` | `self.num_rows` |
  | `obj` | length-`num_cols` dense array: packed lower-tri `C` entries × `_lp_scale()` (off-diag ×2) |
  | `obj_sense` | `self.obj_sense` (1.0 = min, −1.0 = max) |
  | `obj_offset` | `self.obj_offset` |
  | `A` | `(num_rows × num_cols)` CSR: `self.A.multiply(_lp_scale())` |
  | `rhs` | `self.row_rhs` copy |
  | `senses` | `self.row_senses` copy (`'<'`/`'>'`/`'='`) |
  | `lb` | length-`num_cols` array: 0.0 if `nonneg` (LP analog of SDPNAL's `L=0`) else −inf |
  | `ub` | length-`num_cols` array: +inf |
  | `col_names` | `['X_i_j', ...]` from `_packed_ij()` (SDPModel.py:70) |
  | `row_names` | `self.row_names` (may be empty) |

  Documented in the docstring like the other exports. Any LP solver front-end can be built from this dict without touching `SDPModel`.
- New script `scripts/solve_gurobi.py` mirroring `scripts/solve_mosek.py`: build the model from a `.lp` (LS or K×K mode), call `to_lp()`, assemble a `gurobipy.Model` from the dict (`addMVar(num_cols, lb, ub)`, `setObjective(obj @ y + obj_offset)`, one `addMConstr` per sense mask), solve, print the bound (bound = −objective, same as solve_mosek.py:168). An optional `--write-lp` flag writes the `.lp`/`.mps` file via Gurobi from the same assembled model. `gurobipy` is imported only in this script.

### Plan persistence (first step of implementation)

Copy this plan to `docs/plans/m-plus-kk-lifting.md` in the repo and commit it, so it survives across sessions and is visible to collaborators. Keep it updated as decisions are made during implementation. Also add a memory pointer (plan path + status) to the project memory index.

### Files to modify

- `src/pyModules/SDPLifting.py` — sense fix; new `G`-based vectorized builder; `lift_mode='ls'|'kk'` + `include_squares`; remove `push`; keep `mat_idx`, parser, wrappers; `'ls'` cut_classes preserved exactly, `'kk'` sets `cut_classes=None`.
- `src/pyModules/SDPModel.py` — `_lp_scale`, `to_lp` (dict export, no solver deps); `cut_classes` docstring; `to_adal` direction audit.
- `scripts/solve_mosek.py` — ineq counting; sanity after sense fix.
- `scripts/solve_gurobi.py` — new; only place with a Gurobi dependency for LP solving.
- `scripts/model_building.py` — pass-through for `lift_mode` (default unchanged).
- `scripts/parameters.py` — note on `SDP_LIFTING_STEP`'s new meaning.
- `tests/test_sdp_lifting.py` — new tests below.

### Tests

1. **Integral-feasibility invariant** (catches the sense bug and any future sign error): for every stable set S of the triangle (and a 5-cycle), build `y = (1, χ_S)`, `Y = yyᵀ`, and assert every model row is satisfied under its sense — for both `'ls'` and `'kk'` modes.
2. **Regression/equivalence**: vectorized `'ls'` builder reproduces the pre-rewrite model exactly (same `A` rows, rhs, senses, **and `cut_classes`**) on the triangle LP and a small clique LP — implement by snapshotting the old builder's output before deleting it. This is the guarantee that the existing cut-class pipeline still works.
3. **K×K hand-check**: for one LP row pair on a 2-variable instance, compare the generated packed row against the manually expanded product.
4. **LP-dict export**: on the triangle FRAC model, check `to_lp()` shape/keys consistency (`A.shape == (num_rows, num_cols)`, off-diag obj/constraint coefficients doubled vs. packed `C`/`A`, `nonneg=True` sets `lb=0`); then assemble a Gurobi model from the dict and solve — the M-operator LP bound must be 1.0 (N(K) of an odd cycle satisfies the odd-cycle inequality). Skip the solve part if `gurobipy` is unavailable.

### Verification (end-to-end)

- `pytest tests/` — all old + new tests.
- `python scripts/solve_mosek.py <triangle .lp>` — after the sense fix, the M+ bound must be ≈ 1.0 (today it would report 3.0, confirming the bug).
- `python scripts/solve_gurobi.py <triangle .lp>` — M-operator LP bound 1.0; repeat on a C₅ FRAC LP (expected N-operator bound 2.0) with both `'ls'` and `'kk'` modes, checking `'kk'` ≤ `'ls'` bound.
- Regenerate one small `.mat` and diff against the pre-change output to confirm the SDPNAL pipeline (including `cut_classes`) is unaffected.
