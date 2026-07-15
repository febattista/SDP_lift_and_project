# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
"""Numerical experiment: solve the M(K,K) LP relaxation with Gurobi.

Builds each model in memory from the existing LP files (no dataset
duplication, no .mat serialization) with lift_mode='kk' + squares, and
solves the PSD-free LP via SDPModel.to_lp() using barrier without
crossover (see docs/plans/kk-lp-gurobi-solving.md for the measurements
behind that choice). The 7200 s time limit is a Gurobi parameter, so it
bounds the solve phase only; build time is recorded separately. Each
lifted M(K,K) LP is also written to results/.../lp/ before solving.

Instance selection (the manifest) lives in this file: smallDIMACS (all),
a DIMACS subset, and Random n in {150, 175, 200}; relaxations nod_gamma /
nod_theta / nod_alpha / edge / cov. Models whose estimated size exceeds
--rows-cap / --nnz-cap are recorded as SKIP rows with the estimates —
those are exactly the instances that need the deferred row-generation
path (dense-graph products have ~|support|^2 nonzeros per row).

Restartable: (instance, relax) pairs already present in the results CSV
are skipped, so it can be resumed after a crash or run incrementally.

Usage
-----
  python experiment_m_kk.py --dry-run          # size map, no solves
  python experiment_m_kk.py                    # full run
  python experiment_m_kk.py --only DSJC125     # filter by substring
  nohup python experiment_m_kk.py > run.log 2>&1 &      # server
"""

import sys, os, csv, time, argparse

# BLAS thread count must be set before numpy (and, transitively, ADMMsolver)
# loads — see model_building.py. Respect an existing setting.
os.environ.setdefault('OMP_NUM_THREADS', '4')

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyModules.SDPLifting import _compute_m_plus_model, read_constr_from_lp
from pyModules.Graphs import read_graph_from_dimacs

DATA = os.path.join('..', 'data', 'StableSets')
RESULT_DIR = os.path.join('..', 'results', 'm_kk_test')

RELAXATIONS = ('nod_gamma', 'nod_theta', 'nod_alpha', 'edge', 'cov')

DIMACS_SUBSET = (
    'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4',
    'p_hat300-1', 'p_hat300-2', 'p_hat300-3',
    'keller4', 'sanr400_0.5', 'sanr400_0.7',
)
RANDOM_SIZES = ('150', '175', '200')

CSV_FIELDS = [
    'dataset', 'instance', 'relax', 'n', 'edges', 'graph_density',
    'lp_rows', 'cols', 'est_rows', 'est_nnz',
    'rows', 'nnz', 'A_density', 'avg_row_nnz', 'max_row_nnz',
    'build_time', 'solve_time', 'status', 'objval', 'lp_bound',
    'alpha', 'rel_gap', 'bar_iters', 'simplex_iters', 'threads',
]


def manifest():
    """Yield (dataset, instance) pairs of the experiment, from existing data."""
    lp_dir = os.path.join(DATA, 'smallDIMACS', 'lp')
    names = sorted({f.rsplit('_', 1)[0].replace('_nod', '')
                    for f in os.listdir(lp_dir) if f.endswith('.lp')})
    for name in names:
        yield 'smallDIMACS', name
    for name in DIMACS_SUBSET:
        yield 'DIMACS', name
    lp_dir = os.path.join(DATA, 'Random', 'lp')
    rnames = sorted({f.rsplit('_', 1)[0].replace('_nod', '')
                     for f in os.listdir(lp_dir) if f.endswith('.lp')
                     and f.split('_')[1] in RANDOM_SIZES})
    for name in rnames:
        yield 'Random', name


def refresh_coefficients(only=None):
    """Recompute theta/alpha caches and regenerate the NOD LPs from scratch.

    For each manifest instance (filtered by `only`), recomputes the nodal
    coefficients with refresh=True — overwriting only that instance's JSONs
    in data/.../coeff_theta|coeff_alpha — validates theta >= alpha, and
    rewrites the <name>_nod_theta.lp / <name>_nod_alpha.lp files the lifting
    reads (nod_gamma is degree-derived and needs no cache). Returns the set
    of refreshed instance names.
    """
    from pyModules.LinearFormulations import (NOD, compute_theta, compute_alpha,
                                              check_nodal_coefficients)
    from parameters import (ADMM_SOLVER_TOL, ADMM_SOLVER_MAX_ITER,
                            ADMM_SOLVER_TIMELIMIT, ADMM_SOLVER_DEBUG)
    # Same ADMM options as model_building.py
    admm_options = {
        'tolerance': ADMM_SOLVER_TOL,
        'max_iter':  ADMM_SOLVER_MAX_ITER,
        'timelimit': ADMM_SOLVER_TIMELIMIT,
        'debug':     ADMM_SOLVER_DEBUG,
    }
    refreshed = set()
    for dataset, name in manifest():
        if only and only not in name:
            continue
        stb_path = os.path.join(DATA, dataset, 'graphs', name + '.stb')
        if not os.path.exists(stb_path):
            print('MISSING %s' % stb_path)
            continue
        G = read_graph_from_dimacs(stb_path)
        print('refresh %s/%s (n=%d)' % (dataset, name, G.number_of_nodes()),
              flush=True)
        # compute_theta/compute_alpha write their JSON caches (and cliquer's
        # temp .stb) with plain open(); they don't create the directory
        theta_dir = os.path.join(DATA, dataset, 'coeff_theta')
        alpha_dir = os.path.join(DATA, dataset, 'coeff_alpha')
        os.makedirs(theta_dir, exist_ok=True)
        os.makedirs(alpha_dir, exist_ok=True)
        theta = compute_theta(G, name, refresh=True, use_patience=True,
                              options=admm_options, model_out_dir=theta_dir)
        alpha = compute_alpha(G, name, refresh=True,
                              model_out_dir=alpha_dir)
        check_nodal_coefficients(theta, alpha, name)
        lp_dir = os.path.join(DATA, dataset, 'lp')
        NOD(G, name + '_nod_theta', theta, lp_dir)
        NOD(G, name + '_nod_alpha', alpha, lp_dir)
        refreshed.add(name)
    return refreshed


def load_alpha(dataset):
    """{instance: alpha} from the dataset's aux_data csv (whitespace-separated)."""
    path = os.path.join(DATA, dataset,
                        'aux_data_%s.csv' % dataset.lower())
    table = {}
    if not os.path.exists(path):
        return table
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                table[parts[0]] = float(parts[1])
    return table


def estimate(lp_path):
    """(n, m, est_rows, est_nnz) of the kk+squares model, from LP supports.

    Supports of the homogenized forms: |a|+1 per LP row, 1 and 2 for the
    bounds. Product row (k,l) has at most s_k*s_l packed nonzeros, so
    sum_{k<=l} s_k s_l = ((sum s)^2 + sum s^2) / 2 bounds the matrix mass.
    """
    constrs = list(read_constr_from_lp(lp_path))
    n, lp_rows = len(constrs[0]), constrs[1:]
    s = np.array([len(c) + 1 for c, _ in lp_rows] + [1] * n + [2] * n,
                 dtype=np.float64)
    forms = len(s)
    est_rows = forms * (forms + 1) // 2
    est_nnz = (s.sum() ** 2 + (s ** 2).sum()) / 2
    return n, len(lp_rows), est_rows, int(est_nnz)


def solve_barrier(model, time_limit, threads, log_file, write_lp=None):
    import gurobipy as gp
    from gurobipy import GRB

    d = model.to_lp()
    lp = gp.Model('m_kk')
    lp.Params.OutputFlag = 0
    lp.Params.LogFile = log_file
    lp.Params.Method = 2          # barrier
    lp.Params.Crossover = 0       # bound only: no basic solution needed
    lp.Params.TimeLimit = time_limit
    lp.Params.Threads = threads

    y = lp.addMVar(d['num_cols'], lb=d['lb'], ub=d['ub'], name='y')
    lp.setObjective(d['obj'] @ y, GRB.MINIMIZE)
    for sense in ('<', '>', '='):
        mask = d['senses'] == sense
        if mask.any():
            lp.addMConstr(d['A'][mask], y, sense, d['rhs'][mask])

    if write_lp:
        lp.write(write_lp)

    t0 = time.time()
    lp.optimize()
    solve_time = time.time() - t0

    status = {GRB.OPTIMAL: 'optimal', GRB.TIME_LIMIT: 'timelimit',
              GRB.NUMERIC: 'numeric'}.get(lp.Status, 'status_%d' % lp.Status)
    try:
        objval = lp.ObjVal
    except gp.GurobiError:
        objval = float('nan')
    return dict(status=status, objval=objval, solve_time=solve_time,
                bar_iters=int(lp.BarIterCount),
                simplex_iters=int(lp.IterCount))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--time-limit', type=float, default=7200.,
                    help='Gurobi TimeLimit in seconds — solve phase only '
                         '(default 7200)')
    ap.add_argument('--threads', type=int, default=0,
                    help='Gurobi Threads (default 0 = all cores)')
    ap.add_argument('--rows-cap', type=float, default=2e7,
                    help='skip models with more estimated rows (default 2e7)')
    ap.add_argument('--nnz-cap', type=float, default=8e8,
                    help='skip models with more estimated matrix nonzeros '
                         '(default 8e8, ~10 GB peak)')
    ap.add_argument('--only', default=None,
                    help='process only instances whose name contains this')
    ap.add_argument('--relax', default=None, choices=RELAXATIONS,
                    help='process only this relaxation')
    ap.add_argument('--refresh-coeffs', action='store_true',
                    help='recompute theta/alpha coefficients from scratch for '
                         'the selected instances (overwriting their JSON '
                         'caches) and regenerate their nod_theta/nod_alpha '
                         'LP files before running')
    ap.add_argument('--dry-run', action='store_true',
                    help='print the size map (build/skip decisions), no solves')
    args = ap.parse_args()

    os.makedirs(RESULT_DIR, exist_ok=True)
    log_dir = os.path.join(RESULT_DIR, 'gurobi_logs')
    os.makedirs(log_dir, exist_ok=True)
    # lifted M(K,K) LPs are exported here; the source relaxation LPs stay
    # in data/StableSets/<dataset>/lp
    lp_out_dir = os.path.join(RESULT_DIR, 'lp')
    os.makedirs(lp_out_dir, exist_ok=True)
    csv_path = os.path.join(RESULT_DIR, 'results_m_kk.csv')

    # ------------------------------------------------------------------
    # PHASE 1: Coefficient refresh (only with --refresh-coeffs)
    # ------------------------------------------------------------------
    # Recompute theta/alpha caches from scratch and regenerate the
    # nod_theta/nod_alpha LP files for the selected instances
    refreshed = set()
    if args.refresh_coeffs:
        refreshed = refresh_coefficients(args.only)
        print('Refreshed coefficients + NOD LPs for %d instance(s)'
              % len(refreshed))

    # ------------------------------------------------------------------
    # PHASE 2: Resume state
    # ------------------------------------------------------------------
    # (instance, relax) pairs already in the results CSV are skipped, so
    # the run can be resumed after a crash or extended incrementally.
    # Rows solved before a coefficient refresh are stale: they are
    # reported but never deleted automatically.
    done = set()
    if os.path.exists(csv_path) and not args.dry_run:
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                done.add((row['instance'], row['relax']))
        print('Resuming: %d results already in %s' % (len(done), csv_path))
        stale = sorted((n, r) for n, r in done
                       if n in refreshed and r in ('nod_theta', 'nod_alpha'))
        if stale:
            print('WARNING: %d result row(s) predate the coefficient refresh '
                  'and will be skipped by resume; delete them from the CSV to '
                  're-solve: %s' % (len(stale), ', '.join(map(str, stale))))

    # graph-level alpha values (aux_data csv), for the rel_gap column
    alphas = {ds: load_alpha(ds) for ds in ('smallDIMACS', 'DIMACS', 'Random')}

    # ------------------------------------------------------------------
    # PHASE 3: Task gathering with size estimates
    # ------------------------------------------------------------------
    # One task per (instance, relaxation) LP in the manifest, with the
    # estimated lifted-model size; sorted cheapest first so results land
    # early. Oversize models are kept as tasks and recorded as SKIP rows.
    tasks = []
    for dataset, name in manifest():
        if args.only and args.only not in name:
            continue
        for relax in RELAXATIONS:
            if args.relax and relax != args.relax:
                continue
            if (name, relax) in done:
                continue
            lp_path = os.path.join(DATA, dataset, 'lp',
                                   '%s_%s.lp' % (name, relax))
            if not os.path.exists(lp_path):
                print('MISSING %s' % lp_path)
                continue
            n, m, est_rows, est_nnz = estimate(lp_path)
            tasks.append((est_nnz, dataset, name, relax, lp_path,
                          n, m, est_rows))
    tasks.sort()
    print('%d tasks (rows cap %.1e, nnz cap %.1e)'
          % (len(tasks), args.rows_cap, args.nnz_cap))

    # append to the results CSV (write the header only on first creation)
    new_file = not os.path.exists(csv_path)
    out = None
    if not args.dry_run:
        out = open(csv_path, 'a', newline='')
        writer = csv.DictWriter(out, fieldnames=CSV_FIELDS)
        if new_file:
            writer.writeheader()

    # ------------------------------------------------------------------
    # PHASE 4: Build + solve (or dry-run size map)
    # ------------------------------------------------------------------
    # For each task: build the M+(K,K) lifting of the LP in memory, solve
    # the M(K,K) LP with Gurobi barrier (no crossover), and append one
    # CSV row. With --dry-run, only print the build/skip map.
    for est_nnz, dataset, name, relax, lp_path, n, m, est_rows in tasks:
        oversize = est_rows > args.rows_cap or est_nnz > args.nnz_cap
        if args.dry_run:
            print('%-6s %-12s %-24s %-9s n=%-4d m=%-6s est_rows=%.2e est_nnz=%.2e'
                  % ('SKIP' if oversize else 'build', dataset, name, relax,
                     n, m, est_rows, est_nnz))
            continue

        rec = dict.fromkeys(CSV_FIELDS, '')
        alpha = alphas[dataset].get(name)
        rec.update(dataset=dataset, instance=name, relax=relax, n=n,
                   lp_rows=m, est_rows=est_rows, est_nnz=est_nnz,
                   alpha='' if alpha is None else alpha,
                   threads=args.threads)

        if oversize:
            rec['status'] = 'skip_size'
            print('SKIP  %-24s %-9s est_rows=%.2e est_nnz=%.2e'
                  % (name, relax, est_rows, est_nnz), flush=True)
        else:
            log_file = os.path.join(log_dir, '%s_%s.log' % (name, relax))
            edges = gd = ''
            t0 = time.time()
            model = _compute_m_plus_model(lp_path, lift_mode='kk',
                                          include_squares=True)
            build_time = time.time() - t0
            # graph density metric
            stb = os.path.join(DATA, dataset, 'graphs', name + '.stb')
            if os.path.exists(stb):
                G = read_graph_from_dimacs(stb)
                edges = G.number_of_edges()
                gd = 2. * edges / (n * (n - 1))
            res = solve_barrier(model, args.time_limit, args.threads,
                                log_file,
                                write_lp=os.path.join(lp_out_dir,
                                                      '%s_%s_m_kk.lp' % (name, relax)))

            row_nnz = np.diff(model.A.indptr)
            lp_bound = -res['objval'] if res['status'] == 'optimal' else ''
            rel_gap = ''
            if alpha and lp_bound != '':
                rel_gap = (lp_bound - alpha) / alpha
            rec.update(edges=edges, graph_density=gd,
                       cols=model.packed_dim, rows=model.num_rows,
                       nnz=model.A.nnz,
                       A_density=model.A.nnz / (model.num_rows * model.packed_dim),
                       avg_row_nnz='%.1f' % row_nnz.mean(),
                       max_row_nnz=int(row_nnz.max()),
                       build_time='%.2f' % build_time,
                       solve_time='%.2f' % res['solve_time'],
                       status=res['status'], objval=res['objval'],
                       lp_bound=lp_bound, rel_gap=rel_gap,
                       bar_iters=res['bar_iters'],
                       simplex_iters=res['simplex_iters'])
            print('%-8s %-24s %-9s rows=%-9d nnz=%.2e build=%ss solve=%ss bound=%s'
                  % (res['status'], name, relax, model.num_rows,
                     model.A.nnz, rec['build_time'], rec['solve_time'],
                     lp_bound), flush=True)
            del model

        writer.writerow(rec)
        out.flush()

    if out:
        out.close()
    print('Done. Results: %s' % csv_path)


if __name__ == '__main__':
    main()
