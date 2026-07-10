# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
"""Numerical experiment: solve the M(K,K) LP relaxation with Gurobi.

Builds each model in memory from the existing LP files (no dataset
duplication, no .mat serialization) with lift_mode='kk' + squares, and
solves the PSD-free LP via SDPModel.to_lp() using barrier without
crossover (see docs/plans/kk-lp-gurobi-solving.md for the measurements
behind that choice). The 7200 s time limit is a Gurobi parameter, so it
bounds the solve phase only; build time is recorded separately.

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


def solve_barrier(model, time_limit, threads, log_file):
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

    y = lp.addMVar(d['num_cols'], lb=d['lb'], ub=d['ub'])
    lp.setObjective(d['obj'] @ y, GRB.MINIMIZE)
    for sense in ('<', '>', '='):
        mask = d['senses'] == sense
        if mask.any():
            lp.addMConstr(d['A'][mask], y, sense, d['rhs'][mask])

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
    ap.add_argument('--dry-run', action='store_true',
                    help='print the size map (build/skip decisions), no solves')
    args = ap.parse_args()

    os.makedirs(RESULT_DIR, exist_ok=True)
    log_dir = os.path.join(RESULT_DIR, 'gurobi_logs')
    os.makedirs(log_dir, exist_ok=True)
    csv_path = os.path.join(RESULT_DIR, 'results_m_kk.csv')

    done = set()
    if os.path.exists(csv_path) and not args.dry_run:
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                done.add((row['instance'], row['relax']))
        print('Resuming: %d results already in %s' % (len(done), csv_path))

    alphas = {ds: load_alpha(ds) for ds in ('smallDIMACS', 'DIMACS', 'Random')}

    # gather tasks with size estimates, cheapest first so results land early
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

    new_file = not os.path.exists(csv_path)
    out = None
    if not args.dry_run:
        out = open(csv_path, 'a', newline='')
        writer = csv.DictWriter(out, fieldnames=CSV_FIELDS)
        if new_file:
            writer.writeheader()

    for est_nnz, dataset, name, relax, lp_path, n, m, est_rows in tasks:
        oversize = est_rows > args.rows_cap or est_nnz > args.nnz_cap
        if args.dry_run:
            print('%-6s %-12s %-24s %-9s n=%-4d m=%-6d est_rows=%.2e est_nnz=%.2e'
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
            t0 = time.time()
            model = _compute_m_plus_model(lp_path, lift_mode='kk',
                                          include_squares=True)
            build_time = time.time() - t0

            row_nnz = np.diff(model.A.indptr)
            # graph density metric
            stb = os.path.join(DATA, dataset, 'graphs', name + '.stb')
            edges = gd = ''
            if os.path.exists(stb):
                G = read_graph_from_dimacs(stb)
                edges = G.number_of_edges()
                gd = 2. * edges / (n * (n - 1))

            res = solve_barrier(model, args.time_limit, args.threads,
                                os.path.join(log_dir, '%s_%s.log' % (name, relax)))
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
