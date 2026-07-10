#!/usr/bin/env python3
# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
"""Solve the M-operator LP relaxation (no PSD constraint) using Gurobi.

Builds the M+ lifting from an LP relaxation file, drops the PSD requirement
via SDPModel.to_lp(), and solves the resulting LP over the packed entries
of the matrix variable.  This is the only script with a Gurobi dependency
for LP solving; SDPModel itself stays solver-free.

Usage
-----
  # Classic Lovász-Schrijver pairs (LP row x bound, bound x bound)
  python solve_gurobi.py path/to/formulation.lp --mode ls

  # Full M(K,K): all distinct pairs, squares included by default
  # (without squares the LP is weaker than M(K,K) proper)
  python solve_gurobi.py path/to/formulation.lp
"""

import sys
import os
import argparse
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import gurobipy as gp
    from gurobipy import GRB
except ImportError:
    sys.exit("Gurobi is not installed.  Install it with:  pip install gurobipy")

from pyModules.SDPLifting import _compute_m_plus_model


def solve(model, verbose=False, write_lp=None):
    d = model.to_lp()

    lp = gp.Model('m_operator')
    if not verbose:
        lp.Params.OutputFlag = 0

    y = lp.addMVar(d['num_cols'], lb=d['lb'], ub=d['ub'], name='y')
    lp.setObjective(d['obj'] @ y + d['obj_offset'],
                    GRB.MINIMIZE if d['obj_sense'] > 0 else GRB.MAXIMIZE)

    A, rhs, senses = d['A'], d['rhs'], d['senses']
    for sense in ('<', '>', '='):
        mask = senses == sense
        if mask.any():
            lp.addMConstr(A[mask], y, sense, rhs[mask])

    if write_lp:
        lp.write(write_lp)

    t0 = time.time()
    lp.optimize()
    elapsed = time.time() - t0

    if lp.Status != GRB.OPTIMAL:
        sys.exit('Gurobi did not reach optimality (status %d)' % lp.Status)
    return lp.ObjVal, lp.Status, elapsed


def main():
    parser = argparse.ArgumentParser(
        description='Solve the M-operator LP relaxation using Gurobi.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'examples:\n'
            '  python solve_gurobi.py edge.lp\n'
            '  python solve_gurobi.py edge.lp --mode ls\n'
            '  python solve_gurobi.py edge.lp --write-lp m_operator.lp\n'
        ),
    )
    parser.add_argument('file', help='Path to a .lp LP relaxation file.')
    parser.add_argument('--mode', choices=['ls', 'kk'], default='kk',
                        help="Lifting pair selection (default: kk).")
    parser.add_argument('--no-squares', action='store_true',
                        help="kk mode: drop the squared-constraint products "
                             "(the LP becomes weaker than M(K,K) proper).")
    parser.add_argument('--write-lp', metavar='PATH', default=None,
                        help='Write the assembled model via Gurobi '
                             '(.lp or .mps, by extension).')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print Gurobi solver log.')
    args = parser.parse_args()

    path = os.path.abspath(args.file)
    if not os.path.exists(path):
        sys.exit('File not found: %s' % path)
    if os.path.splitext(path)[1].lower() != '.lp':
        sys.exit('Expected a .lp file.')

    print('LP file     : %s' % os.path.basename(path))
    print('Lift mode   : %s' % args.mode)

    t0 = time.time()
    model = _compute_m_plus_model(
        path, lift_mode=args.mode,
        include_squares=(args.mode == 'kk' and not args.no_squares))
    t_build = time.time() - t0

    eq_rows   = int((model.row_senses == '=').sum())
    ineq_rows = int((model.row_senses != '=').sum())
    print('LP columns  : %d' % model.packed_dim)
    print('Constraints : %d  (eq=%d  ineq=%d)' % (model.num_rows, eq_rows, ineq_rows))
    print('Build time  : %.2fs' % t_build)
    print('Solving ...')

    obj, status, t_solve = solve(model, verbose=args.verbose,
                                 write_lp=args.write_lp)

    # C has -1 on diagonal entries, so the objective is -sum(X_ii) and the
    # LP bound on the stability number is -obj (same as solve_mosek.py).
    print('Objective   : %.6f' % obj)
    print('LP bound    : %.6f' % -obj)
    print('Solve time  : %.2fs' % t_solve)


if __name__ == '__main__':
    main()
