#!/usr/bin/env python3
# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
"""Solve a Theta, Theta+, or M+ SDP using MOSEK.

Usage
-----
  # Theta+ (default) or Theta from a graph file
  python solve_mosek.py path/to/graph.stb
  python solve_mosek.py path/to/graph.stb --formulation theta

  # M+ lift-and-project from an LP relaxation file
  python solve_mosek.py path/to/formulation.lp
"""

import sys
import os
import argparse
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import mosek
except ImportError:
    sys.exit("MOSEK is not installed.  Install it with:  pip install mosek")

import numpy as np

from pyModules.Graphs import read_graph_from_dimacs
from pyModules.SDPLifting import Theta_SDP, Theta_plus_SDP, _compute_m_plus_model


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def build_model(path, formulation):
    ext = os.path.splitext(path)[1].lower()

    if ext == '.stb':
        G = read_graph_from_dimacs(path)
        print('Graph       : %s' % os.path.basename(path))
        print('Nodes / edges: %d / %d' % (len(G.nodes()), len(G.edges())))
        if formulation == 'theta':
            return Theta_SDP(G)
        else:
            return Theta_plus_SDP(G)

    elif ext == '.lp':
        if formulation == 'theta' or formulation == 'theta_plus':
            sys.exit("Theta / Theta+ require a .stb graph file, not an LP file.")
        print('LP file     : %s' % os.path.basename(path))
        return _compute_m_plus_model(path)

    else:
        sys.exit("Unsupported file extension %r.  Expected .stb or .lp." % ext)


# ---------------------------------------------------------------------------
# MOSEK solve
# ---------------------------------------------------------------------------

def solve(model, verbose=False):
    d = model.to_mosek()

    with mosek.Task() as task:
        if not verbose:
            task.putintparam(mosek.iparam.log, 0)

        task.appendcons(model.num_rows)
        task.appendbarvars([d['dim']])

        # Objective matrix C (lower-triangular triplets)
        n_c = len(d['C_val'])
        if n_c > 0:
            task.putbarcblocktriplet(
                [0] * n_c,
                d['C_rows'].tolist(),
                d['C_cols'].tolist(),
                d['C_val'].tolist(),
            )

        # Constraint matrices (lower-triangular triplets)
        n_a = len(d['A_val'])
        if n_a > 0:
            task.putbarablocktriplet(
                d['A_cons'].tolist(),
                [0] * n_a,
                d['A_rows'].tolist(),
                d['A_cols'].tolist(),
                d['A_val'].tolist(),
            )

        blc = d['blc'].tolist()
        buc = d['buc'].tolist()
        bk  = [mosek.boundkey.fx if s == '=' else
               mosek.boundkey.lo if s == '>' else
               mosek.boundkey.up for s in d['row_senses']]
        task.putconboundlist(list(range(model.num_rows)), bk, blc, buc)

        # Objective sense
        task.putobjsense(
            mosek.objsense.minimize if d['obj_sense'] > 0
            else mosek.objsense.maximize
        )

        t0 = time.time()
        task.optimize()
        elapsed = time.time() - t0

        sol_sta = task.getsolsta(mosek.soltype.itr)
        obj     = task.getprimalobj(mosek.soltype.itr)

    return obj, sol_sta, elapsed


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Solve a Theta, Theta+, or M+ SDP using MOSEK.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'examples:\n'
            '  python solve_mosek.py graph.stb\n'
            '  python solve_mosek.py graph.stb --formulation theta\n'
            '  python solve_mosek.py edge.lp\n'
        ),
    )
    parser.add_argument('file',
                        help='Path to a .stb graph file or .lp LP file.')
    parser.add_argument('--formulation', '-f',
                        choices=['theta', 'theta_plus', 'm_plus'],
                        default=None,
                        help=('SDP formulation.  Defaults to theta_plus for .stb '
                              'files and m_plus for .lp files.'))
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print MOSEK solver log.')
    args = parser.parse_args()

    path = os.path.abspath(args.file)
    if not os.path.exists(path):
        sys.exit('File not found: %s' % path)

    ext = os.path.splitext(path)[1].lower()
    formulation = args.formulation or ('m_plus' if ext == '.lp' else 'theta_plus')

    print('Formulation : %s' % formulation)

    t0 = time.time()
    model = build_model(path, formulation)
    t_build = time.time() - t0

    eq_rows   = int((model.row_senses == '=').sum())
    ineq_rows = int((model.row_senses == '>').sum())
    print('PSD dim     : %d' % model.dim)
    print('Constraints : %d  (eq=%d  ineq=%d)' % (model.num_rows, eq_rows, ineq_rows))
    print('Build time  : %.2fs' % t_build)
    print('Solving ...')

    obj, sol_sta, t_solve = solve(model, verbose=args.verbose)

    # C has -1 on diagonal entries, so <C,X> = -sum(X_ii) and the SDP bound
    # on the stability number is -obj.
    sdp_bound = -obj

    print('Status      : %s' % sol_sta)
    print('Objective   : %.6f' % obj)
    print('SDP bound   : %.6f' % sdp_bound)
    print('Solve time  : %.2fs' % t_solve)


if __name__ == '__main__':
    main()
