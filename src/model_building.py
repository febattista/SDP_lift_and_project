# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Main script: computes LP coefficients, builds LP formulations, and generates
# SDP models for all datasets specified in parameters.py.

import os, zipfile, itertools

# Import parameters first so ADMM_SOLVER_THREADS is available before
# LinearFormulations (and thus ADMMsolver) is imported.
from parameters import (
    datasets, results_datasets,
    MAKE_COEFFICIENTS, MAKE_LP_FORMULATION, MAKE_SDP_MODELS, MAKE_TH_PLUS_MODELS,
    DO_SPLIT, ADMM_SOLVER_THREADS, SDP_LIFTING_STEP,
    ADMM_SOLVER_TOL, ADMM_SOLVER_MAX_ITER, ADMM_SOLVER_TIMELIMIT, ADMM_SOLVER_DEBUG,
    USE_CLIQUE_COVER_FROM_LETCHFORD_ET_AL,
)

# Set thread count before importing ADMMsolver (transitively imported by LinearFormulations)
os.environ["OMP_NUM_THREADS"] = str(ADMM_SOLVER_THREADS)

from pyModules.Graphs import read_graph_from_dimacs
from pyModules.SDPLifting import m_plus_lifting, Theta_plus_SDP, lovasz_schrijver_filter
from pyModules.LinearFormulations import (
    FRAC, QSTABC, NOD, compute_gamma, compute_theta, compute_alpha,
)


def _banner(title):
    print("\n--- %-44s" % (title + ' ') + '-' * max(0, 46 - len(title)))


def add_files_to_zip(zip_filename, files):
    with zipfile.ZipFile(zip_filename, 'a', compression=zipfile.ZIP_DEFLATED) as zip_file:
        for file in files:
            print('    Deflating: %s' % os.path.basename(file))
            zip_file.write(file, arcname=os.path.basename(file))
            if os.path.exists(file):
                os.remove(file)


admm_options = {
    'tolerance': ADMM_SOLVER_TOL,
    'max_iter':  ADMM_SOLVER_MAX_ITER,
    'timelimit': ADMM_SOLVER_TIMELIMIT,
    'debug':     ADMM_SOLVER_DEBUG,
}

for d in datasets:
    print("\n" + "=" * 50)
    print("  Dataset: %s" % d)
    print("=" * 50)

    data_path    = datasets[d]
    results_path = results_datasets[d]

    graph_dir       = os.path.join(data_path,    'graphs')
    lp_dir          = os.path.join(data_path,    'lp')
    coeff_alpha_dir = os.path.join(data_path,    'coeff_alpha')
    coeff_theta_dir = os.path.join(data_path,    'coeff_theta')
    model_dir       = os.path.join(results_path, 'models')

    if not os.path.exists(data_path):
        print("  [SKIP] Dataset path not found: %s" % data_path)
        continue
    if not os.path.exists(graph_dir) or not any(
            f.endswith('.stb') for f in os.listdir(graph_dir)):
        print("  [SKIP] No .stb graphs found in %s" % graph_dir)
        continue

    os.makedirs(model_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # PHASE 1: Coefficient computation (theta, alpha)
    # ------------------------------------------------------------------
    _banner("Phase 1: Coefficient computation")
    if not MAKE_COEFFICIENTS:
        print("  [SKIP] MAKE_COEFFICIENTS = False")
    else:
        for path in [coeff_alpha_dir, coeff_theta_dir]:
            os.makedirs(path, exist_ok=True)

        with os.scandir(graph_dir) as inst_it:
            for instance in sorted(inst_it, key=lambda e: e.name):
                if instance.name.lower().endswith('.stb'):
                    graphname = os.path.splitext(instance.name)[0]
                    print("  > %s" % graphname)
                    G = read_graph_from_dimacs(os.path.join(graph_dir, instance.name))
                    compute_theta(G, graphname, model_out_dir=coeff_theta_dir,
                                  use_patience=True, options=admm_options)
                    compute_alpha(G, graphname, model_out_dir=coeff_alpha_dir)

    # ------------------------------------------------------------------
    # PHASE 2: LP formulation
    # ------------------------------------------------------------------
    _banner("Phase 2: LP formulation")
    if not MAKE_LP_FORMULATION:
        print("  [SKIP] MAKE_LP_FORMULATION = False")
    else:
        for path in [lp_dir, coeff_alpha_dir, coeff_theta_dir]:
            os.makedirs(path, exist_ok=True)

        with os.scandir(graph_dir) as inst_it:
            for instance in sorted(inst_it, key=lambda e: e.name):
                if instance.name.lower().endswith('.stb'):
                    graphname = os.path.splitext(instance.name)[0]
                    print("  > %s" % graphname)
                    G = read_graph_from_dimacs(os.path.join(graph_dir, instance.name))

                    gamma = compute_gamma(G)
                    NOD(G, graphname + '_nod_gamma', gamma, lp_dir)

                    theta = compute_theta(G, graphname,
                                         model_out_dir=coeff_theta_dir,
                                         use_patience=True, options=admm_options)
                    NOD(G, graphname + '_nod_theta', theta, lp_dir)

                    alpha = compute_alpha(G, graphname, model_out_dir=coeff_alpha_dir)
                    NOD(G, graphname + '_nod_alpha', alpha, lp_dir)

                    QSTABC(G, graphname + '_cov', dir_path=lp_dir,
                           use_letchford=USE_CLIQUE_COVER_FROM_LETCHFORD_ET_AL)
                    FRAC(G, graphname + '_edge', dir_path=lp_dir)

    # ------------------------------------------------------------------
    # PHASE 3: SDP model generation
    # ------------------------------------------------------------------
    _banner("Phase 3: SDP model generation")
    if not MAKE_SDP_MODELS:
        print("  [SKIP] MAKE_SDP_MODELS = False")
    else:
        out_path = os.path.join(model_dir, 'models_%s.zip' % d.lower())
        if os.path.exists(out_path):
            print("  [SKIP] Archive already exists: %s" % out_path)
        else:
            with zipfile.ZipFile(out_path, 'a'):
                pass

            if MAKE_TH_PLUS_MODELS:
                print("  > Theta+ SDPs ...")
                with os.scandir(graph_dir) as inst_it:
                    for instance in sorted(inst_it, key=lambda e: e.name):
                        if instance.name.lower().endswith('.stb'):
                            filename = os.path.splitext(instance.name)[0] + '_th+'
                            G = read_graph_from_dimacs(os.path.join(graph_dir, instance.name))
                            Theta_plus_SDP(G, filename, model_out_dir=model_dir,
                                           model_out='sdpnal', debug=True)
                            add_files_to_zip(out_path,
                                             [os.path.join(model_dir, filename + '.mat')])
            else:
                print("  [SKIP] Theta+ models (MAKE_TH_PLUS_MODELS = False)")

            if not os.path.exists(lp_dir) or not any(
                    f.endswith('.lp') for f in os.listdir(lp_dir)):
                print("  [SKIP] No .lp files found in %s — run Phase 2 first" % lp_dir)
            else:
                with os.scandir(lp_dir) as inst_it:
                    for instance in sorted(inst_it, key=lambda e: e.name):
                        if instance.name.lower().endswith('.lp'):
                            filename = os.path.splitext(instance.name)[0]
                            m_plus_lifting(
                                instance.name,
                                model_out_dir=model_dir,
                                work_dir=lp_dir,
                                step=SDP_LIFTING_STEP,
                                skip_func=(lovasz_schrijver_filter
                                           if 'edge' in instance.name else None),
                                do_split=DO_SPLIT,
                            )
                            files_to_zip = [
                                os.path.join(model_dir, ''.join(s))
                                for s in itertools.product(
                                    [filename],
                                    ['_1.mat', '_2.mat', '_3.mat', '_4.mat', '_5.mat']
                                    if DO_SPLIT else ['.mat'],
                                )
                            ]
                            add_files_to_zip(out_path, files_to_zip)

    print("\n  Done: %s" % d)
