# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# LP relaxation formulations for the Maximum Stable Set Problem using Gurobi.
#
# Indexing conventions
# --------------------
# LP variables are named x[1]..x[n]: graph node i (0-based, see Graphs.py)
# becomes variable x[i+1]. The coefficient dicts returned by compute_gamma /
# compute_theta / compute_alpha — and their JSON caches in coeff_theta/ and
# coeff_alpha/ — are keyed by the 0-based node.

import os, json, subprocess, time
import gurobipy as gb
from tqdm import tqdm

from pyModules.Graphs import write_graph_to_dimacs
from pyModules.ADMMsolver import ADMMsolver, Problem
from pyModules.SDPLifting import Theta_SDP
import networkx as nx
import numpy as np


def _add_obj_label(path):
    """Post-process a Gurobi-written LP file to add the 'obj:' label."""
    with open(path) as f:
        lines = f.readlines()

    found = False
    new_lines = []
    for line in lines:
        if not found and 'x[' in line:
            new_lines.append("obj: " + line)
            found = True
        else:
            new_lines.append(line)

    with open(path, 'w') as f:
        f.writelines(new_lines)


def FRAC(G, model_name, dir_path='', write_lp=True):
    """Build the edge LP relaxation (FRAC) for the stable set on G.

    Each edge (i, j) contributes the constraint x_i + x_j <= 1.
    """
    path = os.path.join(dir_path, model_name + '.lp')
    frac = gb.Model()
    # 0-based graph node i -> 1-based LP variable x[i+1] (module convention)
    x = frac.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    frac.setObjective(x.sum(), gb.GRB.MAXIMIZE)
    k = 0
    for i, j in G.edges():
        k += 1
        frac.addConstr(x[i+1] + x[j+1] <= 1, name=('edge%d' % k))
    frac.update()

    if write_lp:
        print('Saving LP at: %s' % path)
        frac.write(path)
        _add_obj_label(path)

    return frac


def QSTAB(G, model_name, dir_path='', write_lp=True):
    """Build the clique LP relaxation (QSTAB) using all maximal cliques."""
    cliques = nx.find_cliques(G)
    path = os.path.join(dir_path, model_name + '.lp')
    cov = gb.Model()
    x = cov.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    cov.setObjective(x.sum(), gb.GRB.MAXIMIZE)
    for i, K in enumerate(cliques):
        K_ = [k + 1 for k in K]
        if len(K) > 1:
            cov.addConstr(x.sum(K_) <= 1, name='clq' + str(i + 1))
    cov.update()

    if write_lp:
        print('Saving LP at: %s' % path)
        cov.write(path)
        _add_obj_label(path)

    return cov


def QSTABC(G, model_name, dir_path='', write_lp=True, use_letchford=True):
    """Build the clique cover LP relaxation (QSTABC) using a greedy cover.

    Parameters
    ----------
    use_letchford : bool
        If True (default), use greedy_clique_cover_letchford_et_al (tighter
        cover, slower). If False, use greedy_clique_cover (faster, looser).
    """
    from pyModules.Graphs import greedy_clique_cover, greedy_clique_cover_letchford_et_al
    cliques = greedy_clique_cover_letchford_et_al(G) if use_letchford else greedy_clique_cover(G)
    path = os.path.join(dir_path, model_name + '.lp')
    cov = gb.Model()
    x = cov.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    cov.setObjective(x.sum(), gb.GRB.MAXIMIZE)
    for i, K in enumerate(cliques):
        K_ = [k + 1 for k in K]
        if len(K) > 1:
            cov.addConstr(x.sum(K_) <= 1, name='clq' + str(i + 1))
    cov.update()

    if write_lp:
        print('Saving LP at: %s' % path)
        cov.write(path)
        _add_obj_label(path)

    return cov


def NOD(G, model_name, coeff, dir_path=''):
    """Build the nodal LP relaxation for the stable set on G.

    For each node i: sum_{j in N(i)} x_j + coeff[i] * x_i <= coeff[i].
    The coefficient dict (gamma, theta, or alpha) controls the tightness.
    """
    path = os.path.join(dir_path, model_name + '.lp')
    nod = gb.Model()
    x = nod.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    nod.setObjective(x.sum(), gb.GRB.MAXIMIZE)

    for i in G.nodes():
        neighbors = [j+1 for j in G.neighbors(i)]
        # coeff is keyed by the 0-based node i; the LP variable is x[i+1]
        nod.addConstr(x.sum(neighbors) + coeff[i]*x[i + 1] <= coeff[i], name='nod' + str(i + 1))

    nod.update()
    print('Saving LP at: %s' % path)
    nod.write(path)
    _add_obj_label(path)

    return nod


def compute_gamma(G):
    """Return the degree of each node as a dict {node: degree}."""
    return dict(G.degree())


def compute_theta(G, graphname, model_out_dir='', use_patience=False, options=None,
                  refresh=False):
    """Compute the theta coefficient for each node of G via ADMM.

    For each node v, solves the Theta SDP on the subgraph induced by N(v).
    Results are cached as a JSON file in model_out_dir; refresh=True ignores
    an existing cache and overwrites it.

    Returns a dict {node: theta_value}.
    """
    if not options:
        options = {
            'tolerance': 1e-4,
            'max_iter': 1000000,
            'timelimit': 3600,
            'print_it': 1000,
            'debug': False,
        }

    json_path = os.path.join(model_out_dir, "%s_theta.json" % graphname)
    if not refresh and os.path.exists(json_path):
        with open(json_path) as f:
            j = json.load(f)
            return {int(i): j[i][0] for i in j}

    theta = {}
    nodes = list(G.nodes())
    with tqdm(total=len(nodes), desc='    theta', unit='node') as pbar:
        for j in nodes:
            H = G.subgraph(G.neighbors(j))
            H = nx.convert_node_labels_to_integers(H)
            adal = Theta_SDP(H).to_adal()
            P = Problem(adal['A'], adal['b'], adal['C'], c=0)
            # The safe dual bound needs norm_bound >= ||X*||_F. For the
            # homogenized Theta_SDP (dim n+1), ||X*||_F <= trace(X*) =
            # 1 + theta(H) <= n + 1
            opt_sol, _, _, secs = ADMMsolver.solve(P, norm_bound=len(H.nodes()) + 1, options=options, use_patience=use_patience)
            theta[j] = [np.floor(-opt_sol.safe_dual), secs]
            pbar.set_postfix(node=j, theta='%.2f' % theta[j][0])
            pbar.update(1)

    with open(json_path, 'w') as f:
        json.dump(theta, f)

    return {i: theta[i][0] for i in theta}


_CLIQUER = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 '..', '..', 'ThirdParty', 'cliquer-1.21', 'cl'))


def cliquer_max_clique(filepath):
    """Run cliquer on a DIMACS file and return (alpha, cpu_time)."""
    command = [_CLIQUER, filepath, '-q', '-q']
    start = time.time()
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    cltime = time.time() - start
    alpha = float(result.stdout.split(",")[0].split("=")[1])
    return alpha, cltime


def compute_alpha(G, graphname, model_out_dir='', refresh=False):
    """Compute the alpha (clique number) coefficient for each node of G.

    For each node v, finds the maximum clique in the complement of N(v)
    using cliquer. Results are cached as a JSON file in model_out_dir;
    refresh=True ignores an existing cache and overwrites it.

    Returns a dict {node: alpha_value}.
    """
    json_path = os.path.join(model_out_dir, "%s_alpha.json" % graphname)
    if not refresh and os.path.exists(json_path):
        with open(json_path) as f:
            j = json.load(f)
            return {int(i): j[i][0] for i in j}

    tmp_path = os.path.join(model_out_dir, 'subgraph.stb')
    alpha = {}
    nodes = list(G.nodes())
    with tqdm(total=len(nodes), desc='    alpha', unit='node') as pbar:
        for j in nodes:
            H = G.subgraph(G.neighbors(j))
            H = nx.convert_node_labels_to_integers(H)
            write_graph_to_dimacs(nx.complement(H), tmp_path)
            a, cltime = cliquer_max_clique(tmp_path)
            alpha[j] = [a, cltime]
            pbar.set_postfix(node=j, alpha='%.2f' % alpha[j][0])
            pbar.update(1)

    if os.path.exists(tmp_path):
        os.remove(tmp_path)

    with open(json_path, 'w') as f:
        json.dump(alpha, f)

    return {i: alpha[i][0] for i in alpha}


def check_nodal_coefficients(theta, alpha, graphname=''):
    """Raise ValueError if any node has theta < alpha.

    theta(N(v)) >= alpha(N(v)) holds for every graph, so a violation proves
    at least one coefficient is wrong and the resulting NOD_theta relaxation
    would be invalid.
    """
    bad = {v: (theta[v], alpha[v]) for v in theta
           if v in alpha and theta[v] < alpha[v]}
    if bad:
        detail = '; '.join('node %s: theta=%g < alpha=%g' % (v, t, a)
                           for v, (t, a) in sorted(bad.items())[:10])
        raise ValueError(
            '%s: %d invalid nodal coefficient(s) (theta < alpha): %s%s'
            % (graphname or 'coefficients', len(bad), detail,
               ' ...' if len(bad) > 10 else ''))
