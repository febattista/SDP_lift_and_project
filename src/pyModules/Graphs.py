# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Graph generation, I/O utilities, and clique cover algorithms.

import networkx as nx
import scipy.io as io
import numpy as np
import re, itertools


def web(p, q):
    """Return a (p, q)-web graph. Assumes p > 2q + 1."""
    G = nx.Graph()
    edges = []
    for i in range(p):
        edges += [(i % p, j % p) for j in range(i + q , p + i - q +1)]
    G.add_edges_from(edges)
    return G


def antiweb(p, q):
    """Return a (p, q)-antiweb graph. Assumes p > 2q + 1."""
    G = nx.complement(web(p, q))
    return G


def wheel(n):
    """Return a n-wheel graph."""
    return nx.wheel_graph(n)


def cycle(n):
    """Return a chordless n-cycle graph."""
    return nx.cycle_graph(n)


def random_graph(n, p, seed=1234):
    """Return a random graph G(n, p) in the Erdos-Renyi model."""
    return nx.erdos_renyi_graph(n, p, seed=seed)


def random_perfect_graph(n, p, seed=1234):
    """Return a random perfect graph as line graph of a random bipartite graph.

    Parameters
    ----------
    n : int
        Expected number of nodes in the line graph (number of edges in
        the underlying bipartite graph).
    p : float
        Edge density of the random bipartite graph.
    """
    return nx.convert_node_labels_to_integers(
        nx.line_graph(nx.bipartite.random_graph(n, n, p, seed=seed))
    )


def greedy_clique_cover(G):
    """Return a greedy clique cover as a list of node lists.

    At each step the highest-degree node (in the remaining edge subgraph)
    seeds a new clique, which is grown greedily through its neighbours.
    Covered edges are removed from the working copy until no edges remain.
    """
    Gcopy = G.copy()
    cliques = dict()
    i = 0
    while len(Gcopy.edges()):
        node = sorted(Gcopy.nodes(), key=Gcopy.degree, reverse=True)[0]
        cliques[i] = set([node])
        frontier = set(Gcopy[node])
        done = False
        while not done:
            if frontier:
                node = sorted(frontier, key=Gcopy.degree, reverse=True)[0]
                cliques[i].add(node)
                frontier = frontier.intersection(set(Gcopy[node]))
            else:
                done = True
        Gcopy.remove_edges_from(G.subgraph(cliques[i]).edges)
        i += 1
    return [list(clique) for i, clique in cliques.items()]


def find_maximal_clique(G, u, v, nodeweight):
    """Yield nodes of a maximal clique containing edge (u, v).

    Extends the clique greedily by picking the highest-weight candidate
    that is adjacent to all current clique members.
    """
    yield u
    yield v

    candidates = set(nx.common_neighbors(G, u, v))

    while candidates:
        inclique = max(candidates, key=lambda x: nodeweight[x])
        yield inclique
        candidates = candidates.intersection(G.neighbors(inclique))


def greedy_clique_cover_letchford_et_al(inG):
    """Return a clique cover using the heuristic of Letchford et al.

    Iteratively selects the uncovered edge incident to the highest-weight
    node, builds a maximal clique through that edge, marks all its edges
    covered, and updates node weights accordingly.
    """
    G = inG.copy()

    cliques = set()

    nodeweight = {i: G.degree(i) for i in G.nodes()}

    nx.set_edge_attributes(G, False, name='covered')

    find_clique = True

    while find_clique:

        root_u = max(nodeweight, key=lambda key: nodeweight[key])

        if nodeweight[root_u] > 0:
            max_weight = -1
            root_v = None

            for u, v, data in filter(lambda x: not x[2]['covered'], G.edges(root_u, data=True)):
                if nodeweight[v] > max_weight:
                    max_weight = nodeweight[v]
                    root_v = v

            if root_v is not None:
                clique = list(find_maximal_clique(G, root_u, root_v, nodeweight))

                cliques.add(tuple(sorted(clique)))

                for u, v in itertools.combinations(clique, 2):
                    if not G[u][v]['covered']:
                        G[u][v]['covered'] = True
                        nodeweight[u] -= 1.0
                        nodeweight[v] -= 1.0
            else:
                find_clique = False
        else:
            find_clique = False

    return [list(clique) for clique in cliques]


def G_hat():
    G = nx.Graph()
    G.add_nodes_from(range(0, 8))
    edges = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 5), (1, 6), (2, 3),
             (2, 6), (2, 7), (3, 4), (3, 7), (4, 5), (4, 6), (4, 7), (5, 6),
             (5, 7), (6, 7)]
    G.add_edges_from(edges)
    return G


def G_LT():
    G = nx.wheel_graph(6)
    G.remove_edge(0, 1)
    G.remove_edge(0, 5)
    return G


def G_EMN():
    G = nx.wheel_graph(6)
    G.remove_edge(0, 1)
    return G


def write_graph_to_dimacs(G, path):
    """Write the NetworkX graph G to a DIMACS-format file at path."""
    with open(path, "w") as f:
        f.write("p edge {} {}\n".format(G.number_of_nodes(), G.number_of_edges()))
        for u, v in G.edges():
            f.write("e {} {}\n".format(u + 1, v + 1))


def read_graph_from_dimacs(path):
    """Read a DIMACS-format file and return a NetworkX graph."""
    try:
        with open(path, 'r') as f:
            edges = []
            for line in f:
                if line[0] in ['c', '', '\n']:
                    continue
                if line[0] in ['p']:
                    n, m = re.split(r'\s+', line)[2:4]
                    continue
                edge = re.split(r'\s+', line)
                edge = int(edge[1]) - 1, int(edge[2]) - 1
                edges.append(edge)
    except IOError:
        raise

    G = nx.Graph()
    G.add_edges_from(edges)
    return G


def read_graph_from_mat(path, var, type='adjacency'):
    """Read a graph stored in a .mat file."""
    try:
        G = nx.Graph()
        data = io.loadmat(path)
        if 'adjacency' in type:
            G = nx.from_numpy_matrix(np.array(data[var]))
        elif 'edgelist' in type:
            G.add_edges_from(np.array(data[var]))
    except IOError:
        raise

    return nx.convert_node_labels_to_integers(G)


def read_alpha_from_lp(filename):
    """Read nodal coefficients from a Nodal LP file.

    Returns a dict mapping node index (0-based) to coefficient value.
    """
    alphas = {}
    pattern = re.compile(r"\d+\sx\d+")
    with open(filename) as f:
        for line in f:
            match = re.findall(pattern, line)
            if match:
                match = match[0].split(' ')
                alphas[int(match[1][1:]) - 1] = int(match[0])
    return alphas
