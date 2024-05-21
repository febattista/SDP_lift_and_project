import networkx as nx
import scipy.io as io
import numpy as np
import re

'''
Return a (p, q)-web graph. Assumes p > 2q + 1
'''
def web(p, q):
    G = nx.Graph()
    edges = []
    for i in range(p):
        edges += [(i % p, j % p) for j in range(i + q , p + i - q +1)]
    G.add_edges_from(edges)
    return G


'''
Return a (p, q)-antiweb graph. Assumes p > 2q + 1
'''
def antiweb(p, q):
    G = nx.complement(web(p, q))
    return G

'''
Return a n-wheel graph. 
'''
def wheel(n):
    return nx.wheel_graph(n)


'''
Return a chordless n-cycle graph. 
'''
def cycle(n):
    return nx.cycle_graph(n)


'''
Return a random graph G(n, p) in the Erdos-Renyi model. 
'''
def random_graph(n, p, seed=1234):
    return nx.erdos_renyi_graph(n, p, seed=seed)


'''
Return a "random" perfect graph. 
'''
def random_perfect_graph(n, p, seed=1234):
    # Returns random perfect graph obtained as line graph of a random bipartite graph
    # n: the expected number of nodes in the line graph (i.e. number of edges in the bipartite graph)
    # p: density of the random bipartite graph
    return nx.convert_node_labels_to_integers(nx.line_graph(nx.bipartite.random_graph(n, n, p, seed=seed)))


'''
Given a graph G, returns a greedy clique cover as collection of lists of nodes. 
'''
def greedy_clique_cover(G):
    Gcopy = G.copy()
    cliques = dict()
    frontier = set()
    i = 0
    while len(Gcopy.edges()):
        node = sorted(G.nodes(), key=Gcopy.degree, reverse=True)[0]
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
    # print(cliques)
    return [list(clique) for i, clique in cliques.items()]


def G_hat():
    G_hat = nx.Graph()
    G_hat.add_nodes_from(range(0,8))
    edges = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 5), (1, 6), (2, 3), (2, 6), (2, 7), (3, 4), (3, 7), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)]
    G_hat.add_edges_from(edges)
    return G_hat


def G_LT():
    G_LT = nx.wheel_graph(6)
    G_LT.remove_edge(0, 1)
    G_LT.remove_edge(0, 5)
    return G_LT
    

def G_EMN():
    G_EMN = nx.wheel_graph(6)
    G_EMN.remove_edge(0, 1)
    return G_EMN


'''
Write the Networkx graph G into DIMACS format file
'''
def write_graph_to_dimacs(G, path):
    with open(path, "w") as f:
        # write the header
        f.write("p edge {} {}\n".format(G.number_of_nodes(), G.number_of_edges()))
        # now write all edges
        for u, v in G.edges():
            f.write("e {} {}\n".format(u + 1, v + 1))

'''
Read DIMACS format file in a Networkx graph G
'''
def read_graph_from_dimacs(path):
    try:
        with open(path, 'r') as f:
            edges = []
            for line in f:
                if line[0] in ['c', '', '\n']: # Skip comments, header and empty
                    continue
                if line[0] in ['p']: 
                    n, m = re.split('\s+', line)[2:4] 
                    continue
                edge = re.split('\s+', line)
                edge = int(edge[1]) - 1, int(edge[2]) - 1
                edges.append(edge)
                
    except IOError as exc:
        raise

    G = nx.Graph()
    #G.add_nodes_from(range(0, int(n)))
    G.add_edges_from(edges)
    return G


def read_graph_from_mat(path, var, type='adjacency'):
    try:
        G = nx.Graph()
        data = io.loadmat(path)
        if 'adjacency' in type:
            G = nx.from_numpy_matrix(np.array(data[var])) 
        elif 'edgelist' in type:
            G.add_edges_from(np.array(data[var]))       
    except IOError as exc:
        raise

    return nx.convert_node_labels_to_integers(G)


'''
Read coefficients from a Nodal LP. The result is a dictionary 
with node labels as keys.
'''
def read_alpha_from_lp(filename):
    alphas = {}
    pattern = re.compile(r"\d+\sx\d+")
    for line in open(filename):
        match = re.findall(pattern, line)
        if match:
            match = match[0].split(' ')
            alphas[int(match[1][1:]) - 1] = int(match[0]) 
    return alphas

