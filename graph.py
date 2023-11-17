import networkx as nx
# import matplotlib.pyplot as plt
import scipy.io as io
import numpy as np
import re

def web(p, q):
    G = nx.Graph()
    edges = []
    for i in range(p):
        edges += [(i % p, j % p) for j in range(i + q , p + i - q +1)]
    G.add_edges_from(edges)
    # print('WEB: ',G.edges(), len(G.edges()))
    return G


def antiweb(p, q):
    G = nx.complement(web(p, q))
    # print('ANTIWEB: ',G.edges())
    return G


def wheel(n):
    return nx.wheel_graph(n)


def cycle(n):
    return nx.cycle_graph(n)


def random_graph(n, p, seed=1234):
    return nx.erdos_renyi_graph(n, p, seed=seed)


def random_perfect_graph(n, p, seed=1234):
    # Returns random perfect graph obtained as line graph of a random bipartite graph
    # n: the expected number of nodes in the line graph (i.e. number of edges in the bipartite graph)
    # p: density of the random bipartite graph
    
    # num_edges_complete_bip = int(n / p)
    # dim_partitions = int(math.floor(math.sqrt(num_edges_complete_bip)))
    # return nx.convert_node_labels_to_integers(nx.line_graph(nx.bipartite.random_graph(dim_partitions, dim_partitions, p, seed=seed)))
    return nx.convert_node_labels_to_integers(nx.line_graph(nx.bipartite.random_graph(n, n, p, seed=seed)))

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
    

def one_edge():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2)])
    return G

def test():
    G = nx.Graph()
    G.add_edges_from([(0,1), (0,2), (0,3), (0,4), (0,5), (0,6), (1, 2), (2,3), (3,4), (2,5), (2,6)])
    return G


def write_graph_to_dimacs(G, path):
    with open(path, "w") as f:
        # write the header
        f.write("p edge {} {}\n".format(G.number_of_nodes(), G.number_of_edges()))
        # now write all edges
        for u, v in G.edges():
            f.write("e {} {}\n".format(u + 1, v + 1))


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


def draw(G, path='', show=False, node_labels=True, figsize=(13, 13), X=[], height=0.05):
    pos = nx.spring_layout(G)
    plt.figure(figsize=figsize)
    nx.draw_networkx(G, with_labels=node_labels, labels={i: i for i in G.nodes()}, pos=pos)
    
    if X.any():                     # A solution is provided
        for p in pos:               # raise text pos
            pos[p][1] += height 
        nx.draw_networkx_labels(G, pos, labels={i: X[i] for i in G.nodes()}, font_color='red')

    if show:
        plt.show()

    if path:
        plt.savefig(path)

    plt.close()


def read_alpha_from_lp(filename):
    alphas = {}
    pattern = re.compile(r"\d+\sx\d+")
    for line in open(filename):
        match = re.findall(pattern, line)
        if match:
            match = match[0].split(' ')
            alphas[int(match[1][1:]) - 1] = int(match[0]) 
    return alphas

