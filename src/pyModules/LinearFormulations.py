import os, json, subprocess
import gurobi as gb
from scipy.io import loadmat

from pyModules.Graphs import *
from pyModules.ADMMsolver import *
from pyModules.SDPLifting import *


def FRAC(G, model_name, dir_path='', write_lp=True):
    path = os.path.join(dir_path, model_name + '.lp')
    frac = gb.Model()
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
        with open(path) as f:
            found = False
            lines = []
            for line in f:
                if not found and 'x[' in line :
                    s = "obj: " + line
                    lines.append(s)
                    found = True
                else:
                    lines.append(line)
        with open(path, 'w') as f:
            for line in lines:
                f.write(line)

    return frac


def QSTAB(G, model_name, dir_path='', write_lp=True):
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
        with open(path) as f:
            found = False
            lines = []
            for line in f:
                if not found and 'x[' in line :
                    s = "obj: " + line
                    lines.append(s)
                    found = True
                else:
                    lines.append(line)
        with open(path, 'w') as f:
            for line in lines:
                f.write(line)

    return cov


def QSTABC(G, model_name, dir_path='', write_lp=True):
    cliques = greedy_clique_cover(G)
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
        with open(path) as f:
            found = False
            lines = []
            for line in f:
                if not found and 'x[' in line :
                    s = "obj: " + line
                    lines.append(s)
                    found = True
                else:
                    lines.append(line)
        with open(path, 'w') as f:
            for line in lines:
                f.write(line)

    return cov


def NOD(G, model_name, coeff, dir_path=''):
    path = os.path.join(dir_path, model_name + '.lp')
    nod = gb.Model()
    x = nod.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    nod.setObjective(x.sum(), gb.GRB.MAXIMIZE)

    for i in G.nodes():
        neighbors = [j+1 for j in G.neighbors(i)]
        nod.addConstr(x.sum(neighbors) + coeff[i]*x[i + 1] <= coeff[i], name='nod' + str(i + 1))

    nod.update()
    print('Saving LP at: %s' % path)
    nod.write(path)
    with open(path) as f:
        found = False
        lines = []
        for line in f:
            if not found and 'x[' in line :
                s = "obj: " + line
                lines.append(s)
                found = True
            else:
                lines.append(line)
    with open(path, 'w') as f:
        for line in lines:
            f.write(line)

    return nod

def compute_gamma(G):
    # Compute \Gamma's
    degs = {}
    for x, y in G.degree:
        degs[x] = y
    return degs


def compute_theta(G, graphname, model_out_dir='', use_patience=False, options=None):
    # Options for the Admm
    if not options:
        options = {}
        options['tolerance'] = 1e-4
        options['max_iter'] = 1000000
        options['timelimit'] = 3600
        options['print_it'] = 1000
        options['debug'] = False

    file_tmp = 'tmp'
    if os.path.exists(os.path.join(model_out_dir, "%s_theta.json" % (graphname))):
        # Read Theta coeff in JSON file
        with open(os.path.join(model_out_dir, graphname + '_theta.json')) as f:
            j = json.load(f)
            theta = {}
            for i in j:
                theta[int(i)] = j[i][0]
        return theta

    tmp_path = os.path.join(model_out_dir, file_tmp + '.mat')
    theta = {}
    n = len(G.nodes())
    for i, j in enumerate(G.nodes()):
        H = G.subgraph(G.neighbors(j))
        H = nx.convert_node_labels_to_integers(H)
        Theta_SDP(H, file_tmp, model_out_dir=model_out_dir, model_out='adal', debug=False)
        h = loadmat(tmp_path)
        P = Problem(h['A'], h['b'], h['C'], c=0)
        opt_sol, _, _, secs = ADMMsolver.solve(P, norm_bound=len(H.nodes()), options=options, use_patience=use_patience)
        theta[j] = [np.floor(-opt_sol.safe_dual), secs]
        print("Node %d of %d : Theta %.2f" % (i + 1, n, theta[j][0]), end='\r')
    print("\nDone!")
    # Clean up tmp file
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    # create json object from dictionary
    to_write = json.dumps(theta)
    f = open(os.path.join(model_out_dir, "%s_theta.json" % (graphname)),"w")
    f.write(to_write)
    f.close()
    for i in theta:
        theta[i] = theta[i][0]
    return theta


def cliquer_max_clique(filepath):
    command = ['../ThirdParty/cliquer-1.21/cl', filepath, '-q', '-q']
    start = time.time()
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    cltime = time.time() - start
    # Get the output as a string
    output = result.stdout
    alpha = float(output.split(",")[0].split("=")[1])
    return alpha, cltime


def compute_alpha(G, graphname, model_out_dir=''):
    if os.path.exists(os.path.join(model_out_dir, "%s_alpha.json" % (graphname))):
        # Read Alpha coeff in JSON file
        with open(os.path.join(model_out_dir, graphname + '_alpha.json')) as f:
            j = json.load(f)
            alpha = {}
            for i in j:
                alpha[int(i)] = j[i][0]
        return alpha
    tmp_path = os.path.join(model_out_dir, 'subgraph.stb')
    alpha = {}
    n = len(G.nodes())
    for i, j in enumerate(G.nodes()):
        H = G.subgraph(G.neighbors(j))
        H = nx.convert_node_labels_to_integers(H)
        write_graph_to_dimacs(nx.complement(H), tmp_path)
        a, cltime = cliquer_max_clique(tmp_path)
        alpha[j] = [a, cltime]
        print("Node %d of %d : Alpha %.2f" % (i + 1, n, alpha[j][0]), end='\r')
    print("\nDone!")
    # Clean up tmp file
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    # create json object from dictionary
    to_write = json.dumps(alpha)
    f = open(os.path.join(model_out_dir, "%s_alpha.json" % (graphname)),"w")
    f.write(to_write)
    f.close()
    for i in alpha:
        alpha[i] = alpha[i][0]
    return alpha