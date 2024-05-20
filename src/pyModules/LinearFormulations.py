import re, os
from Graphs import *
import gurobi as gb


def FRAC(G, model_name, dir_path='', write_lp=True):
    path = os.path.join(dir_path, model_name + '.lp')
    frac = gb.Model()
    # frac.setParam('OutputFlag', 0)
    # frac.setParam('Threads', 80)
    x = frac.addVars([i+1 for i in G.nodes()], vtype=gb.GRB.BINARY, name='x', lb=0.0, ub=1.0)
    frac.setObjective(x.sum(), gb.GRB.MAXIMIZE)
    for i, j in G.edges():
        frac.addConstr(x[i+1] + x[j+1] <= 1, name='edge')
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


def compute_theta(filename, model_out_dir='', use_patience=False):
    # Options for the Admm
    options = {}
    options['tolerance'] = 1e-4
    options['max_iter'] = 1000000
    options['timelimit'] = 3600
    options['print_it'] = 1000
    options['debug'] = True

    file_tmp = 'tmp'

    G = read_graph_from_dimacs(os.path.join(work_dir, filename + '.stb'))

    if os.path.exists(os.path.join(model_out_dir, "%s_theta.json" % (filename))):
        # Read Theta coeff in JSON file
        with open(os.path.join(model_out_dir, filename + '_theta.json')) as f:
            j = json.load(f)
            theta = {}
            for i in j:
                theta[int(i)] = j[i][0]
        return theta

    theta = {}
    for j in G.nodes():
        H = G.subgraph(G.neighbors(j))
        H = nx.convert_node_labels_to_integers(H)
        sdpss.Theta_SDP2(H, file_tmp, model_out_dir=model_out_dir, model_out='adal', debug=False)
        h = loadmat(os.path.join(model_out_dir, file_tmp))
        P = Problem(h['A'], h['b'], h['C'], c=0)
        opt_sol, _, _, secs = ADMMsolver.solve(P, norm_bound=len(H.nodes()), options=options, use_patience=use_patience)
        theta[j] = [np.floor(-opt_sol.safe_dual), secs]
        print("Node %d : Theta %.2f" % (j, theta[j][0]), end='\r')
    # create json object from dictionary
    to_write = json.dumps(theta)
    f = open(os.path.join(model_out_dir, "%s_theta.json" % (filename)),"w")
    f.write(to_write)
    f.close()
    for i in theta:
        theta[i] = theta[i][0]
    return theta


def compute_alpha(filename, model_out_dir=''):
    G = read_graph_from_dimacs(os.path.join(work_dir, filename + '.stb'))
    if os.path.exists(os.path.join(model_out_dir, "%s_alpha.json" % (filename))):
        # Read Alpha coeff in JSON file
        with open(os.path.join(model_out_dir, filename + '_alpha.json')) as f:
            j = json.load(f)
            alpha = {}
            for i in j:
                alpha[int(i)] = j[i][0]
        return alpha
    
    alpha = {}
    for j in G.nodes():
        H = G.subgraph(G.neighbors(j))
        H = nx.convert_node_labels_to_integers(H)
        # cov = FRAC(G, 'tmp', write_lp=False)
        # cov.optimize()
        # alpha[j] = [cov.objVal, cov.Runtime]
        start = time.time()
        a = max(len(c) for c in nx.find_cliques(nx.complement(H)))
        end = time.time() - start
        alpha[j] = [a, end]
        print("Node %d : Alpha %.2f" % (j, alpha[j][0]), end='\r')
    # create json object from dictionary
    to_write = json.dumps(alpha)
    f = open(os.path.join(model_out_dir, "%s_alpha.json" % (filename)),"w")
    f.write(to_write)
    f.close()
    for i in alpha:
        alpha[i] = alpha[i][0]
    return alpha