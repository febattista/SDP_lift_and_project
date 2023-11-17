import os, time, json, zipfile, itertools
from graph import *
from m_plus_building import *
from scipy.io import savemat, loadmat
from ADMMsolver import *

import gurobipy as gb
import SDP_StableSet as sdpss

def FRAC(G, model_name, dir_path='', write_lp=True):
    path = os.path.join(dir_path, model_name + '.lp')
    frac = gb.Model()
    frac.setParam('OutputFlag', 0)
    frac.setParam('Threads', 80)
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


def COV(G, model_name, dir_path='', write_lp=True):
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

def NSTAB(G, model_name, coeff, dir_path=''):
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

def add_files_to_zip(zip_filename, files):
    with zipfile.ZipFile(zip_filename, 'a', compression=zipfile.ZIP_DEFLATED) as zip_file:
        for file in files:
            print('Deflating: %s' % file)
            zip_file.write(file)
            if os.path.exists(file):
                os.remove(file)

model_out_dir = os.path.join('StableSets', 'DIMACS', 'models')
json_dir = os.path.join('StableSets', 'DIMACS', 'coeff_theta')
json_dir_1 = os.path.join('StableSets', 'DIMACS', 'coeff_alpha')
work_dir = os.path.join('StableSets', 'DIMACS')
lp_dir = os.path.join('StableSets', 'DIMACS', 'lp')

for path in [model_out_dir, json_dir, json_dir_1, work_dir, lp_dir]:
    if not os.path.exists(path):
        os.makedirs(path)

out_file = 'models_dimacs.zip'
out_path = os.path.join(model_out_dir, out_file)

if not os.path.exists(out_file):
    with zipfile.ZipFile(out_path, 'a') as zip_file:
        # Create an empty zip file
        pass

else:
    out_path = out_path[:-4] + '_1.zip'
    with zipfile.ZipFile(out_path, 'w') as zip_file:
        # Create an empty zip file
        pass

graphs = [ 'brock200_1' ]
        # , 'brock200_2', 'brock200_3', 'brock200_4', \
        # 'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', \
        # 'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', \
        # 'C125-9', 'C250-9', 'C500-9', \
        # 'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC500-5',\
        # 'MANN_a9', 'MANN_a27',
        # 'johnson32-2-4',\
        # 'keller4', # 'keller5', \
        # 'p_hat300-1', 'p_hat300-2', 'p_hat300-3', \
        # 'p_hat500-1', 'p_hat500-2', 'p_hat500-3', \
        # 'p_hat700-1', 'p_hat700-2', 'p_hat700-3', \
        # 'sanr200_0.7', 'sanr200_0.9',  \
        # 'sanr400_0.5', 'sanr400_0.7']


for g in graphs:
    filename = g
    G = read_graph_from_dimacs(os.path.join(work_dir, filename + '.stb'))
    sdpss.Theta_plus_SDP2(G, filename + '_th+', model_out_dir=model_out_dir, model_out='sdpnal', debug=True)

    # Compute \Gamma's
    degs = {}
    for x, y in G.degree:
        degs[x] = y

    # Compute \alpha's
    # Read the file if pre-computed coefficients are provided
    if os.path.exists(os.path.join(lp_dir, filename + '_empty_nod.lp')):
        alphas = read_alpha_from_lp(os.path.join(lp_dir, filename + '_empty_nod.lp'))
    else:
        alphas = compute_alpha(filename, model_out_dir=json_dir_1)

    # Compute \theta's
    theta = compute_theta(filename, model_out_dir=json_dir, use_patience=True)

    # LP formulations
    FRAC(G, filename + '_edge', dir_path=lp_dir)
    COV(G, filename + '_cov', dir_path=lp_dir)
    NSTAB(G, filename + '_nod_gamma', degs, lp_dir)
    NSTAB(G, filename + '_nod_theta', theta, lp_dir)
    NSTAB(G, filename + '_nod_alpha', alphas, lp_dir)

    # SDP formulations
    do_split = True
    m_plus_lifting(filename + '_nod_gamma.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
    m_plus_lifting(filename + '_nod_theta.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
    m_plus_lifting(filename + '_nod_alpha.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
    m_plus_lifting(filename + '_cov.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
    m_plus_lifting(filename + '_edge.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split, lift_bounds=False, skip_func=lovasz_schrijver_filer)
    
    files_to_zip = [os.path.join(model_out_dir, ''.join(s)) for s in itertools.product([filename], \
                    ['_nod_gamma', '_nod_theta', '_nod_alpha', '_cov', '_edge'], \
                    ['_1.mat', '_2.mat', '_3.mat', '_4.mat', '_5.mat'])]
    
    add_files_to_zip(out_path, files_to_zip)

###################################################################
# -------------------- BUILDING DIMACS MODELS ---------------------
###################################################################
# model_out_dir = os.path.join('StableSets', 'DIMACS', 'models')
# json_dir = os.path.join('StableSets', 'DIMACS', 'coeff_theta')
# json_dir_1 = os.path.join('StableSets', 'DIMACS', 'coeff_alpha')
# work_dir = os.path.join('StableSets', 'DIMACS')
# lp_dir = os.path.join('StableSets', 'DIMACS', 'lp')

# out_file = 'models_dimacs.zip'
# out_path = os.path.join(model_out_dir, out_file)

# if not os.path.exists(out_file):
#     with zipfile.ZipFile(out_path, 'w') as zip_file:
#         # Create an empty zip file
#         pass

# else:
#     out_path = out_path[:-4] + '_1.zip'
#     with zipfile.ZipFile(out_path, 'w') as zip_file:
#         # Create an empty zip file
#         pass

# graphs = ['brock200_1', 'brock200_2', 'brock200_3', 'brock200_4', \
#         'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', \
#         'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', \
#         'C125-9', 'C250-9', 'C500-9', \
#         'DSJC125.1', 'DSJC125.5', \
#         'MANN_a9', 'MANN_a27', \
#         'keller4', 'keller5', \
#         'p_hat300-1', 'p_hat300-2', 'p_hat300-3', \
#         'p_hat500-1', 'p_hat500-2', 'p_hat500-3', \
#         'p_hat700-1', 'p_hat700-2', 'p_hat700-3', \
#         'sanr200_0.7', 'sanr200_0.9',  \
#         'sanr400_0.5', 'sanr400_0.7']

# for g in graphs:
#     filename = g
#     G = read_graph_from_dimacs(os.path.join(work_dir, filename + '.stb'))
#     # Compute \Gamma's
#     degs = {}
#     for x, y in G.degree:
#         degs[x] = y

#     # Compute \alpha's
#     alphas = compute_alpha(filename, model_out_dir=json_dir_1)

#     # Compute \theta's
#     theta = compute_theta(filename, model_out_dir=json_dir, use_patience=True)

#     # LP formulations
#     FRAC(G, filename + '_edge', dir_path=lp_dir)
#     COV(G, filename + '_cov', dir_path=lp_dir)
#     NSTAB(G, filename + '_nod_gamma', degs, lp_dir)
#     NSTAB(G, filename + '_nod_theta', theta, lp_dir)
#     NSTAB(G, filename + '_nod_alpha', alphas, lp_dir)

#     # SDP formulations
#     do_split = True
#     m_plus_lifting(filename + '_nod_gamma.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
#     m_plus_lifting(filename + '_nod_theta.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
#     m_plus_lifting(filename + '_nod_alpha.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
#     m_plus_lifting(filename + '_cov.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split)
#     m_plus_lifting(filename + '_edge.lp', model_out_dir=model_out_dir, work_dir=lp_dir, do_split=do_split, lift_bounds=False, skip_func=lovasz_schrijver_filer)
    
#     files_to_zip = [os.path.join(model_out_dir, ''.join(s)) for s in itertools.product([filename], \
#                     ['_nod_gamma', '_nod_theta', '_nod_alpha', '_cov', '_edge'], \
#                     ['_1.mat', '_2.mat', '_3.mat', '_4.mat', '_5.mat'])]
    
#     add_files_to_zip(out_path, files_to_zip)
###################################################################

###################################################################
# -------------------- BUILDING RANDOM MODELS ---------------------
###################################################################

# model_out_dir = os.path.join('StableSets', 'Random-np', 'models')
# json_dir = os.path.join('StableSets', 'Random-np', 'coeff_theta')
# json_dir_1 = os.path.join('StableSets', 'Random-np', 'coeff_alpha')
# work_dir = os.path.join('StableSets', 'Random-np')

# out_file = 'models_random.zip'
# out_path = os.path.join(model_out_dir, out_file)

# if not os.path.exists(out_path):
#     with zipfile.ZipFile(out_path, 'w') as zip_file:
#         # Create an empty zip file
#         pass

# else:
#     out_path = out_path[:-4] + '_1.zip'
#     with zipfile.ZipFile(out_path, 'w') as zip_file:
#         # Create an empty zip file
#         pass

# for n in [300, 275, 250, 225, 200, 175, 150]:
#     for d in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
#         for i in [1, 2, 3, 4, 5]:
#             work_dir_i = os.path.join(work_dir, '%d' % n)
#             filename = 'G_%d_%1.1f_%d' % (n, d, i)
#             G = read_graph_from_dimacs(os.path.join(work_dir_i, filename + '.stb'))
#             # Compute \Gamma's
#             degs = {}
#             for x, y in G.degree:
#                 degs[x] = y

#             # Compute \alpha's
#             alphas = {}
#             for i in G.nodes():
#                 H = nx.complement(G.subgraph(G.neighbors(i)))
#                 alphas[i] = float(nx.graph_clique_number(H))

#             # Read Theta coeff in JSON file
#             with open(os.path.join(json_dir, filename + '_theta.json')) as f:
#                 j = json.load(f)
#                 theta = {}
#                 tot_time = 0
#                 for i in j:
#                     theta[int(i)] = j[i][0]
#                     tot_time += j[i][1]

#             # LP formulations
#             NSTAB(G, filename + '_nod_gamma', degs,   work_dir_i)
#             NSTAB(G, filename + '_nod_theta', theta,  work_dir_i)
#             NSTAB(G, filename + '_nod_alpha', alphas, work_dir_i)

#             # SDP formulations
#             m_plus_lifting(filename + '_nod_gamma.lp', model_out_dir=model_out_dir, work_dir=work_dir_i)
#             m_plus_lifting(filename + '_nod_theta.lp', model_out_dir=model_out_dir, work_dir=work_dir_i)
#             m_plus_lifting(filename + '_nod_alpha.lp', model_out_dir=model_out_dir, work_dir=work_dir_i)
#             if n == 275:
#                 m_plus_lifting(filename + '_clique.lp', model_out_dir=model_out_dir, work_dir=work_dir_i)
#             else:
#                 m_plus_lifting(filename + '_cov.lp', model_out_dir=model_out_dir, work_dir=work_dir_i)
#             m_plus_lifting(filename + '_edge.lp', model_out_dir=model_out_dir, work_dir=work_dir_i, lift_bounds=False, skip_func=lovasz_schrijver_filer)
            
#             files_to_zip = [os.path.join(model_out_dir, ''.join(s)) for s in itertools.product([filename], \
#                     ['_nod_gamma.mat', '_nod_theta.mat', '_nod_alpha.mat', \
#                      '_clique.mat' if n == 275 else '_cov.mat', '_edge.mat'])]
    
#             add_files_to_zip(out_path, files_to_zip)

###################################################################