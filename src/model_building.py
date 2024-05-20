import os, time, json, zipfile, itertools
from pyModules.Graphs import *
from pyModules.SDPLifting import *
from pyModules.ADMMsolver import *

from parameters import *


def add_files_to_zip(zip_filename, files):
    with zipfile.ZipFile(zip_filename, 'a', compression=zipfile.ZIP_DEFLATED) as zip_file:
        for file in files:
            print('Deflating: %s' % file)
            zip_file.write(file, arcname=os.path.basename(file))
            if os.path.exists(file):
                os.remove(file)


for d in datasets:
    # Set up paths
    data_path = datasets[d]
    model_dir = os.path.join(data_path, 'models')
    graph_dir = os.path.join(data_path, 'graphs')
    lp_dir    = os.path.join(data_path, 'lp')
    coeff_alpha_dir = os.path.join(data_path, 'coeff_alpha')
    coeff_theta_dir = os.path.join(data_path, 'coeff_theta')
    
    if not os.path.exists(data_path):
        print("WARNING: path to Dataset %s does not exist! Skipping this step ..." % d)
        continue

    # Create directories
    for path in [model_dir, graph_dir, lp_dir, coeff_alpha_dir, coeff_theta_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    if MAKE_LP_FORMULATION:
        # TODO: CREATE LP FORMULATIONS
        pass
    
    if MAKE_SDP_MODELS:
        print("--------------------------------------------------")
        print("  BUILDING SDP MODELS ")
        print("--------------------------------------------------")
        # Create an empty zip file
        out_path = os.path.join(model_dir, 'models_%s.zip' % d.lower())
        if not os.path.exists(out_path):
            with zipfile.ZipFile(out_path, 'a') as zip_file:
                pass
        else:
            print("WARNING: %s already exists! Skipping this step ..." 
                  % out_path)
            continue

        if MAKE_TH_PLUS_MODELS:
            print("> Creating Theta_plus SDPs ...")
            with os.scandir(graph_dir) as inst_it: 
                for instance in inst_it:
                    if instance.name.lower().endswith('.stb'):
                        filename = os.path.splitext(instance.name)[0] + '_th+'
                        G = read_graph_from_dimacs(os.path.join(graph_dir, instance.name))
                        
                        Theta_plus_SDP2(G, filename, 
                                        model_out_dir=model_dir, 
                                        model_out='sdpnal', debug=True)
                        
                        files_to_zip = [os.path.join(model_dir, filename + '.mat')]
                        
                        add_files_to_zip(out_path, files_to_zip)
            print("> Done!")

        with os.scandir(lp_dir) as inst_it: 
            for instance in inst_it:
                if instance.name.lower().endswith('.lp'):
                    filename = os.path.splitext(instance.name)[0]
                    # SDP formulations
                    m_plus_lifting(instance.name, 
                                    model_out_dir=model_dir, 
                                    work_dir=lp_dir, 
                                    skip_func=(lovasz_schrijver_filter if 'edge' in instance.name else None),
                                    do_split=DO_SPLIT)
                    
                    files_to_zip = [os.path.join(model_dir, ''.join(s)) \
                                    for s in itertools.product([filename], \
                                    ['_1.mat', '_2.mat', '_3.mat', '_4.mat', '_5.mat'] \
                                    if DO_SPLIT else ['.mat'])]
                    
                    add_files_to_zip(out_path, files_to_zip)
        print("> Done!")
    print("> Task completed for Dataset: %s" % d)
                    

# model_out_dir = os.path.join('StableSets', 'DIMACS', 'models')
# json_dir = os.path.join('StableSets', 'DIMACS', 'coeff_theta')
# json_dir_1 = os.path.join('StableSets', 'DIMACS', 'coeff_alpha')
# work_dir = os.path.join('StableSets', 'DIMACS')
# lp_dir = os.path.join('StableSets', 'DIMACS', 'lp')

# for path in [model_out_dir, json_dir, json_dir_1, work_dir, lp_dir]:
#     if not os.path.exists(path):
#         os.makedirs(path)

# out_file = 'models_dimacs.zip'
# out_path = os.path.join(model_out_dir, out_file)

# if not os.path.exists(out_file):
#     with zipfile.ZipFile(out_path, 'a') as zip_file:
#         # Create an empty zip file
#         pass

# graphs = [ 'brock200_1' ]
#         # , 'brock200_2', 'brock200_3', 'brock200_4', \
#         # 'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', \
#         # 'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', \
#         # 'C125-9', 'C250-9', 'C500-9', \
#         # 'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC500-5',\
#         # 'MANN_a9', 'MANN_a27',
#         # 'johnson32-2-4',\
#         # 'keller4', # 'keller5', \
#         # 'p_hat300-1', 'p_hat300-2', 'p_hat300-3', \
#         # 'p_hat500-1', 'p_hat500-2', 'p_hat500-3', \
#         # 'p_hat700-1', 'p_hat700-2', 'p_hat700-3', \
#         # 'sanr200_0.7', 'sanr200_0.9',  \
#         # 'sanr400_0.5', 'sanr400_0.7']


# for g in graphs:
#     filename = g
#     G = read_graph_from_dimacs(os.path.join(work_dir, filename + '.stb'))
#     sdpss.Theta_plus_SDP2(G, filename + '_th+', model_out_dir=model_out_dir, model_out='sdpnal', debug=True)

#     # Compute \Gamma's
#     degs = {}
#     for x, y in G.degree:
#         degs[x] = y

#     # Compute \alpha's
#     # Read the file if pre-computed coefficients are provided
#     if os.path.exists(os.path.join(lp_dir, filename + '_empty_nod.lp')):
#         alphas = read_alpha_from_lp(os.path.join(lp_dir, filename + '_empty_nod.lp'))
#     else:
#         alphas = compute_alpha(filename, model_out_dir=json_dir_1)

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