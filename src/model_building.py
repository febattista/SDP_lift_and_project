import os, time, json, zipfile, itertools, sys
from pyModules.Graphs import *
from pyModules.SDPLifting import *
from pyModules.LinearFormulations import *

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
        print("--------------------------------------------------")
        print("  BUILDING LP FORMULATIONS ")
        print("--------------------------------------------------")
        # Create directories
        for path in [lp_dir, coeff_alpha_dir, coeff_theta_dir]:
            if not os.path.exists(path):
                os.makedirs(path)

        with os.scandir(graph_dir) as inst_it: 
            for instance in inst_it:
                if instance.name.lower().endswith('.stb'):
                    print("> Processing %s ..." %  instance.name)
                    graphname = os.path.splitext(instance.name)[0]
                    G = read_graph_from_dimacs(os.path.join(graph_dir, instance.name))

                    # NOD(G, gamma)
                    gamma = compute_gamma(G)
                    NOD(G, graphname + '_nod_gamma', gamma, lp_dir)

                    # NOD(G, theta)
                    theta = compute_theta(G, graphname, 
                                            model_out_dir=coeff_theta_dir, 
                                            use_patience=True)
                    NOD(G, graphname + '_nod_theta', theta, lp_dir)

                    # NOD(G, alpha)
                    alpha = compute_alpha(G, graphname, 
                                            model_out_dir=coeff_alpha_dir)
                    NOD(G, graphname + '_nod_alpha', alpha, lp_dir)

                    # QSTAB(G, C)
                    QSTABC(G, graphname + '_cov', dir_path=lp_dir)

                    # FRAC(G)
                    FRAC(G, graphname + '_edge', dir_path=lp_dir)
        

    
    if MAKE_SDP_MODELS:
        print("--------------------------------------------------")
        print("  BUILDING SDP MODELS ")
        print("--------------------------------------------------")
        # Create directory
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
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
                        
                        Theta_plus_SDP(G, filename, 
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
            