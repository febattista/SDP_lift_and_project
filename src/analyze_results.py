import os
import re
import zipfile
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import savemat, loadmat
import time
import itertools, collections
import json

from pyModules.SDPLifting import *


def count_values(input_dict, idx):
    output_dict = {}
    for i in idx:
        value = input_dict[i[0]]
        if value not in output_dict:
            output_dict[value] = 1
        else:
            output_dict[value] += 1
    return output_dict

data_dir   = os.path.join('StableSets', 'DIMACS', 'coeff_theta')
# data_dir   = os.path.join('StableSets', 'Alpha')
# data_dir   = os.path.join('StableSets', 'DIMACS', 'violated_cuts')
# data_dir   = os.path.join('StableSets', 'Random-np', 'violated_cuts')
lp_dir     = os.path.join('StableSets', 'DIMACS', 'lp')
# lp_dir     = os.path.join('StableSets', 'Random-np', '%s')
work_dir   = os.path.join('StableSets', 'DIMACS')
# work_dir   = os.path.join('StableSets', 'Random-np')
model_path = os.path.join('StableSets', 'DIMACS', 'models')
# model_path = os.path.join('StableSets', 'Random-np', 'models')

# csv_filename = 'DIMACS_cutting_planes'
csv_filename = 'random_cutting_planes'

zip_name = 'models_dimacs.zip'

pieces = 5

results = collections.defaultdict(list)

graphs = [ 
        # 'brock200_1',
        # 'brock200_2', 'brock200_3', 'brock200_4', \
        # 'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', \
        'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', \
        # 'C125-9', 'C250-9', 
        'C500-9', \
        # 'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC500-5',\
        # 'MANN_a9', 'MANN_a27',
        # 'johnson32-2-4',\
        # 'keller4', # 'keller5', \
        'p_hat300-1', 'p_hat300-2', 'p_hat300-3', \
        'p_hat500-1', 'p_hat500-2', 'p_hat500-3', \
        'p_hat700-1', 'p_hat700-2', 'p_hat700-3', \
        # 'sanr200_0.7', 'sanr200_0.9',  \
        'sanr400_0.5', 'sanr400_0.7']

for graph in graphs:
    results["instance"].append(graph)
    print(graph)
    data_dir   = os.path.join('StableSets', 'DIMACS', 'coeff_theta')
    filepath = os.path.join(data_dir, graph + '_theta.json')
    # Open the JSON file for reading
    with open(filepath, 'r') as file:
        # # Load JSON data
        data = json.load(file)
        times = []
        for _, item in data.items():
            times.append(item[1])
        
        times = np.array(times)
        results["min_theta"].append(np.min(times))
        results["max_theta"].append(np.max(times))
        results["average_theta"].append(np.average(times))

    data_dir   = os.path.join('StableSets', 'Alpha')
    filepath = os.path.join(data_dir, graph + '.stb.log')
    with open(filepath, 'r') as file:
        times = []

        lines = file.readlines()
        start = next(i for i, line in enumerate(lines) if "Rank" in line)
        
        for line in lines[start + 1:]:
            if "----------------------------------------" in line:
                break
            times.append(float(line.split("|")[5].strip()))
        
        times = np.array(times)
        results["min_alpha"].append(np.min(times))
        results["max_alpha"].append(np.max(times))
        results["average_alpha"].append(np.average(times))

df_result = pd.DataFrame(results)
df_result.set_index('instance', inplace=True)
df_result = df_result.transpose()

def format_value(val):
    if val < 1e-2:
        return '$<0.01$'
    else:
        return '{:.2f}'.format(val)

print(df_result.style.format(format_value, precision=2).map_index(
    lambda v: "rotatebox:{80}--rwrap;", axis=1
    ).to_latex())


# n = [300, 275, 250, 225, 200, 175, 150]
# d = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
# i = [1, 2, 3, 4, 5]

# basename = 'G_%d_%1.1f_%d'

# idx = [g for g in itertools.product(n, d, i)]

# graphs = [basename % graph for graph in itertools.product(n, d, i)]

# relaxations = ['edge', 'cov', 'nod_alpha', 'nod_theta', 'nod_gamma']

# # with zipfile.ZipFile(os.path.join(model_path, zip_name), 'r') as zip_ref:
# for graph in graphs: 
#     results["instance"].append(graph)
#     print(graph)
#     for rel in relaxations:
#         filename = '%s_%s_viol_test.mat' % (graph, rel)
#         dirname = graph.split('_')[1]
        
#         if '275' in dirname and 'cov' in rel:
#             lpname = '%s_%s.lp' % (graph, 'clique')
#         else:
#             lpname = '%s_%s.lp' % (graph, rel)

#         filepath = os.path.join(data_dir, filename)
#         lppath = os.path.join(lp_dir % dirname, lpname)
#         if os.path.exists(filepath): # and os.path.exists(lppath):
#             if 'nod' in rel:
#                 maps = map_nod_constraints(lppath)
#             elif 'cov' in rel:
#                 maps = map_clique_constraints(lppath)
#             elif 'edge' in rel:
#                 maps = map_edge_constraints(lppath)
#             else:
#                 print('Unknown relaxation!')
#             # print(set(maps.values()))
#             h = loadmat(filepath)
#             added_cuts = count_values(maps, h['added_cuts_idx'])
#             for i in set(maps.values()):
#                 if i in added_cuts:
#                     results[rel + '_' + i].append(added_cuts[i])
#                 else:
#                     results[rel + '_' + i].append(0)

#             # Check for duplicate constraints added during the cutting plane
#             # model_name = '%s_%s_%d.mat'
#             # Bt = []
#             # for i in range(2, pieces + 1):
#             #     zip_ref.extract(os.path.join(model_path, model_name % (graph, rel, i)), 
#             #                     os.path.join(model_path, model_name % (graph, rel, i)))
#             #     Bt.append(loadmat(os.path.join(model_path, model_name % (graph, rel, i)))['Bt_' + i])
            
#             # Bt = sparse.hstack(Bt, format='csc')
                    
                    


# df_result = pd.DataFrame(results)

# df_result.set_index('instance', inplace=True)

# col_mapping_to_labels = {
#     'nod_gamma_bound2' : r'\eqref{eq:liftbound2}-\eqref{eq:liftbound3}', 
#     'nod_gamma_bound4' : r'\eqref{eq:liftbound4}', 
#     'nod_gamma_lift_nod1:1' : r'\eqref{nod1}',
#     'nod_gamma_lift_nod2:1' : r'\eqref{nod2}',
#     'nod_gamma_lift_nod3:1' : r'\eqref{nod3}',
#     'nod_gamma_lift_nod1:2' : r'\eqref{nod4}',
#     'nod_gamma_lift_nod2:2' : r'\eqref{nod5}',
#     'nod_gamma_lift_nod3:2' : r'\eqref{nod6}',
#     'nod_theta_bound2' : r'\eqref{eq:liftbound2}-\eqref{eq:liftbound3}', 
#     'nod_theta_bound4' : r'\eqref{eq:liftbound4}', 
#     'nod_theta_lift_nod1:1' : r'\eqref{nod1}',
#     'nod_theta_lift_nod2:1' : r'\eqref{nod2}',
#     'nod_theta_lift_nod3:1' : r'\eqref{nod3}',
#     'nod_theta_lift_nod1:2' : r'\eqref{nod4}',
#     'nod_theta_lift_nod2:2' : r'\eqref{nod5}',
#     'nod_theta_lift_nod3:2' : r'\eqref{nod6}',
#     'nod_alpha_bound2' : r'\eqref{eq:liftbound2}-\eqref{eq:liftbound3}', 
#     'nod_alpha_bound4' : r'\eqref{eq:liftbound4}', 
#     'nod_alpha_lift_nod1:1' : r'\eqref{nod1}',
#     'nod_alpha_lift_nod2:1' : r'\eqref{nod2}',
#     'nod_alpha_lift_nod3:1' : r'\eqref{nod3}',
#     'nod_alpha_lift_nod1:2' : r'\eqref{nod4}',
#     'nod_alpha_lift_nod2:2' : r'\eqref{nod5}',
#     'nod_alpha_lift_nod3:2' : r'\eqref{nod6}',
#     'cov_bound2' : r'\eqref{eq:liftbound2}-\eqref{eq:liftbound3}', 
#     'cov_bound4' : r'\eqref{eq:liftbound4}', 
#     'cov_lift_clq1:1' : r'\eqref{clq1}',
#     'cov_lift_clq2:1' : r'\eqref{clq2}',
#     'cov_lift_clq1:2' : r'\eqref{clq3}',
#     'cov_lift_clq2:2' : r'\eqref{clq4}',
#     'edge_bound2' : r'\eqref{eq:liftbound2}-\eqref{eq:liftbound3}', 
#     'edge_bound4' : r'\eqref{eq:liftbound4}', 
#     'edge_lift_edge1:1' : r'\eqref{eq:LS1}',
#     'edge_lift_edge1:2' : r'\eqref{eq:LS2}'
# }

# cols_ord = [
#     'nod_gamma_bound2', 
#     'nod_gamma_bound4', 
#     'nod_gamma_lift_nod1:1',
#     'nod_gamma_lift_nod2:1',
#     'nod_gamma_lift_nod3:1',
#     'nod_gamma_lift_nod1:2',
#     'nod_gamma_lift_nod2:2',
#     'nod_gamma_lift_nod3:2',
#     'nod_theta_bound2', 
#     'nod_theta_bound4', 
#     'nod_theta_lift_nod1:1',
#     'nod_theta_lift_nod2:1',
#     'nod_theta_lift_nod3:1',
#     'nod_theta_lift_nod1:2',
#     'nod_theta_lift_nod2:2',
#     'nod_theta_lift_nod3:2',
#     'nod_alpha_bound2', 
#     'nod_alpha_bound4', 
#     'nod_alpha_lift_nod1:1',
#     'nod_alpha_lift_nod2:1',
#     'nod_alpha_lift_nod3:1',
#     'nod_alpha_lift_nod1:2',
#     'nod_alpha_lift_nod2:2',
#     'nod_alpha_lift_nod3:2',
#     'cov_bound2', 
#     'cov_bound4', 
#     'cov_lift_clq1:1',
#     'cov_lift_clq2:1',
#     'cov_lift_clq1:2',
#     'cov_lift_clq2:2',
#     'edge_bound2',
#     'edge_bound4',
#     'edge_lift_edge1:1',
#     'edge_lift_edge1:2'
# ]

# df_result = df_result[cols_ord]

# df_result.to_csv(os.path.join(work_dir, csv_filename + '.csv'), index=True)

# # Identify columns with all zeros
# columns_to_remove = df_result.columns[df_result.eq(0).all()]

# # Remove columns with all zeros
# df_result = df_result.loc[:, ~df_result.columns.isin(columns_to_remove)]

# df_result.index = pd.MultiIndex.from_tuples(idx, names=('n', 'd', 'i'))

# df_result.groupby(by=['n', 'd']).mean().style \
#     .format_index(formatter={'n' : '{:.0f}', 'd': '{:.1f}'}) \
#     .format(precision=0) \
#     .to_latex(buf=os.path.join(work_dir,csv_filename + '_table.tex'))

# print("HEADER FOR LATEX TABLE")
# s = 'Graph & '
# for col in cols_ord:
#     if col in df_result.columns:
#         s += col_mapping_to_labels[col] + "& "

# s = s[:-2] + "\\\\"

# print(s)


# inst_name_pattern = re.compile(r'(.*?)_([a-zA-Z]+_?[a-zA-Z]*)_viol_test\.mat')

# with os.scandir(data_dir) as dataset_it:
#     for d_entry in dataset_it:
#         if d_entry.name.endswith(".mat"):
#             match = inst_name_pattern.match(d_entry.name) 
#             instance =  "" if not match else match.group(1)
#             if instance:
#                 results["instance"].append(instance)

