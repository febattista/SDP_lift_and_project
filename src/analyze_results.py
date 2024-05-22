import os, re, collections, json
import numpy as np
import pandas as pd
import json

from parameters import *

pd.options.mode.copy_on_write = True


# function to compute the percentage gap to the target column ('alpha')
def compute_gaps(dfin, columns, target='alpha'):
    l = 0
    for column in columns:
        dfin['gap' + column] = 100*(dfin[column] - dfin[target])/dfin[target]
        l+=1
    return dfin.iloc[:, l:]


# function to compute the number of violated cuts for each class
def count_values(input_dict, idx):
    output_dict = {}
    for i in idx:
        value = input_dict[i[0]]
        if value not in output_dict:
            output_dict[value] = 1
        else:
            output_dict[value] += 1
    return output_dict


# function to create a style formatter for index
def get_idx_formatter(dataset):
    formatter = {}
    if 'random' in dataset.lower():
        formatter['n'] = '{:.0f}'
        formatter['d'] = '{:.1f}'
        # formatter['i'] = '{:.0f}'
    elif 'dimacs' in dataset.lower():
        formatter['n'] = '{:s}'
    return formatter


# function to create a style formatter for columns
def get_col_formatter(columns):
    formatter = {}
    for c in columns:
        if 'time' in c:
            formatter[c] = '{:.2f}'
        elif 'rounds' in c:
            formatter[c] = '{:.0f}'
        elif 'gap' in c:
            formatter[c] = '{:.3f}'
        elif 'obj' in c or 'TH+' == c:
            formatter[c] = '{:.2f}'
    return formatter


# function to set the index based on dataset
def make_index(df, dataset):
    if 'random' in dataset.lower():
        # Make the triplet (n, d, i) the index
        # Assuming the graph's name is G_n_d_i
        df['n'] = df['n'].str.replace('G_', '')
        df[['n', 'd', 'i']] = df['n'].str.split('_', expand=True).astype(float)
        df.set_index(['n', 'd', 'i'], inplace=True)

    elif 'dimacs' in dataset.lower():
        # Just set the name of the graph as index
        df.set_index(['n'], inplace=True)


# Columns divided by category
columns   = ['TH+', 'NOD1_obj', 'NOD2_obj', 'NOD3_obj', 'COV_obj', 'EDGE_obj']
gaps      = ['gapTH+', 'gapNOD1_obj', 'gapNOD2_obj', 'gapNOD3_obj', 'gapCOV_obj', 'gapEDGE_obj']
cuts      = ['NOD1_cut', 'NOD2_cut', 'NOD3_cut', 'COV_cut', 'EDGE_cut']
rounds    = ['NOD1_rounds', 'NOD2_rounds', 'NOD3_rounds', 'COV_rounds', 'EDGE_rounds']
cpu_times = ['TH+_time', 'NOD1_time', 'NOD2_time', 'NOD3_time', 'COV_time', 'EDGE_time']

# Ordered columns to be displayed
# All 
all_cols = ['TH+', 'gapTH+', 'TH+_time', \
        'NOD1_obj', 'gapNOD1_obj', 'NOD1_rounds', 'NOD1_time', \
        'NOD2_obj', 'gapNOD2_obj', 'NOD2_rounds', 'NOD2_time', \
        'NOD3_obj', 'gapNOD3_obj', 'NOD3_rounds', 'NOD3_time', \
        'COV_obj', 'gapCOV_obj', 'COV_rounds', 'COV_time', \
        'EDGE_obj', 'gapEDGE_obj', 'EDGE_rounds', 'EDGE_time']
# UBs and gaps
ub_gaps = ['TH+', 'gapTH+', \
           'EDGE_obj', 'gapEDGE_obj', \
           'COV_obj', 'gapCOV_obj', \
           'NOD1_obj', 'gapNOD1_obj', \
           'NOD2_obj', 'gapNOD2_obj', \
           'NOD3_obj', 'gapNOD3_obj']
# CPU times and rounds
cpu_rounds = ['TH+_time', \
              'EDGE_rounds', 'EDGE_time', \
              'COV_rounds', 'COV_time', \
              'NOD1_rounds', 'NOD1_time', \
              'NOD2_rounds', 'NOD2_time', \
              'NOD3_rounds', 'NOD3_time']

for d in datasets:
    # Set up paths
    print("--------------------------------------------------")
    print("  DATASET %s " % d)
    print("--------------------------------------------------")
    data_path = datasets[d]
    model_dir = os.path.join(data_path, 'models')
    graph_dir = os.path.join(data_path, 'graphs')
    lp_dir    = os.path.join(data_path, 'lp')
    coeff_alpha_dir = os.path.join(data_path, 'coeff_alpha')
    coeff_theta_dir = os.path.join(data_path, 'coeff_theta')
    tables_dir = os.path.join(data_path, 'tables')

    res_csv = "results_%s.csv" % d.lower()
    aux_csv = "aux_data_%s.csv" % d.lower()
    
    if not os.path.exists(data_path):
        print("WARNING: path to Dataset %s does not exist! Skipping this step ..." % d)
        continue

    if not os.path.exists(os.path.join(data_path, res_csv)):
        print("WARNING: %s not found! Skipping this step ..." % res_csv)
        continue

    if not os.path.exists(tables_dir):
        os.makedirs(tables_dir)

    print("> Reading: %s" % res_csv)
    # Read csv and auxiliary files into dataframes
    with open(os.path.join(data_path, res_csv)) as f:
        l = []
        for line in f:
            l.append(re.sub(r'\s+', '\t',line).split('\t')[:-1])

        df = pd.DataFrame(l[1:], columns=l[0])
    
    print("> Reading: %s" % res_csv)
    with open(os.path.join(data_path, aux_csv)) as f:
        l = []
        for line in f:
            l.append(re.sub(r'\s+', '\t',line).split('\t')[:-1])

        aux = pd.DataFrame(l[1:], columns=l[0])

    df = df.merge(aux, on='n', how='left')

    target_cols = list(set(df.columns) - set('n'))
    df[target_cols] = df[target_cols].astype(float)

    # Compute percentage gaps
    compute_gaps(df, (['NOD_obj'] if 'random' in d.lower() else []) + columns)
    # Compute number of cutting plane's rounds
    df[rounds] = np.ceil(df[cuts] / CUTTING_PLANE_MAX_CUTS_PER_ITER)

    # Save the raw processed table as csv
    print("> Saving table: %s" % ('raw_results_%s.csv' % d.lower()))
    target = df[list(aux.columns) + (['gapNOD_obj'] if 'random' in d.lower() else []) + all_cols]
    make_index(target, d)
    filepath = os.path.join(tables_dir, 'raw_results_%s.csv' % d.lower())
    target.to_csv(filepath)
    
    # Format the output for LaTeX tables
    formatter_col = get_col_formatter(list(target.columns))
    formatter_idx = get_idx_formatter(d)
    
    print("> Saving table: %s" % ('ub_gaps_%s.tex' % d.lower()))
    filepath = os.path.join(tables_dir, 'ub_gaps_%s.tex' % d.lower())
    # If random graphs, aggregate over n and d and compute the mean
    if 'random' in d.lower():
        # Save UBs and gaps LaTeX table
        target[['gapNOD_obj'] +  ub_gaps].groupby(by=['n', 'd']).mean() \
            .style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .highlight_min(subset=gaps, axis=1, props="textbf:--rwrap;") \
            .highlight_min(subset=columns[1:], axis=1, props="textbf:--rwrap;") \
            .to_latex(buf=filepath)
    else :
        # Save UBs and gaps LaTeX table
        target[ub_gaps].style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .highlight_min(subset=gaps, axis=1, props="textbf:--rwrap;") \
            .highlight_min(subset=columns[1:], axis=1, props="textbf:--rwrap;") \
            .to_latex(buf=filepath)
    
    print("> Saving table: %s" % ('cpu_rounds_%s.tex' % d.lower()))
    filepath = os.path.join(tables_dir, 'cpu_rounds_%s.tex' % d.lower())
    if 'random' in d.lower():
        # Save CPU Times and rounds
        target[cpu_rounds].groupby(by=['n', 'd']).mean() \
            .style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .to_latex(buf=filepath)
    
    else:
        # Save CPU Times and rounds
        target[cpu_rounds].style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .to_latex(buf=filepath)

    if os.path.exists(coeff_theta_dir) and os.path.exists(coeff_alpha_dir):
        print("> Reading Theta coefficients ...")
        results_th = collections.defaultdict(list)
        results_alpha = collections.defaultdict(list)
        with os.scandir(coeff_theta_dir) as inst_it: 
            for instance in inst_it:
                if instance.name.lower().endswith('.json'):
                    graphname = os.path.splitext(instance.name)[0].replace('_theta', '')
                    results_th["n"].append(graphname)
                    filepath = os.path.join(coeff_theta_dir, instance.name)
                    # Open the JSON file for reading
                    with open(filepath, 'r') as file:
                        # Load JSON data
                        data = json.load(file)
                        times = []
                        for _, item in data.items():
                            times.append(item[1])
                        
                        times = np.array(times)
                        results_th["min_theta"].append(np.min(times))
                        results_th["max_theta"].append(np.max(times))
                        results_th["average_theta"].append(np.average(times))

        print("> Reading Alpha coefficients ...")
        results_th = pd.DataFrame(results_th)
        results_th.set_index('n', inplace=True)
        
        with os.scandir(coeff_alpha_dir) as inst_it: 
            for instance in inst_it:
                if instance.name.lower().endswith('.json'):
                    graphname = os.path.splitext(instance.name)[0].replace('_alpha', '')
                    results_alpha["n"].append(graphname)
                    filepath = os.path.join(coeff_alpha_dir, instance.name)
                    # Open the JSON file for reading
                    with open(filepath, 'r') as file:
                        # Load JSON data
                        data = json.load(file)
                        times = []
                        for _, item in data.items():
                            times.append(item[1])
                        
                        times = np.array(times)
                        results_alpha["min_alpha"].append(np.min(times))
                        results_alpha["max_alpha"].append(np.max(times))
                        results_alpha["average_alpha"].append(np.average(times))
        
        results_alpha = pd.DataFrame(results_alpha)
        results_alpha.set_index('n', inplace=True)

        results = results_th.merge(results_alpha, on='n', how='left')
        print("> Saving table: %s" % ('cpu_rounds_%s.tex' % d.lower()))
        filepath = os.path.join(tables_dir, 'cpu_coeff_%s.csv' % d.lower())
        results.round(2).to_csv(filepath)
