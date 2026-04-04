# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Post-process experimental results and generate LaTeX tables.

import os, re, collections, json
import numpy as np
import pandas as pd
from scipy.io import loadmat

from parameters import datasets, results_datasets
from pyModules.SDPLifting import map_nod_constraints, map_clique_constraints, map_edge_constraints

pd.options.mode.copy_on_write = True


def compute_gaps(dfin, columns, target='alpha'):
    """Add percentage-gap columns for each column relative to target."""
    for column in columns:
        dfin['gap' + column] = 100 * (dfin[column] - dfin[target]) / dfin[target]
    return dfin


def count_values(input_dict, idx):
    """Count how many added cuts belong to each constraint class."""
    output_dict = {}
    for i in idx:
        value = input_dict[i[0]]
        output_dict[value] = output_dict.get(value, 0) + 1
    return output_dict


def get_idx_formatter(dataset):
    """Return a pandas style formatter for the index columns."""
    if 'random' in dataset.lower():
        return {'n': '{:.0f}', 'd': '{:.1f}'}
    return {'n': '{:s}'}


def get_col_formatter(columns):
    """Return a pandas style formatter for result columns."""
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


def make_index(df, dataset):
    """Set the DataFrame index based on dataset type."""
    if 'random' in dataset.lower():
        df['n'] = df['n'].str.replace('G_', '')
        df[['n', 'd', 'i']] = df['n'].str.split('_', expand=True).astype(float)
        df.set_index(['n', 'd', 'i'], inplace=True)
    else:
        df.set_index(['n'], inplace=True)


def escape_special_char(tex_file_path):
    """Escape underscores in a LaTeX file."""
    if os.path.exists(tex_file_path):
        with open(tex_file_path, 'r') as f:
            lines = [line.replace('_', r'\_') for line in f]
        with open(tex_file_path, 'w') as f:
            f.writelines(lines)


# Column groups
columns   = ['TH+', 'NOD1_obj', 'NOD2_obj', 'NOD3_obj', 'COV_obj', 'EDGE_obj']
gaps      = ['gapTH+', 'gapNOD1_obj', 'gapNOD2_obj', 'gapNOD3_obj', 'gapCOV_obj', 'gapEDGE_obj']
cuts      = ['NOD1_cut', 'NOD2_cut', 'NOD3_cut', 'COV_cut', 'EDGE_cut']
rounds    = ['NOD1_rounds', 'NOD2_rounds', 'NOD3_rounds', 'COV_rounds', 'EDGE_rounds']
cpu_times = ['TH+_time', 'NOD1_time', 'NOD2_time', 'NOD3_time', 'COV_time', 'EDGE_time']

all_cols = ['TH+', 'gapTH+', 'TH+_time',
            'NOD1_obj', 'gapNOD1_obj', 'NOD1_rounds', 'NOD1_time',
            'NOD2_obj', 'gapNOD2_obj', 'NOD2_rounds', 'NOD2_time',
            'NOD3_obj', 'gapNOD3_obj', 'NOD3_rounds', 'NOD3_time',
            'COV_obj',  'gapCOV_obj',  'COV_rounds',  'COV_time',
            'EDGE_obj', 'gapEDGE_obj', 'EDGE_rounds', 'EDGE_time']

ub_gaps = ['TH+', 'gapTH+',
           'EDGE_obj', 'gapEDGE_obj',
           'COV_obj',  'gapCOV_obj',
           'NOD1_obj', 'gapNOD1_obj',
           'NOD2_obj', 'gapNOD2_obj',
           'NOD3_obj', 'gapNOD3_obj']

cpu_rounds = ['TH+_time',
              'EDGE_rounds', 'EDGE_time',
              'COV_rounds',  'COV_time',
              'NOD1_rounds', 'NOD1_time',
              'NOD2_rounds', 'NOD2_time',
              'NOD3_rounds', 'NOD3_time']

cols_ord = [
    'edge_bound2', 'edge_bound4', 'edge_lift_edge1:1', 'edge_lift_edge1:2',
    'cov_bound2',  'cov_bound4',  'cov_lift_clq1:1',   'cov_lift_clq2:1',
    'cov_lift_clq1:2', 'cov_lift_clq2:2',
    'nod_gamma_bound2', 'nod_gamma_bound4',
    'nod_gamma_lift_nod1:1', 'nod_gamma_lift_nod2:1', 'nod_gamma_lift_nod3:1',
    'nod_gamma_lift_nod1:2', 'nod_gamma_lift_nod2:2', 'nod_gamma_lift_nod3:2',
    'nod_theta_bound2', 'nod_theta_bound4',
    'nod_theta_lift_nod1:1', 'nod_theta_lift_nod2:1', 'nod_theta_lift_nod3:1',
    'nod_theta_lift_nod1:2', 'nod_theta_lift_nod2:2', 'nod_theta_lift_nod3:2',
    'nod_alpha_bound2', 'nod_alpha_bound4',
    'nod_alpha_lift_nod1:1', 'nod_alpha_lift_nod2:1', 'nod_alpha_lift_nod3:1',
    'nod_alpha_lift_nod1:2', 'nod_alpha_lift_nod2:2', 'nod_alpha_lift_nod3:2',
]


for d in datasets:
    print("--------------------------------------------------")
    print("  DATASET %s " % d)
    print("--------------------------------------------------")

    data_path    = datasets[d]
    results_path = results_datasets[d]

    graph_dir       = os.path.join(data_path, 'graphs')
    lp_dir          = os.path.join(data_path, 'lp')
    coeff_alpha_dir = os.path.join(data_path, 'coeff_alpha')
    coeff_theta_dir = os.path.join(data_path, 'coeff_theta')

    viol_cuts_dir = os.path.join(results_path, 'violated_cuts')
    tables_dir    = os.path.join(results_path, 'tables')
    res_csv       = "results_%s.csv" % d.lower()
    aux_csv       = "aux_data_%s.csv" % d.lower()

    if not os.path.exists(data_path):
        print("WARNING: path to Dataset %s does not exist! Skipping ..." % d)
        continue

    if not os.path.exists(os.path.join(results_path, res_csv)):
        print("WARNING: %s not found in %s! Skipping ..." % (res_csv, results_path))
        continue

    os.makedirs(tables_dir, exist_ok=True)

    print("> Reading: %s" % res_csv)
    with open(os.path.join(results_path, res_csv)) as f:
        l = []
        for line in f:
            l.append(re.sub(r'\s+', '\t', line).split('\t')[:-1])
    df = pd.DataFrame(l[1:], columns=l[0])

    print("> Reading: %s" % aux_csv)
    with open(os.path.join(data_path, aux_csv)) as f:
        l = []
        for line in f:
            l.append(re.sub(r'\s+', '\t', line).split('\t')[:-1])
    aux = pd.DataFrame(l[1:], columns=l[0])

    df = df.merge(aux, on='n', how='left')

    target_cols = list(set(df.columns) - set('n'))
    df[target_cols] = df[target_cols].astype(float)

    compute_gaps(df, (['NOD_obj'] if 'random' in d.lower() else []) + columns)

    print("> Saving table: %s" % ('raw_results_%s.csv' % d.lower()))
    target = df[list(aux.columns) + (['gapNOD_obj'] if 'random' in d.lower() else []) + all_cols]
    make_index(target, d)
    target.to_csv(os.path.join(tables_dir, 'raw_results_%s.csv' % d.lower()))

    formatter_col = get_col_formatter(list(target.columns))
    formatter_idx = get_idx_formatter(d)

    print("> Saving table: %s" % ('ub_gaps_%s.tex' % d.lower()))
    filepath = os.path.join(tables_dir, 'ub_gaps_%s.tex' % d.lower())
    if 'random' in d.lower():
        target[['gapNOD_obj'] + ub_gaps].groupby(by=['n', 'd']).mean() \
            .style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .highlight_min(subset=gaps, axis=1, props="textbf:--rwrap;") \
            .highlight_min(subset=columns[1:], axis=1, props="textbf:--rwrap;") \
            .to_latex(buf=filepath)
        escape_special_char(filepath)
    else:
        target[ub_gaps].style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .highlight_min(subset=gaps, axis=1, props="textbf:--rwrap;") \
            .highlight_min(subset=columns[1:], axis=1, props="textbf:--rwrap;") \
            .to_latex(buf=filepath)
        escape_special_char(filepath)

    print("> Saving table: %s" % ('cpu_rounds_%s.tex' % d.lower()))
    filepath = os.path.join(tables_dir, 'cpu_rounds_%s.tex' % d.lower())
    if 'random' in d.lower():
        target[cpu_rounds].groupby(by=['n', 'd']).mean() \
            .style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .to_latex(buf=filepath)
        escape_special_char(filepath)
    else:
        target[cpu_rounds].style.format(formatter=formatter_col) \
            .format_index(formatter=formatter_idx) \
            .to_latex(buf=filepath)
        escape_special_char(filepath)

    if os.path.exists(coeff_theta_dir) and os.path.exists(coeff_alpha_dir):
        print("> Reading Theta coefficients ...")
        results_th    = collections.defaultdict(list)
        results_alpha = collections.defaultdict(list)

        with os.scandir(coeff_theta_dir) as inst_it:
            for instance in inst_it:
                if instance.name.lower().endswith('.json'):
                    graphname = os.path.splitext(instance.name)[0].replace('_theta', '')
                    results_th["n"].append(graphname)
                    with open(os.path.join(coeff_theta_dir, instance.name)) as f:
                        data = json.load(f)
                        times = np.array([item[1] for item in data.values()])
                        results_th["min_theta"].append(np.min(times))
                        results_th["max_theta"].append(np.max(times))
                        results_th["average_theta"].append(np.average(times))

        results_th = pd.DataFrame(results_th).set_index('n') if results_th["n"] else pd.DataFrame()

        print("> Reading Alpha coefficients ...")
        with os.scandir(coeff_alpha_dir) as inst_it:
            for instance in inst_it:
                if instance.name.lower().endswith('.json'):
                    graphname = os.path.splitext(instance.name)[0].replace('_alpha', '')
                    results_alpha["n"].append(graphname)
                    with open(os.path.join(coeff_alpha_dir, instance.name)) as f:
                        data = json.load(f)
                        times = np.array([item[1] for item in data.values()])
                        results_alpha["min_alpha"].append(np.min(times))
                        results_alpha["max_alpha"].append(np.max(times))
                        results_alpha["average_alpha"].append(np.average(times))

        results_alpha = pd.DataFrame(results_alpha).set_index('n') if results_alpha["n"] else pd.DataFrame()

        if not results_th.empty and not results_alpha.empty:
            results = results_th.merge(results_alpha, on='n', how='left')
            print("> Saving table: %s" % ('cpu_coeff_%s.csv' % d.lower()))
            results.round(2).to_csv(os.path.join(tables_dir, 'cpu_coeff_%s.csv' % d.lower()))
        else:
            print("> Skipping table: Alpha and/or Theta coefficient datasets are empty.")

    print("> Reading Violated Cuts ...")
    results = collections.defaultdict(list)
    relaxations = ['edge', 'cov', 'nod_alpha', 'nod_theta', 'nod_gamma']

    if os.path.exists(viol_cuts_dir):
        with os.scandir(graph_dir) as inst_it:
            for instance in inst_it:
                graphname = os.path.splitext(instance.name)[0]
                results["n"].append(graphname)

                for rel in relaxations:
                    filename  = '%s_%s_viol_test.mat' % (graphname, rel)
                    lpname    = graphname + '_' + rel + '.lp'
                    filepath  = os.path.join(viol_cuts_dir, filename)
                    lppath    = os.path.join(lp_dir, lpname)

                    if os.path.exists(filepath) and os.path.exists(lppath):
                        if 'nod' in rel:
                            maps = map_nod_constraints(lppath)
                        elif 'cov' in rel:
                            maps = map_clique_constraints(lppath)
                        elif 'edge' in rel:
                            maps = map_edge_constraints(lppath)
                        else:
                            print('Unknown relaxation: %s' % rel)
                            continue

                        h = loadmat(filepath)
                        added_cuts = count_values(maps, h['added_cuts_idx'])
                        for i in set(maps.values()):
                            results[rel + '_' + i].append(added_cuts.get(i, 0))

        df_viol_cuts = pd.DataFrame(results)
        make_index(df_viol_cuts, d)
        df_viol_cuts = df_viol_cuts.sort_index()
        df_viol_cuts = df_viol_cuts[cols_ord]

        filepath = os.path.join(tables_dir, 'raw_violated_cuts_%s.csv' % d.lower())
        print("> Saving table: %s" % ('raw_violated_cuts_%s.csv' % d.lower()))
        df_viol_cuts.to_csv(filepath)

        filepath = os.path.join(tables_dir, 'violated_cuts_%s.tex' % d.lower())
        print("> Saving table: %s" % ('violated_cuts_%s.tex' % d.lower()))
        columns_to_remove = df_viol_cuts.columns[df_viol_cuts.eq(0).all()]
        df_viol_cuts = df_viol_cuts.loc[:, ~df_viol_cuts.columns.isin(columns_to_remove)]

        if 'random' in d.lower():
            df_viol_cuts.groupby(by=['n', 'd']).mean() \
                .style.format(precision=0) \
                .format_index(formatter=formatter_idx) \
                .to_latex(buf=filepath)
            escape_special_char(filepath)
        else:
            df_viol_cuts \
                .style.format(precision=0) \
                .format_index(formatter=formatter_idx) \
                .to_latex(buf=filepath)
            escape_special_char(filepath)
