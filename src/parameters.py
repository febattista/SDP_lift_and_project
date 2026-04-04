# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Configuration file for experiments: datasets, paths, and solver parameters.

# -----------------------------------------------------------------
#   DATASET PATHS
# -----------------------------------------------------------------
# e.g. datasets = {'dataset_name' : 'dataset_path'}
# Input data for each dataset. Expected subdirectory structure:
#   path/to/dataset/
#   ├── coeff_alpha/  : Alpha coefficients for Nodal LP formulation
#                       (created during coefficient computation)
#   ├── coeff_theta/  : Theta coefficients for Nodal LP formulation
#                       (created during coefficient computation)
#   ├── graphs/       : Graphs in DIMACS file format .stb (Required)
#   ├── lp/           : Linear formulations in .lp format
#                       (created during LP formulations)
#   └── aux_data_DATASETNAME.csv : Known optimal stability numbers
datasets = {
     'testGraphs'  : '../data/StableSets/testGraphs',
     'smallDIMACS' : '../data/StableSets/smallDIMACS',
     # 'smallRandom' : '../data/StableSets/smallRandom',
          # 'DIMACS' : '../data/StableSets/DIMACS',
          # 'Random' : '../data/StableSets/Random'
}

# -----------------------------------------------------------------
#   RESULTS DIRECTORY
# -----------------------------------------------------------------
# Base directory for all generated outputs: SDP model archives,
# violated cuts, result CSVs, and LaTeX tables.
# The subdirectory layout mirrors that of data/.
# Output subdirectory structure per dataset:
#   results_dir/StableSets/DATASET/
#   ├── models/        : SDP model archives (.zip of .mat files)
#   ├── violated_cuts/ : Cut indices written by the cutting plane (MATLAB)
#   └── tables/        : LaTeX tables and summary CSVs (analyze_results.py)
results_dir = '../results'

# Auto-derived: replaces the '../data' prefix in each dataset path with results_dir.
results_datasets = {k: v.replace('../data', results_dir, 1) for k, v in datasets.items()}

# -----------------------------------------------------------------
#   PYTHON PARAMETERS
# -----------------------------------------------------------------

# Compute LP coefficients (theta, alpha)?
# Set to False to skip and reuse existing .json cache files.
MAKE_COEFFICIENTS = True

# Create LP formulations?
MAKE_LP_FORMULATION = True

# Clique cover algorithm used in QSTABC:
#   True  — greedy_clique_cover_letchford_et_al (tighter cover, slower)
#   False — greedy_clique_cover (faster, looser cover)
USE_CLIQUE_COVER_FROM_LETCHFORD_ET_AL = True

# Create SDP models?
MAKE_SDP_MODELS = True

# Create Theta_plus models?
MAKE_TH_PLUS_MODELS = True

# Splitting .mat files in 5 parts?
# scipy savemat has a limit on the size of .mat files it can save.
# MATLAB also uses this parameter.
DO_SPLIT = False

# ADMMsolver.py parameters
# Used only if MAKE_COEFFICIENTS == True (for NOD_theta coefficient computation)
ADMM_SOLVER_TOL = 1e-4
ADMM_SOLVER_MAX_ITER = 1000000
ADMM_SOLVER_TIMELIMIT = 3600
ADMM_SOLVER_DEBUG = False
ADMM_SOLVER_THREADS = 20       # Number of threads for parallel BLAS operations

# SDPLifting.py parameters
SDP_LIFTING_STEP = 100000      # Batch size for building the sparse inequality
                                # constraint matrix in m_plus_lifting

# -----------------------------------------------------------------
#   MATLAB PARAMETERS
# -----------------------------------------------------------------

# Kelley cutting plane parameters
CUTTING_PLANE_EPSILON = 1e-3            # Threshold for a cut to be
                                        # considered violated

CUTTING_PLANE_TAILOFF = 1e-1            # Tail off tolerance during
                                        # cutting plane iterations

CUTTING_PLANE_MAX_CUTS_PER_ITER = 1000  # Max num of cuts to be added in
                                        # a single cutting plane iteration

CUTTING_PLANE_TIMELIMIT = 7200          # Time limit in seconds


# SDPNALplus_parameters
# (see ./ThirdParty/SDPNAL+v1.0/SDPNALplus_parameters.m)
SDPNAL_TOL = 1e-6
SDPNAL_MAXITER = 2000000
SDPNAL_MAXTIME = 250000
SDPNAL_ADM_MAXITER = 20000
SDPNAL_STOP_OPTION = 0
SDPNAL_PRINT_LEVEL = 0
