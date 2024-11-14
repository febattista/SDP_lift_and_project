# Script to set up parameters for computational results

# -----------------------------------------------------------------
#   DATASET PATHS
# -----------------------------------------------------------------
# e.g. datasets = {'dataset_name' : 'dataset_path'}
# Each dataset dir is assumed to have the following structure
#   path/to/dataset/
#   ├── coeff_alpha/  : Alpha coefficients in Nodal LP formulation
#                       (created during LP formulations)
#   ├── coeff_theta/  : Theta coefficients in Nodal LP formulation
#                       (created during LP formulations)
#   ├── graphs/       : Graphs in DIMACS file format .stb (Required)
#   ├── lp/           : Linear formulations in .lp format
#                       (created during LP formulations)
#   ├── tables/       : Results and LaTeX tables
#                       (created by analyze_results.py)
#   └── models/       : SDP models in .mat files
#                       (created during SDP formulations)
datasets = {
     'testGraphs'  : '../data/StableSets/testGraphs',
     'smallDIMACS' : '../data/StableSets/smallDIMACS',
     'smallRandom' : '../data/StableSets/smallRandom',
          # 'DIMACS' : '../data/StableSets/DIMACS',
          # 'Random' : '../data/StableSets/Random'
}

# -----------------------------------------------------------------
#   PYTHON PARAMETERS
# -----------------------------------------------------------------

# Create LP formulations?
MAKE_LP_FORMULATION = False

# Create SDP models?
MAKE_SDP_MODELS = True

# Create Theta_plus models?
MAKE_TH_PLUS_MODELS = True

# Splitting .mat files in 5 parts?
# scipy savemat has a limit on the size of .mat files it can save
# MATLAB also uses this parameter
DO_SPLIT = True

# ADMMsolver.py parameters
# Used only if MAKE_LP_FORMULATION == True
# in NOD_theta LP formulation
ADMM_SOLVER_TOL = 1e-4
ADMM_SOLVER_MAX_ITER = 1000000
ADMM_SOLVER_TIMELIMIT = 3600
ADMM_SOLVER_DEBUG = False

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

