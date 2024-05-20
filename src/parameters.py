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
#   ├── graphs/       : Graphs in DIMACS file format .stb
#   ├── lp/           : Linear formulations in .lp format
#                       (created during LP formulations)
#   └── models/       : SDP models in .mat files
#                       (created during SDP formulations)
datasets = {
     'smallDIMACS' : '/Users/feb223/projects/SDP_MSSP_GCP/data/StableSets/smallDIMACS',
     'smallRandom' : '/Users/feb223/projects/SDP_MSSP_GCP/data/StableSets/smallRandom'
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

# -----------------------------------------------------------------
#   MATLAB PARAMETERS
# -----------------------------------------------------------------

# Kelley cutting plane parameters
EPSILON = 1e-3           # Threshold for a cut to be considered violated
TAILOFF = 1e-1           # Tail off tolerance during cutting plane iterations
MAX_CUTS_PER_ITER = 1000 # Max num of cuts to be added in a single cutting plane iteration


