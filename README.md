# SDP Lift-and-Project software
This software package includes Python libraries for formulating Lov&aacute;sz and Schrijver's [1] Lift-and-Project semidefinite relaxations, along with MATLAB functions implementing a Cutting Plane method for solving them relying on a state-of-the-art Semidefinite Programming (SDP) solver [3]. Although this software is designed to be as general as possible, allowing the application of the Lov&aacute;sz and Schrijver's operator to any Linear Programming (LP) formulation, a strong emphasis has been placed on LP relaxations of the Maximum Stable Set Problem (MSSP).


## Features
 * Formulation of Lov&aacute;sz and Schrijver's [1] Lift-and-Project SDP relaxations
 * Formulation of a variety of LP relaxations for the MSSP
 * Solution of SDP relaxations of the MSSP via a Cutting Plane method


# Installation Guide
To install and use the software, UNIX platforms (i.e. Linux and MacOS) are strongly recommended. The easiest way to get this package on Windows is using the Windows Subsystem for Linux ([WSL](https://learn.microsoft.com/en-us/windows/wsl/install)).


## Prerequisites
- [Anaconda](https://www.anaconda.com)
- [MATLAB](https://matlab.mathworks.com)


## Optional requirements
This repo comes with all LP formulations used in computations already available. However, tools to replicate them are provided and rely on C routines `cliquer` [2]. In order to use these functionalities, the following packages are required:
- [Make](https://www.gnu.org/software/make/)
- [GCC](https://gcc.gnu.org)


## Installation Steps
1. Create `conda` environment `sdplift` (by default)
```
conda env create -f environment.yml
```
2. Activate the `conda` environment
```
conda activate sdplift
```
3. Run `install.py` to download and build dependencies [2, 3]
```
python install.py
```


# Usage
All experiments are configured in `src/parameters.py`. The key parameters are described below. Run all commands from the `src/` directory.

## Quick start (test graphs)
The `testGraphs` dataset is enabled by default and contains a small set of graphs suitable for verifying the setup.

1. **Build LP and SDP formulations** (Python):
```
cd src
python model_building.py
```
This performs three phases (each can be enabled/disabled independently in `parameters.py`):
- **Coefficient computation** (`MAKE_COEFFICIENTS`): computes theta (ADMM) and alpha (cliquer) nodal coefficients and caches them as JSON files in `data/StableSets/DATASET/coeff_theta/` and `coeff_alpha/`.
- **LP formulation** (`MAKE_LP_FORMULATION`): builds five LP relaxations per graph (NOD_gamma, NOD_theta, NOD_alpha, COV, EDGE) and writes them to `data/StableSets/DATASET/lp/`.
- **SDP model generation** (`MAKE_SDP_MODELS`): applies the M+ lift-and-project operator to each LP and archives the resulting `.mat` models in `data/StableSets/DATASET/models/models_DATASET.zip`.

2. **Solve the SDP models** (MATLAB): open a MATLAB terminal with `src/` as the working directory and run:
```
run_experiments
```
Results are written to `results/StableSets/DATASET/results_DATASET.csv` and violated-cut indices to `results/StableSets/DATASET/violated_cuts/`.

3. **Analyse results and generate LaTeX tables** (Python, optional):
```
python analyze_results.py
```
LaTeX tables and summary CSVs are written to `results/StableSets/DATASET/tables/`.


## Directory layout
```
SDP_lift_and_project/
├── data/StableSets/          input data (graphs, LP files, SDP models)
│   └── DATASET/
│       ├── graphs/           DIMACS .stb graph files (required)
│       ├── lp/               generated LP formulations
│       ├── coeff_alpha/      cached alpha coefficients (JSON)
│       ├── coeff_theta/      cached theta coefficients (JSON)
│       ├── models/           SDP models archive (.zip of .mat)
│       └── aux_data_DATASET.csv    known optimal stability numbers
│
├── results/StableSets/       experiment outputs (gitignored)
│   └── DATASET/
│       ├── results_DATASET.csv     solver output (written by MATLAB)
│       ├── violated_cuts/          added cut indices (.mat)
│       └── tables/                 LaTeX tables and summary CSVs
│
└── src/                      all source code
    ├── parameters.py         single configuration file
    ├── model_building.py     builds LP and SDP formulations
    ├── analyze_results.py    post-processes results
    ├── run_experiments.m     MATLAB solver driver
    ├── pyModules/            Python modules
    └── MatlabModules/        MATLAB modules
```


## Key parameters (`src/parameters.py`)
| Parameter | Description |
|-----------|-------------|
| `datasets` | Dict mapping dataset name to input data path |
| `results_dir` | Output directory |
| `MAKE_COEFFICIENTS` | Compute theta/alpha coefficients (skipped if JSON cache exists) |
| `MAKE_LP_FORMULATION` | Build LP relaxations |
| `MAKE_SDP_MODELS` | Apply M+ operator and generate SDP archives |
| `MAKE_TH_PLUS_MODELS` | Include Theta+ baseline SDP |
| `DO_SPLIT` | Split large `.mat` files into 5 parts (needed for very large graphs) |
| `ADMM_SOLVER_THREADS` | BLAS thread count for theta coefficient computation |
| `SDP_LIFTING_STEP` | Batch size for building sparse SDP constraint matrices |
| `CUTTING_PLANE_EPSILON` | Violation threshold to add a cut |
| `CUTTING_PLANE_TAILOFF` | Minimum improvement per cutting-plane round |
| `CUTTING_PLANE_TIMELIMIT` | Wall-clock time limit per instance (seconds) |


## Adding custom instances
Place your graph file(s) in `data/StableSets/testGraphs/graphs/` in the DIMACS edge-list format (`.stb`):
```
p edge <num_nodes> <num_edges>
e <u> <v>
...
```
Then add the known stability number for each graph to `data/StableSets/testGraphs/aux_data_testgraphs.csv`:
```
example  5
```

## Running the tests
```
conda activate sdplift
pytest tests/
```


## Datasets
A collection of 5 datasets are provided:
  * `DIMACS`: 36 graphs from the Second DIMACS Implementation Challenge benchmark
  * `smallDIMACS`: compact subset of `DIMACS` (graphs with ≤ 250 nodes, 11 instances)
  * `Random`: 315 Erdös–Rényi random graphs at various sizes and densities
  * `smallRandom`: compact subset of `Random` (200-node graphs only, 45 instances)
  * `testGraphs`: small toy graphs for testing and development


## References
 1. L&aacute;szl&oacute; Lov&aacute;sz and Alexander Schrijver. *Cones of matrices and set-functions and 0–1 optimization*. SIAM journal on optimization, 1(2):166–190, 1991.
 2. Sampo Niskanen and Patric R. J. Östergård, *Cliquer User's Guide, Version 1.0*, Communications Laboratory, Helsinki University of Technology, Espoo, Finland, Tech. Rep. T48, 2003.
 3. Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, and Xin-Yuan Zhao. *SDPNAL+: A Matlab software for semidefinite programming with bound constraints*. Optimization Methods and Software, 35(1):87–115, 2020.
