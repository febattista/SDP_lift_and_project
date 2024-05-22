# SDP Lift-and-Project software

This software package includes Python libraries for formulating Lov&aacute;sz and Schrijver's [1] Lift-and-Project semidefinite relaxations, along with MATLAB functions implementing a Cutting Plane method for solving them relying on a state-of-the-art Semidefinite Programming (SDP) solver [3]. Although this software is designed to be as general as possible, allowing the application of the Lov&aacute;sz and Schrijver's operator to any Linear Programming (LP) formulation, a strong emphasisy has been places on LP relaxations of the Maximum Stable Set Problem (MSSP).


## Features

 * Formulation of Lov&aacute;sz and Schrijver's [1] Lift-and-Project SDP relaxations
 * Formulation of a variety of LP relaxations for the MSSP
 * Solution of SDP relaxations of the MSSP via a Cutting Plane method


# Installation Guide

To install and use the software, UNIX platforms (i.e. Linux and MacOS) are strongly recommended. The easiest way to get this package on Windows is using the Windows Subsystem for Linux ([WSL](https://learn.microsoft.com/en-us/windows/wsl/install)).

## Prerequisites
- [Anaconda](https://www.anaconda.com)
- [MATLAB](https://matlab.mathworks.com)
- [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses/) License (free for academics)

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

## Usage
Instructions on how to configure the project can be found in `parameters.py`.
1. Create LP and/or SDP formulations in Python
```
cd src/
python model_building.py
```
2. Solve them in MATLAB
```
# Current folder in MATLAB should be src
run_experiments
```
3. Collect and analyze results by creating `LaTeX` tables in Python
```
python analyze_results.py
```


## References
 1. L&aacute;szl&oacute; Lov&aacute;sz and Alexander Schrijver. *Cones of matrices and set-functions and 0–1 optimization*. SIAM journal on optimization, 1(2):166–190, 1991.
 2. Sampo Niskanen and Patric R. J. Östergård, *Cliquer User's Guide, Version 1.0*, Communications Laboratory, Helsinki University of Technology, Espoo, Finland, Tech. Rep. T48, 2003.
 3. Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, and Xin-Yuan Zhao. *SDPNAL+: A Matlab software for semidefinite programming with bound constraints*. Optimization Methods and Software, 35(1):87–115, 2020.
