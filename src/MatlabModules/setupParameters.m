% Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
% SPDX-License-Identifier: GPL-3.0-or-later
% Load experiment parameters from parameters.py into the MATLAB workspace.

clear

matParams = struct();

addpath(genpath(fullfile('..', 'ThirdParty', 'SDPNAL+v1.0')), path);

% Import the parameters from Python module
Params = py.importlib.import_module('parameters');
py.importlib.reload(Params);

% Convert Python dicts to MATLAB structs (string values only)
tmp = struct(Params.datasets);
keys = fieldnames(tmp);
for i = 1:length(keys)
    tmp.(keys{i}) = char(tmp.(keys{i}));
end
matParams.datasets = tmp;

tmp = struct(Params.results_datasets);
keys = fieldnames(tmp);
for i = 1:length(keys)
    tmp.(keys{i}) = char(tmp.(keys{i}));
end
matParams.results_datasets = tmp;

matParams.DO_SPLIT = Params.DO_SPLIT;

% Cutting plane parameters
matParams.CUTTING_PLANE_EPSILON           = Params.CUTTING_PLANE_EPSILON;
matParams.CUTTING_PLANE_TAILOFF           = Params.CUTTING_PLANE_TAILOFF;
matParams.CUTTING_PLANE_MAX_CUTS_PER_ITER = double(Params.CUTTING_PLANE_MAX_CUTS_PER_ITER);
matParams.CUTTING_PLANE_TIMELIMIT         = double(Params.CUTTING_PLANE_TIMELIMIT);

% SDPNAL+ options
matParams.sdpnal            = SDPNALplus_parameters;
matParams.sdpnal.tol        = Params.SDPNAL_TOL;
matParams.sdpnal.maxiter    = double(Params.SDPNAL_MAXITER);
matParams.sdpnal.maxtime    = double(Params.SDPNAL_MAXTIME);
matParams.sdpnal.ADMmaxiter = double(Params.SDPNAL_ADM_MAXITER);
matParams.sdpnal.stopoption = double(Params.SDPNAL_STOP_OPTION);
matParams.sdpnal.printlevel = double(Params.SDPNAL_PRINT_LEVEL);

clear Params tmp keys i ans RESTOREDEFAULTPATH_EXECUTED

Params = matParams;
clear matParams
