%%% Options
clear

matParams = struct();

addpath(genpath(fullfile('..', 'ThirdParty', 'SDPNAL+v1.0')),path);

% Import the parameters from Python module
Params = py.importlib.import_module('parameters');
py.importlib.reload(Params);

% Convert the Python dictionary to a MATLAB struct
datasets = struct(Params.datasets);
keys = fieldnames(datasets);
for i = 1:length(keys)
    key = keys{i};
    datasets.(key) = char(datasets.(key));
end

matParams.datasets = datasets;
matParams.DO_SPLIT = Params.DO_SPLIT;

% CUTTING PLANE PARAMETERS
matParams.CUTTING_PLANE_EPSILON = Params.CUTTING_PLANE_EPSILON;
matParams.CUTTING_PLANE_TAILOFF = Params.CUTTING_PLANE_TAILOFF;
matParams.CUTTING_PLANE_MAX_CUTS_PER_ITER = ...
        double(Params.CUTTING_PLANE_MAX_CUTS_PER_ITER);
matParams.CUTTING_PLANE_TIMELIMIT = double(Params.CUTTING_PLANE_TIMELIMIT);

% SDPNAL Options
matParams.sdpnal = SDPNALplus_parameters;
matParams.sdpnal.tol = Params.SDPNAL_TOL;
matParams.sdpnal.maxiter = double(Params.SDPNAL_MAXITER);
matParams.sdpnal.maxtime = double(Params.SDPNAL_MAXTIME);
matParams.sdpnal.ADMmaxiter = double(Params.SDPNAL_ADM_MAXITER);
matParams.sdpnal.stopoption = double(Params.SDPNAL_STOP_OPTION);
matParams.sdpnal.printlevel = double(Params.SDPNAL_PRINT_LEVEL);

clear Params paramDir key keys datasets i ans RESTOREDEFAULTPATH_EXECUTED

Params = matParams;

clear matParams





