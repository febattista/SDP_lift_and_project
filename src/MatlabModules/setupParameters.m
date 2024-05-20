%%% Options
clear

matParams = struct();

addpath(genpath(fullfile('..', 'SDPNAL+v1.0')),path);
% SDPNAL Options
matParams.sdpnal = SDPNALplus_parameters;

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
matParams.EPSILON = Params.EPSILON;
matParams.TAILOFF = Params.TAILOFF;
matParams.MAX_CUTS_PER_ITER = double(Params.MAX_CUTS_PER_ITER);

clear Params paramDir key keys datasets i ans RESTOREDEFAULTPATH_EXECUTED

Params = matParams;

clear matParams





