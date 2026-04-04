% Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
% SPDX-License-Identifier: GPL-3.0-or-later
% Main MATLAB script: solves SDP models via the Kelley cutting plane method
% using SDPNAL+ for each dataset specified in parameters.py.

addpath 'MatlabModules';
%%% Read parameters from parameters.py
setupParameters;

datasets        = fieldnames(Params.datasets);
results_datasets = fieldnames(Params.results_datasets);

% For each dataset solve all models we can find
for d=1:length(datasets)
    % Set up paths
    dataset      = char(datasets(d));
    data_path    = Params.datasets.(dataset);
    results_path = Params.results_datasets.(dataset);
    model_dir     = fullfile(results_path, "models");
    graph_dir     = fullfile(data_path,   "graphs");
    viol_cuts_dir = fullfile(results_path, "violated_cuts");

    if ~exist(data_path,  'dir') || ~exist(graph_dir, 'dir') ...
    || ~exist(model_dir, 'dir')
        fprintf('WARNING: path to Dataset %s does not exist!\n', dataset);
        continue
    end

    if ~exist(results_path, 'dir')
        mkdir(results_path);
    end

    if ~exist(viol_cuts_dir, 'dir')
        mkdir(viol_cuts_dir);
    end

    % SDP Models zip file
    zipFilename = sprintf("models_%s.zip", lower(dataset));

    if ~exist(fullfile(model_dir, zipFilename), 'file')
        fprintf('WARNING: Cannot find %s archive in %s.\n', zipFilename, model_dir);
        continue
    end

    zipJavaFile = java.io.File(fullfile(model_dir, zipFilename));
    zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);

    files = dir(fullfile(graph_dir, '*.stb'));

    file_h = ['n\tTH+\tTH+_time\t' ...
              'NOD1_obj\tNOD1_cut\tNOD1_time\tNOD1_rounds\t' ...
              'NOD2_obj\tNOD2_cut\tNOD2_time\tNOD2_rounds\t' ...
              'NOD3_obj\tNOD3_cut\tNOD3_time\tNOD3_rounds\t' ...
              'COV_obj\tCOV_cut\tCOV_time\tCOV_rounds\t' ...
              'EDGE_obj\tEDGE_cut\tEDGE_time\tEDGE_rounds\n'];

    file_l = ['%s\t%10.6f\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\t%d\t' ...
              '%10.6f\t%d\t%10.2f\t%d\t' ...
              '%10.6f\t%d\t%10.2f\t%d\t' ...
              '%10.6f\t%d\t%10.2f\t%d\t' ...
              '%10.6f\t%d\t%10.2f\t%d\n'];

    res_csv = fullfile(results_path, sprintf('results_%s.csv', lower(dataset)));
    f = fopen(res_csv, 'w');
    fprintf(f, file_h);

    for i = 1:length(files)
        file = files(i).name(1:end-4);
        fprintf('Processing Graph: %s\n', file);
        obj_th = [-Inf, -Inf]; theta_plus_time = 0.0;
        nod1_obj = -Inf; nod1_cuts = 0; nod1_time = 0.0; nod1_rounds = 0;
        nod2_obj = -Inf; nod2_cuts = 0; nod2_time = 0.0; nod2_rounds = 0;
        nod3_obj = -Inf; nod3_cuts = 0; nod3_time = 0.0; nod3_rounds = 0;
        cov_obj  = -Inf; cov_cuts  = 0; cov_time  = 0.0; cov_rounds  = 0;
        edge_obj = -Inf; edge_cuts = 0; edge_time = 0.0; edge_rounds = 0;

        % Extract and solve Theta+
        model = fullfile(strcat(file, '_th+.mat'));
        termcode = extract_file(zipFile, model, fullfile(model_dir, model), false);
        if ~termcode
            fprintf('WARNING: Cannot find Theta plus model for %s in %s.\n', file, zipFilename);
        else
            model_path = fullfile(model_dir, model);
            load(model_path);
            At_th = At;
            b_th  = b;
            blk   = {'s', s};
            [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = ...
                sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], Params.sdpnal);
            theta_plus_time = info.totaltime;

            clean_files(model_dir, model(1:end-4), false);

            % Nodal_gamma
            model = fullfile(strcat(file, '_nod_gamma.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT);
            if termcode
                [nod1_obj, nod1_cuts, nod1_time, nod1_rounds] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_gamma', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Nodal_theta
            model = fullfile(strcat(file, '_nod_theta.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT);
            if termcode
                [nod2_obj, nod2_cuts, nod2_time, nod2_rounds] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_theta', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Nodal_alpha
            model = fullfile(strcat(file, '_nod_alpha.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT);
            if termcode
                [nod3_obj, nod3_cuts, nod3_time, nod3_rounds] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_alpha', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Greedy clique cover
            model = fullfile(strcat(file, '_cov.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT);
            if termcode
                [cov_obj, cov_cuts, cov_time, cov_rounds] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_cov', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Edge formulation
            model = fullfile(strcat(file, '_edge.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT);
            if termcode
                [edge_obj, edge_cuts, edge_time, edge_rounds] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_edge', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

        end

        fprintf(f, file_l, file, ...
            -obj_th(1), theta_plus_time, ...
            -nod1_obj, nod1_cuts, nod1_time, nod1_rounds, ...
            -nod2_obj, nod2_cuts, nod2_time, nod2_rounds, ...
            -nod3_obj, nod3_cuts, nod3_time, nod3_rounds, ...
            -cov_obj,  cov_cuts,  cov_time,  cov_rounds, ...
            -edge_obj, edge_cuts, edge_time, edge_rounds);
    end

    zipFile.close();

end
