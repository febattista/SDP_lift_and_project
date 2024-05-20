
addpath 'MatlabModules';
display(path);
%%% Read parameters from parameters.py
setupParameters;

datasets = fieldnames(Params.datasets);

% For each dataset and
% For each '.stb' graphs ...
% Solve all models we can find 
for d=1:length(datasets)
    % Set up paths
    dataset = char(datasets(d));
    data_path = Params.datasets.(dataset);
    model_dir = fullfile(data_path, "models");
    graph_dir = fullfile(data_path, "graphs");
    viol_cuts_dir = fullfile(data_path, "violated_cuts");

    if ~exist(data_path, 'dir') || ~exist(model_dir, 'dir') ... 
    || ~exist(model_dir, 'dir') || ~exist(graph_dir, 'dir')
        fprintf('WARNING: path to Dataset %s does not exist!\n', dataset);
        continue
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

    num_pieces = 5;

    files = dir(fullfile(graph_dir, '*.stb'));

    file_h = ['n\tTH+\tTH+_time\t' ...
              'NOD1_obj\tNOD1_cut\tNOD1_time\t' ...
              'NOD2_obj\tNOD2_cut\tNOD2_time\t' ...
              'NOD3_obj\tNOD3_cut\tNOD3_time\t' ...
              'COV_obj\tCOV_cut\tCOV_time\t' ...
              'EDGE_obj\tEDGE_cut\tEDGE_time\n'];
    
    file_l = ['%s\t%10.6f\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\t' ...
              '%10.6f\t%d\t%10.2f\n'];
    
    f = fopen(fullfile(data_path, sprintf('results_%s.csv', lower(dataset))), 'w');
    fprintf(f, file_h);

    for i = 1:length(files)
        file = files(i).name(1:end-4);
        fprintf('Processing Graph: %s\n', file);
        % Initialize values
        obj_th = [-Inf, -Inf]; theta_plus_time = 0.0; 
        nod1_obj = -Inf; nod1_cuts = 0; nod1_time = 0.0;
        nod2_obj = -Inf; nod2_cuts = 0; nod2_time = 0.0;
        nod3_obj = -Inf; nod3_cuts = 0; nod3_time = 0.0;
        cov_obj = -Inf; cov_cuts = 0; cov_time = 0.0;
        edge_obj = -Inf; edge_cuts = 0; edge_time = 0.0;

        % Extract SDP Theta plus
        model = fullfile(strcat(file, '_th+.mat'));
        termcode = extract_file(zipFile, model, fullfile(model_dir, model), false); 
        if ~termcode
            % Theta model not found
            fprintf('WARNING: Cannot find Theta plus model for %s in %s.\n', file, zipFilename);
        else
            % Solve Theta plus
            model_path = fullfile(model_dir, model);
            load(model_path);
            % Load theta model
            At_th = At;
            b_th = b;
            blk = {'s', s};
            % Solve theta model (i.e. no cuts)
            [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = ...
                sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], Params.sdpnal);
            theta_plus_time = info.totaltime;

            clean_files(model_dir, model(1:end-4), false);

            % Now solve other SDP models using the cutting plane
            % Nodal_gamma
            model = fullfile(strcat(file, '_nod_gamma.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT); 
            if termcode
                [nod1_obj, nod1_cuts, nod1_time] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_gamma', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Nodal_theta
            model = fullfile(strcat(file, '_nod_theta.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT); 
            if termcode
                [nod2_obj, nod2_cuts, nod2_time] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_theta', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Nodal_alpha
            model = fullfile(strcat(file, '_nod_alpha.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT); 
            if termcode
                [nod3_obj, nod3_cuts, nod3_time] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_nod_alpha', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Greedy clique cover
            model = fullfile(strcat(file, '_cov.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT); 
            if termcode
                [cov_obj, cov_cuts, cov_time] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_cov', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end

            % Edge formulation
            model = fullfile(strcat(file, '_edge.mat'));
            termcode = extract_file(zipFile, model, fullfile(model_dir, model), Params.DO_SPLIT); 
            if termcode
                [edge_obj, edge_cuts, edge_time] = ...
                    kelley_cutting_plane(model_dir, viol_cuts_dir, ...
                    file, '_edge', At_th, b_th, X_th, blk, ...
                    obj_th, theta_plus_time, Params);
                
                clean_files(model_dir, model(1:end-4), Params.DO_SPLIT);
            end
            
        end
        
        % Print in the CSV file
        fprintf(f, file_l, file, ...
        -obj_th(1), theta_plus_time, ...
        -nod1_obj, nod1_cuts, nod1_time, ...
        -nod2_obj, nod2_cuts, nod2_time, ...
        -nod3_obj, nod3_cuts, nod3_time, ...
        -cov_obj, cov_cuts, cov_time, ...
        -edge_obj, edge_cuts, edge_time);
        
    end

    zipFile.close();

end
