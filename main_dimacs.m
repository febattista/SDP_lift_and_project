%%% Options
clear
restoredefaultpath;
addpath(genpath(strcat(pwd, '/SDPNAL+v1.0')),path);
OPTIONS = SDPNALplus_parameters;

work_dir = fullfile('StableSets', 'DIMACS');
model_folder = fullfile('StableSets', 'DIMACS', 'models');
model_th_folder = fullfile('StableSets', 'DIMACS', 'models');
sol_folder = fullfile('StableSets', 'DIMACS', 'violated_cuts');

if exist(sol_folder, 'dir') == 0
    mkdir(sol_folder);
end

zipFilename = 'models_dimacs.zip';
zipJavaFile  = java.io.File(fullfile(model_folder, zipFilename));
zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);

num_pieces = 5;


graphs = {'brock200_1'  , 'brock200_2', 'brock200_3', 'brock200_4', ...
        'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', ...
        'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', ...
        'C125-9', 'C250-9', 'C500-9', ...
        'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC500-5', ...
        'MANN_a9', 'MANN_a27', 'johnson32-2-4', ...
        'keller4',  ...
        'p_hat300-1', 'p_hat300-2', 'p_hat300-3', ...
        'p_hat500-1', 'p_hat500-2', 'p_hat500-3', ...
        'p_hat700-1', 'p_hat700-2', 'p_hat700-3', ...
        'sanr200_0.7', 'sanr200_0.9',  ...
        'sanr400_0.5', 'sanr400_0.7'};

% log_file = 'matlab_cut_selection_DIMACS.log';
% if exist(log_file, 'file')
%     delete(log_file);
% end
% diary(log_file);

file_h = 'n\tTH+\tTH+_time\tNOD1_obj\tNOD1_cut\tNOD1_time\tNOD2_obj\tNOD2_cut\tNOD2_time\tNOD3_obj\tNOD3_cut\tNOD3_time\tCOV_obj\tCOV_cut\tCOV_time\tEDGE_obj\tEDGE_cut\tEDGE_time\n';
file_l = '%s\t%10.6f\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\n';

f = fopen(fullfile(work_dir, 'DIMACS_separation_cut_selection.txt'), 'w');
fprintf(f, file_h);

for i=1:length(graphs)
    splitted = true;
    g = char(graphs(i));
    fprintf(strcat('File: ', g, '\n'));
    extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_gamma.mat')), fullfile(model_folder, strcat(g, '_nod_gamma.mat')), splitted); 
    extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_theta.mat')), fullfile(model_folder, strcat(g, '_nod_theta.mat')), splitted); 
    extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_alpha.mat')), fullfile(model_folder, strcat(g, '_nod_alpha.mat')), splitted);
    extract_file(zipFile, fullfile(model_folder, strcat(g, '_cov.mat')), fullfile(model_folder, strcat(g, '_cov.mat')), splitted); 
    extract_file(zipFile, fullfile(model_folder, strcat(g, '_edge.mat')), fullfile(model_folder, strcat(g, '_edge.mat')), splitted); 
    
    model_path = fullfile(model_th_folder, strcat(g, '_th+.mat'));
    load(model_path);
    % Load theta model
    At_th = At;
    b_th = b;
    blk = {'s', s};
    % Solve theta model (i.e. no cuts)
    [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], OPTIONS);
    theta_plus_time = info.totaltime;
    % Nod_gamma
    [nod1_obj, nod1_cuts, nod1_time] = kelley_cutting_plane(model_folder,sol_folder, g, '_nod_gamma', At_th, b_th, X_th, blk, obj_th, theta_plus_time, OPTIONS, splitted);
    % Nod_theta
    [nod2_obj, nod2_cuts, nod2_time] = kelley_cutting_plane(model_folder,sol_folder, g, '_nod_theta', At_th, b_th, X_th, blk, obj_th, theta_plus_time, OPTIONS, splitted);
    % Nod_alpha
    [nod3_obj, nod3_cuts, nod3_time] = kelley_cutting_plane(model_folder,sol_folder, g, '_nod_alpha', At_th, b_th, X_th, blk, obj_th, theta_plus_time, OPTIONS, splitted);
    % Cov
    [cov_obj, cov_cuts, cov_time] = kelley_cutting_plane(model_folder,sol_folder, g, '_cov', At_th, b_th, X_th, blk, obj_th, theta_plus_time, OPTIONS, splitted);
    % Edge
    [edge_obj, edge_cuts, edge_time] = kelley_cutting_plane(model_folder,sol_folder, g, '_edge', At_th, b_th, X_th, blk, obj_th, theta_plus_time, OPTIONS, splitted);

    if splitted
        for k=1:num_pieces
            delete(fullfile(model_folder, strcat(g, '_nod_gamma',sprintf('_%d.mat', k))), ...
            fullfile(model_folder, strcat(g, '_nod_theta',sprintf('_%d.mat', k))), ...
            fullfile(model_folder, strcat(g, '_nod_alpha',sprintf('_%d.mat', k))), ...
            fullfile(model_folder, strcat(g, '_cov',sprintf('_%d.mat', k))), ...
            fullfile(model_folder, strcat(g, '_edge',sprintf('_%d.mat', k))));
        end
    else
        delete(fullfile(model_folder, strcat(g, '_nod_gamma.mat')), ...
            fullfile(model_folder, strcat(g, '_nod_theta.mat')), ...
            fullfile(model_folder, strcat(g, '_nod_alpha.mat')), ...
            fullfile(model_folder, strcat(g, '_cov.mat')), ...
            fullfile(model_folder, strcat(g, '_edge.mat')));
    end

    fprintf(f, file_l, g, ...
        -obj_th(1), theta_plus_time, ...
        -nod1_obj, nod1_cuts, nod1_time, ...
        -nod2_obj, nod2_cuts, nod2_time, ...
        -nod3_obj, nod3_cuts, nod3_time, ...
        -cov_obj, cov_cuts, cov_time, ...
        -edge_obj, edge_cuts, edge_time);
end
% diary off;
zipFile.close();
fclose(f);
% entries = zipFile.getEntries; % to get all entries then iterate
