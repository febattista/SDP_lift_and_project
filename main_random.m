%%% Options
clear
restoredefaultpath;
addpath(genpath(strcat(pwd, '/SDPNAL+v1.0')),path);
OPTIONS = SDPNALplus_parameters;

work_dir = fullfile('StableSets', 'Random-np');
model_folder = fullfile('StableSets', 'Random-np', 'models');
model_th_folder = fullfile('StableSets', 'Random-np', 'models_th');
sol_folder = fullfile('StableSets', 'Random-np', 'violated_cuts');
zipFilename = 'models_random.zip';
zipJavaFile  = java.io.File(fullfile(model_folder, zipFilename));
zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);

if exist(sol_folder, 'dir') == 0
    mkdir(sol_folder);
end

num_pieces = 4;

ns = [150 175 200 225 250 275 300];
ds = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
is = [1 2 3 4 5];
filename = 'G_%d_%1.1f_%d';

% log_file = 'matlab_cut_selection_2.log';
% if exist(log_file, 'file')
%     delete(log_file);
% end
% diary(log_file);

file_h = 'n\td\ti\tTH+\tTH+_time\tNOD1_obj\tNOD1_cut\tNOD1_time\tNOD2_obj\tNOD2_cut\tNOD2_time\tNOD3_obj\tNOD3_cut\tNOD3_time\tCOV_obj\tCOV_cut\tCOV_time\tEDGE_obj\tEDGE_cut\tEDGE_time\n';
file_l = '%d\t%1.1f\t%d\t%10.6f\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\n';

f = fopen(fullfile(work_dir, 'random_separation_cut_selection_3.txt'), 'w');
fprintf(f, file_h);

for n=1:length(ns)
    for d=1:length(ds)
        for i=1:length(is)
            splitted = false;
            g = sprintf(filename, ns(n), ds(d), is(i));
            fprintf(strcat('File: ', g, '\n'));
           extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_gamma.mat')), fullfile(model_folder, strcat(g, '_nod_gamma.mat')), splitted); 
           extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_theta.mat')), fullfile(model_folder, strcat(g, '_nod_theta.mat')), splitted); 
           extract_file(zipFile, fullfile(model_folder, strcat(g, '_nod_alpha.mat')), fullfile(model_folder, strcat(g, '_nod_alpha.mat')), splitted);
           extract_file(zipFile, fullfile(model_folder, strcat(g, '_edge.mat')), fullfile(model_folder, strcat(g, '_edge.mat')), splitted); 
           if ns(n) == 275
               extract_file(zipFile, fullfile(model_folder, strcat(g, '_clique.mat')), fullfile(model_folder, strcat(g, '_cov.mat')), splitted); 
           else    
               extract_file(zipFile, fullfile(model_folder, strcat(g, '_cov.mat')), fullfile(model_folder, strcat(g, '_cov.mat')), splitted); 
           end
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
            
            delete(fullfile(model_folder, strcat(g, '_nod_gamma.mat')), ...
                   fullfile(model_folder, strcat(g, '_nod_theta.mat')), ...
                   fullfile(model_folder, strcat(g, '_nod_alpha.mat')), ...
                   fullfile(model_folder, strcat(g, '_cov.mat')), ...
                   fullfile(model_folder, strcat(g, '_edge.mat')));

            if splitted
                for k=1:num_pieces
                    delete(fullfile(model_folder, strcat(g, '_nod_gamma',sprintf('_%d.mat', k))), ...
                   fullfile(model_folder, strcat(g, '_nod_theta',sprintf('_%d.mat', k))), ...
                   fullfile(model_folder, strcat(g, '_nod_alpha',sprintf('_%d.mat', k))), ...
                   fullfile(model_folder, strcat(g, '_cov',sprintf('_%d.mat', k))), ...
                   fullfile(model_folder, strcat(g, '_edge',sprintf('_%d.mat', k))));
                end
            end
        
           fprintf(f, file_l, ns(n), ds(d), is(i), ...
               -obj_th(1), theta_plus_time, ...
               -nod1_obj, nod1_cuts, nod1_time, ...
               -nod2_obj, nod2_cuts, nod2_time, ...
               -nod3_obj, nod3_cuts, nod3_time, ...
               -cov_obj, cov_cuts, cov_time, ...
               -edge_obj, edge_cuts, edge_time);
        end
    end
end
% diary off;
zipFile.close();
fclose(f);
% entries = zipFile.getEntries; % to get all entries then iterate
