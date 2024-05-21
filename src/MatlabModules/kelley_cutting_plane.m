%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KELLEY CUTTING PLANE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each iteration:
%    1. Find violated cuts from the classes
%    2. When tailOff or no violated cuts are found, pass to the next class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nod_curr_obj, nod_tot_cuts, nod_time] = kelley_cutting_plane(model_folder,sol_folder, g, model_ext, At_th, b_th, X_th, blk, obj_th, theta_plus_time, Params)
 
splitted = Params.DO_SPLIT;
OPTIONS = Params.sdpnal;

if isfield(Params, 'CUTTING_PLANE_EPSILON')
    epsilon = Params.CUTTING_PLANE_EPSILON;
else
    epsilon = 1e-3;
end

if isfield(Params, 'CUTTING_PLANE_TAILOFF')
    delta = Params.CUTTING_PLANE_TAILOFF;
else
    delta = 1e-1;
end

if isfield(Params, 'CUTTING_PLANE_MAX_CUTS_PER_ITER')
    max_c = Params.CUTTING_PLANE_MAX_CUTS_PER_ITER;
else
    max_c = 1000;
end

num_classes = 4;
class_to_select = 1;
num_pieces = 5;
% Find cuts from M+
if splitted
    for i=1:num_pieces
        piece = fullfile(model_folder, strcat(g, model_ext, sprintf('_%d.mat', i)));
        load(piece);
    end
    Bt = [Bt_1 Bt_2 Bt_3 Bt_4];
else
    model_path = fullfile(model_folder, strcat(g, model_ext, '.mat'));
    fprintf('\nLoading: %s\n', model_path);
    % Load nodal model
    load(model_path);
end
start = tic;
Bt_nod = Bt;
u_nod = u;
cut_classes = cut_classes';

curr_X = X_th;
curr_obj = obj_th(1);
new_obj = Inf;
cutsBt = [ ];
cuts_u = [ ];
tot_cuts = 0;
not_done = 1;
n_cuts = 1;
init_obj = curr_obj;
added_cuts_idx = [];
while not_done
    to_add = [];
    X_vec = svec(blk, curr_X{1,1});
    % NOD(G) constraint violation from Theta X* opt sol
    violations = AXfun(blk, Bt_nod, X_vec) - u_nod;
    [ord_viol, I] = sort(violations, 'descend');
    % violations > 0 means that the constraint isn t satisfied by X*
    idx = find(ord_viol>epsilon);
    I = I(idx);
    violated_classes = cut_classes(I);
    I_i = I(violated_classes==class_to_select);
    to_add = [to_add; I_i];
    % ord_viol = ord_viol(idx);
    % take the most violated constr and add to theta
    n_cuts = min(max_c, size(to_add, 1));
    if n_cuts > 0
        tot_cuts = tot_cuts + n_cuts;
        cutsBt = [cutsBt Bt_nod(:, to_add(1:n_cuts))];
        cuts_u = [cuts_u; u_nod(to_add(1:n_cuts))];
        %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
        [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS);
        new_obj = new_obj(1);
        added_cuts_idx = [added_cuts_idx; to_add(1:n_cuts)];
        fprintf('\n');
    else
        fprintf('No violated cuts \n')
    end
    improve = abs(curr_obj - new_obj);
    if improve < delta || n_cuts == 0
        fprintf('Finished class %d. Tot cuts %d\n', class_to_select, tot_cuts);
        if class_to_select == num_classes
            not_done = 0;
        else
            class_to_select = class_to_select + 1;
        end
    end
    curr_obj = new_obj;
    fprintf('Obj changed of %13.5e\n', improve)
end
tend = toc(start);
tot_impro = abs(curr_obj - init_obj);
fprintf('Total Improvement from theta: %13.5e\n', tot_impro);
fprintf('Total Cuts added: %13d\n', tot_cuts);
fprintf('Total time: %13.2f\n', tend);

th_plus = init_obj;
nod_curr_obj = curr_obj;
if curr_obj == inf
    nod_curr_obj = th_plus;
end
nod_tot_cuts = tot_cuts;
nod_time = tend + theta_plus_time;

save(fullfile(sol_folder, strcat(g, model_ext, '_viol_test.mat')) , 'added_cuts_idx');
end


