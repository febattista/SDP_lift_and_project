%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KELLEY CUTTING PLANE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each iteration:
%    1. Find violated cuts from the classes
%    2. When tailOff or no violated cuts are found, pass to the next class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_curr_obj, R_tot_cuts, R_time, R_iter] = kelley_cutting_plane(model_folder,sol_folder, g, model_ext, At_th, b_th, X_th, blk, obj_th, theta_plus_time, Params)

fprintf('\n');

% Initialize parameters
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

if isfield(Params, 'CUTTING_PLANE_TIMELIMIT')
    timelimit = Params.CUTTING_PLANE_TIMELIMIT;
else
    timelimit = 7200;
end

num_classes = 4;        % Total number of cut classes
class_to_select = 1;    % Class to start to select the cuts from
num_pieces = 5;         % Num of pieces .mat files are splitted

% Load the SDP formulation R that is a restriction of Theta_plus and
% cuts will be selected from
fprintf('Model: %s%s\n', g, model_ext);
if splitted
    for i=1:num_pieces
        piece = fullfile(model_folder, strcat(g, model_ext, sprintf('_%d.mat', i)));
        load(piece);
    end
    Bt = [Bt_1 Bt_2 Bt_3 Bt_4];
else
    model_path = fullfile(model_folder, strcat(g, model_ext, '.mat'));
    load(model_path);
end

iter = 0;                   % Iterations count: an iteration is counted when  
                            % vioated cuts are identified and the SDP is resolved

start = tic;                % Start timer
Bt_R = Bt;                  % LHS matrix of R
u_R = u;                    % RHS vector of R
cut_classes = cut_classes'; % Class for each ineq in R
curr_X = X_th;              % Current opt sol matrix
curr_obj = obj_th(1);       % Current opt obj val
init_obj = curr_obj;        % Obj val before cuts have been added
new_obj = Inf;              % Obj val after cuts have been added
cuts_Bt = [ ];              % LHS of cuts added so far
cuts_u = [ ];               % RHS of cuts added so far
added_cuts_idx = [];        % Idx of cuts added so far
tot_cuts = 0;               % Num cuts added so far
not_done = 1;               % Are we done?
n_cuts = 1;                 % Num of cuts we are going to add in this iter

while not_done
    fprintf('Iter %d -----------------------------------\n', iter);
    % compute violations and look for cuts to be added
    % violations[i] > 0 means that i-th constraint is not satisfied by X*
    % take the most violated constr and add to theta
    to_add = [];
    X_vec = svec(blk, curr_X{1,1});
    violations = AXfun(blk, Bt_R, X_vec) - u_R;
    [ord_viol, I] = sort(violations, 'descend');
    idx = find(ord_viol>epsilon);
    I = I(idx);
    violated_classes = cut_classes(I);
    I_i = I(violated_classes==class_to_select);
    to_add = [to_add; I_i];
    % ord_viol = ord_viol(idx);
    n_cuts = min(max_c, size(to_add, 1));

    if n_cuts > 0
        fprintf('Found %d violated cuts from class %d\n', n_cuts, class_to_select);
        iter = iter + 1;
        % include the new cuts and re-optimize the SDP
        tot_cuts = tot_cuts + n_cuts;
        cuts_Bt = [cuts_Bt Bt_R(:, to_add(1:n_cuts))];
        cuts_u = [cuts_u; u_R(to_add(1:n_cuts))];
        % Warm starting SDPNALplus does not improve sol time
        %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
        [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cuts_Bt},[],cuts_u, OPTIONS);
        new_obj = new_obj(1);
        added_cuts_idx = [added_cuts_idx; to_add(1:n_cuts)];
        fprintf('\n');
    else
        fprintf('No violated cuts \n')
        new_obj = curr_obj;
    end

    % check tailOff and termination criteria
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
    tend = toc(start);
    fprintf('Time at this iter: %13.2f\n', tend);
    if tend > timelimit
        not_done = 0;
        fprintf('Time limit reached\n')
    end
end
tot_impro = abs(curr_obj - init_obj);
fprintf('========================================================\n')
fprintf('Total Improvement from theta: %13.5e\n', tot_impro);
fprintf('Total Cuts added: %13d\n', tot_cuts);
fprintf('Total time: %13.2f\n', tend);
fprintf('========================================================\n')

th_plus = init_obj;
R_curr_obj = curr_obj;
if curr_obj == inf
    R_curr_obj = th_plus;
end
R_tot_cuts = tot_cuts;
R_time = tend + theta_plus_time;
R_iter = iter;

save(fullfile(sol_folder, strcat(g, model_ext, '_viol_test.mat')) , 'added_cuts_idx');
end


