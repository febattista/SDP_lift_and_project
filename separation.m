%%% Options
restoredefaultpath;
addpath(genpath(strcat(pwd, '/SDPNAL+v1.0')),path);
OPTIONS = SDPNALplus_parameters;
work_dir = fullfile('StableSets','Random-np');
file_out_header = 'n\td\ti\tTH+\tN-th_last_obj\ttot_cuts\ttime\n';
file_out_line = '%s\t%s\t%s\t%10.6f\t%10.6f\t%d\t%10.2f\n';

   models_name = {...
'G_150_0.1_1', 'G_150_0.1_2', 'G_150_0.1_3', 'G_150_0.1_4', 'G_150_0.1_5', ...
'G_150_0.2_1', 'G_150_0.2_2', 'G_150_0.2_3', 'G_150_0.2_4', 'G_150_0.2_5', ...
'G_150_0.3_1', 'G_150_0.3_2', 'G_150_0.3_3', 'G_150_0.3_4', 'G_150_0.3_5', ...
'G_150_0.4_1', 'G_150_0.4_2', 'G_150_0.4_3', 'G_150_0.4_4', 'G_150_0.4_5', ...
'G_150_0.5_1', 'G_150_0.5_2', 'G_150_0.5_3', 'G_150_0.5_4', 'G_150_0.5_5', ...
'G_150_0.6_1', 'G_150_0.6_2', 'G_150_0.6_3', 'G_150_0.6_4', 'G_150_0.6_5', ...
'G_150_0.7_1', 'G_150_0.7_2', 'G_150_0.7_3', 'G_150_0.7_4', 'G_150_0.7_5', ...
'G_150_0.8_1', 'G_150_0.8_2', 'G_150_0.8_3', 'G_150_0.8_4', 'G_150_0.8_5', ...
'G_150_0.9_1', 'G_150_0.9_2', 'G_150_0.9_3', 'G_150_0.9_4', 'G_150_0.9_5', ...
'G_175_0.1_1', 'G_175_0.1_2', 'G_175_0.1_3', 'G_175_0.1_4', 'G_175_0.1_5', ...
'G_175_0.2_1', 'G_175_0.2_2', 'G_175_0.2_3', 'G_175_0.2_4', 'G_175_0.2_5', ...
'G_175_0.3_1', 'G_175_0.3_2', 'G_175_0.3_3', 'G_175_0.3_4', 'G_175_0.3_5', ...
'G_175_0.4_1', 'G_175_0.4_2', 'G_175_0.4_3', 'G_175_0.4_4', 'G_175_0.4_5', ...
'G_175_0.5_1', 'G_175_0.5_2', 'G_175_0.5_3', 'G_175_0.5_4', 'G_175_0.5_5', ...
'G_175_0.6_1', 'G_175_0.6_2', 'G_175_0.6_3', 'G_175_0.6_4', 'G_175_0.6_5', ...
'G_175_0.7_1', 'G_175_0.7_2', 'G_175_0.7_3', 'G_175_0.7_4', 'G_175_0.7_5', ...
'G_175_0.8_1', 'G_175_0.8_2', 'G_175_0.8_3', 'G_175_0.8_4', 'G_175_0.8_5', ...
'G_175_0.9_1', 'G_175_0.9_2', 'G_175_0.9_3', 'G_175_0.9_4', 'G_175_0.9_5', ...
'G_200_0.1_1', 'G_200_0.1_2', 'G_200_0.1_3', 'G_200_0.1_4', 'G_200_0.1_5', ...
'G_200_0.2_1', 'G_200_0.2_2', 'G_200_0.2_3', 'G_200_0.2_4', 'G_200_0.2_5', ...
'G_200_0.3_1', 'G_200_0.3_2', 'G_200_0.3_3', 'G_200_0.3_4', 'G_200_0.3_5', ...
'G_200_0.4_1', 'G_200_0.4_2', 'G_200_0.4_3', 'G_200_0.4_4', 'G_200_0.4_5', ...
'G_200_0.5_1', 'G_200_0.5_2', 'G_200_0.5_3', 'G_200_0.5_4', 'G_200_0.5_5', ...
'G_200_0.6_1', 'G_200_0.6_2', 'G_200_0.6_3', 'G_200_0.6_4', 'G_200_0.6_5', ...
'G_200_0.7_1', 'G_200_0.7_2', 'G_200_0.7_3', 'G_200_0.7_4', 'G_200_0.7_5', ...
'G_200_0.8_1', 'G_200_0.8_2', 'G_200_0.8_3', 'G_200_0.8_4', 'G_200_0.8_5', ...
'G_200_0.9_1', 'G_200_0.9_2', 'G_200_0.9_3', 'G_200_0.9_4', 'G_200_0.9_5', ...
'G_225_0.1_1', 'G_225_0.1_2', 'G_225_0.1_3', 'G_225_0.1_4', 'G_225_0.1_5', ...
'G_225_0.2_1', 'G_225_0.2_2', 'G_225_0.2_3', 'G_225_0.2_4', 'G_225_0.2_5', ...
'G_225_0.3_1', 'G_225_0.3_2', 'G_225_0.3_3', 'G_225_0.3_4', 'G_225_0.3_5', ...
'G_225_0.4_1', 'G_225_0.4_2', 'G_225_0.4_3', 'G_225_0.4_4', 'G_225_0.4_5', ...
'G_225_0.5_1', 'G_225_0.5_2', 'G_225_0.5_3', 'G_225_0.5_4', 'G_225_0.5_5', ...
'G_225_0.6_1', 'G_225_0.6_2', 'G_225_0.6_3', 'G_225_0.6_4', 'G_225_0.6_5', ...
'G_225_0.7_1', 'G_225_0.7_2', 'G_225_0.7_3', 'G_225_0.7_4', 'G_225_0.7_5', ...
'G_225_0.8_1', 'G_225_0.8_2', 'G_225_0.8_3', 'G_225_0.8_4', 'G_225_0.8_5', ...
'G_225_0.9_1', 'G_225_0.9_2', 'G_225_0.9_3', 'G_225_0.9_4', 'G_225_0.9_5', ...
'G_250_0.1_1', 'G_250_0.1_2', 'G_250_0.1_3', 'G_250_0.1_4', 'G_250_0.1_5', ...
'G_250_0.2_1', 'G_250_0.2_2', 'G_250_0.2_3', 'G_250_0.2_4', 'G_250_0.2_5', ...
'G_250_0.3_1', 'G_250_0.3_2', 'G_250_0.3_3', 'G_250_0.3_4', 'G_250_0.3_5', ...
'G_250_0.4_1', 'G_250_0.4_2', 'G_250_0.4_3', 'G_250_0.4_4', 'G_250_0.4_5', ...
'G_250_0.5_1', 'G_250_0.5_2', 'G_250_0.5_3', 'G_250_0.5_4', 'G_250_0.5_5', ...
'G_250_0.6_1', 'G_250_0.6_2', 'G_250_0.6_3', 'G_250_0.6_4', 'G_250_0.6_5', ...
'G_250_0.7_1', 'G_250_0.7_2', 'G_250_0.7_3', 'G_250_0.7_4', 'G_250_0.7_5', ...
'G_250_0.8_1', 'G_250_0.8_2', 'G_250_0.8_3', 'G_250_0.8_4', 'G_250_0.8_5', ...
'G_250_0.9_1', 'G_250_0.9_2', 'G_250_0.9_3', 'G_250_0.9_4', 'G_250_0.9_5', ...
'G_275_0.1_1', 'G_275_0.1_2', 'G_275_0.1_3', 'G_275_0.1_4', 'G_275_0.1_5', ...
'G_275_0.2_1', 'G_275_0.2_2', 'G_275_0.2_3', 'G_275_0.2_4', 'G_275_0.2_5', ...
'G_275_0.3_1', 'G_275_0.3_2', 'G_275_0.3_3', 'G_275_0.3_4', 'G_275_0.3_5', ...
'G_275_0.4_1', 'G_275_0.4_2', 'G_275_0.4_3', 'G_275_0.4_4', 'G_275_0.4_5', ...
'G_275_0.5_1', 'G_275_0.5_2', 'G_275_0.5_3', 'G_275_0.5_4', 'G_275_0.5_5', ...
'G_275_0.6_1', 'G_275_0.6_2', 'G_275_0.6_3', 'G_275_0.6_4', 'G_275_0.6_5', ...
'G_275_0.7_1', 'G_275_0.7_2', 'G_275_0.7_3', 'G_275_0.7_4', 'G_275_0.7_5', ...
'G_275_0.8_1', 'G_275_0.8_2', 'G_275_0.8_3', 'G_275_0.8_4', 'G_275_0.8_5', ...
'G_275_0.9_1', 'G_275_0.9_2', 'G_275_0.9_3', 'G_275_0.9_4', 'G_275_0.9_5', ...
'G_300_0.1_1', 'G_300_0.1_2', 'G_300_0.1_3', 'G_300_0.1_4', 'G_300_0.1_5', ...
'G_300_0.2_1', 'G_300_0.2_2', 'G_300_0.2_3', 'G_300_0.2_4', 'G_300_0.2_5', ...
'G_300_0.3_1', 'G_300_0.3_2', 'G_300_0.3_3', 'G_300_0.3_4', 'G_300_0.3_5', ...
'G_300_0.4_1', 'G_300_0.4_2', 'G_300_0.4_3', 'G_300_0.4_4', 'G_300_0.4_5', ...
'G_300_0.5_1', 'G_300_0.5_2', 'G_300_0.5_3', 'G_300_0.5_4', 'G_300_0.5_5', ...
'G_300_0.6_1', 'G_300_0.6_2', 'G_300_0.6_3', 'G_300_0.6_4', 'G_300_0.6_5', ...
'G_300_0.7_1', 'G_300_0.7_2', 'G_300_0.7_3', 'G_300_0.7_4', 'G_300_0.7_5', ...
'G_300_0.8_1', 'G_300_0.8_2', 'G_300_0.8_3', 'G_300_0.8_4', 'G_300_0.8_5', ...
'G_300_0.9_1', 'G_300_0.9_2', 'G_300_0.9_3', 'G_300_0.9_4', 'G_300_0.9_5', ...
};

%models_name = {...
%'1zc512', '1zc1024', ...
%                'brock200_1', 'brock200_2', 'brock200_3', 'brock200_4', ...
%                'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', ...
%                'brock800_1', 'brock800_2', 'brock800_3', 'brock800_4', ...
%                'C125-9', 'C250-9', 'C500-9',...
%                'DSJC125.1', 'DSJC125.5', 
%                'p_hat300-1', 'p_hat300-2', 'p_hat300-3', ...
%                'p_hat500-1', 'p_hat500-2', 'p_hat500-3', ...
%                'p_hat700-1', 'p_hat700-2', 'p_hat700-3', ...
%                'sanr200_0.7', 'sanr200_0.9',  ...
%                'sanr400_0.5', 'sanr400_0.7'};
%
%
%work_dir = fullfile('StableSets','DIMACS');
%file_out_header = 'n\tTH+\tTH+_time\tNOD_last_obj\tNOD_tot_cuts\tNOD_time\tLS_last_obj\tLS_tot_cuts\tLS_time\tGR_last_obj\tGR_tot_cuts\tGR_time\tDR_last_obj\tDR_tot_cuts\tDR_time\n';
%file_out_line = '%s\t%10.6f\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f\t%10.6f\t%d\t%10.2f \n';
file_out = fullfile(work_dir, 'output_rand_N_th_separation.txt');
f = fopen(file_out, 'w');
fprintf(f, file_out_header);

models = fullfile(work_dir, 'models_th', models_name);

for i=1:length(models)
	filename = char(models_name(i));
    n = filename(3:5);
    d = filename(9);
    ii = filename(end);
    %%%%%% NODAL
    model_path = strcat(models(i), '_th+.mat');
    load(model_path{1});
    % Load theta model
    At_th = At;
    b_th = b;
    blk = {'s', s};
    
    % Solve theta model (i.e. no cuts)
    [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], OPTIONS);
    theta_plus_time = info.totaltime;
    % Find cuts from M+(NOD(G))
    model_path = strcat(models(i), '_nod_th.mat');
    % Load nodal model
    load(model_path{1});
    start = tic;
    Bt_nod = Bt;
    u_nod = u;

    curr_X = X_th;
    curr_obj = obj_th(1);
    new_obj = Inf;
    n_cuts = 1;
    cutsBt = [ ];
    cuts_u = [ ];
    tot_cuts = 0;
    not_done = 1;
    init_obj = curr_obj;
    while not_done
        X_vec = svec(blk, curr_X{1,1});
        % NOD(G) constraint violation from Theta X* opt sol
        violations = AXfun(blk, Bt_nod, X_vec) - u_nod;
        [ord_viol, I] = sort(violations, 'descend');
        % violations > 0 means that the constraint isn t satisfied by X*
        idx = find(ord_viol>1e-3);
        I = I(idx);
        ord_viol = ord_viol(idx);
        % take the most violated constr and add to theta
        n_cuts = min(1000, size(ord_viol, 1));
        if n_cuts > 0
            tot_cuts = tot_cuts + n_cuts;
            cutsBt = [cutsBt Bt_nod(:, I(1:n_cuts))];
            cuts_u = [cuts_u; u_nod(I(1:n_cuts))];
            %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
            [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS);
            new_obj = new_obj(1);
        else
            fprintf('No violated cuts \n')
        end
        improve = abs(curr_obj - new_obj);
        if improve < 1e-1 || n_cuts == 0
            not_done = 0;
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

%    if contains(filename, {'1zc512', ...
%                'brock200_1', 'brock200_2', 'brock200_3', 'brock200_4', ...
%                'brock400_1', 'brock400_2', 'brock400_3', 'brock400_4', ...
%                'C125-9', 'C250-9', 'C500-9',...
%                'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC500-5', ...  
%                'keller4',  ...
%                'MANN_a9', 'MANN_a27', ...
%                'p_hat300-1', 'p_hat300-2', 'p_hat300-3', ...
%                 'p_hat500-2', 'p_hat500-3', ...
%                'sanr200_0.7', 'sanr200_0.9',  ...
%                'sanr400_0.5', 'sanr400_0.7'})
%    
%    %%%%%% LOVASZ-SCHRIJVER
%    model_path = strcat(models(i), '_th+.mat');
%    load(model_path{1});
%    % Load theta model
%    At_th = At;
%    b_th = b;
%    blk = {'s', s};
%    
%    % Solve theta model (i.e. no cuts)
%    [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], OPTIONS);
%
%    % Find cuts from M+(FRAC(G))
%    model_path = strcat(models(i), '_ls.mat');
%    % Load nodal model
%    load(model_path{1});
%    start = tic;
%    Bt_nod = Bt;
%    u_nod = u;
%
%    curr_X = X_th;
%    curr_obj = obj_th(1);
%    new_obj = Inf;
%    n_cuts = 1;
%    cutsBt = [ ];
%    cuts_u = [ ];
%    tot_cuts = 0;
%    not_done = 1;
%    init_obj = curr_obj;
%    while not_done
%        X_vec = svec(blk, curr_X{1,1});
%        % NOD(G) constraint violation from Theta X* opt sol
%        violations = AXfun(blk, Bt_nod, X_vec) - u_nod;
%        [ord_viol, I] = sort(violations, 'descend');
%        % violations > 0 means that the constraint isn t satisfied by X*
%        idx = find(ord_viol>1e-3);
%        I = I(idx);
%        ord_viol = ord_viol(idx);
%        % take the most violated constr and add to theta
%        n_cuts = min(1000, size(ord_viol, 1));
%        if n_cuts > 0
%            tot_cuts = tot_cuts + n_cuts;
%            cutsBt = [cutsBt Bt_nod(:, I(1:n_cuts))];
%            cuts_u = [cuts_u; u_nod(I(1:n_cuts))];
%            %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
%            [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS);
%            new_obj = new_obj(1);
%        else
%            fprintf('No violated cuts \n')
%        end
%        improve = abs(curr_obj - new_obj);
%        if improve < 1e-1 || n_cuts == 0
%            not_done = 0;
%        end
%        curr_obj = new_obj;
%        fprintf('Obj changed of %13.5e\n', improve)
%    end
%    tend = toc(start);
%    tot_impro = abs(curr_obj - init_obj);
%    fprintf('Total Improvement from theta: %13.5e\n', tot_impro);
%    fprintf('Total Cuts added: %13d\n', tot_cuts);
%    fprintf('Total time: %13.2f\n', tend);
%    
%    ls_curr_obj = curr_obj;
%    if curr_obj == inf
%        ls_curr_obj = th_plus;
%    end
%    ls_tot_cuts = tot_cuts;
%    ls_time = tend + theta_plus_time;
%    
%    %%%%%% GRUBER-RENDL
%    model_path = strcat(models(i), '_th+.mat');
%    load(model_path{1});
%    % Load theta model
%    At_th = At;
%    b_th = b;
%    blk = {'s', s};
%    
%    % Solve theta model (i.e. no cuts)
%    [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], OPTIONS);
%
%    % Find cuts from GR
%    model_path = strcat(models(i), '_gr.mat');
%    % Load nodal model
%    load(model_path{1});
%    start = tic;
%    Bt_nod = Bt;
%    u_nod = u;
%
%    curr_X = X_th;
%    curr_obj = obj_th(1);
%    new_obj = Inf;
%    n_cuts = 1;
%    cutsBt = [ ];
%    cuts_u = [ ];
%    tot_cuts = 0;
%    not_done = 1;
%    init_obj = curr_obj;
%    while not_done
%        X_vec = svec(blk, curr_X{1,1});
%        % NOD(G) constraint violation from Theta X* opt sol
%        violations = AXfun(blk, Bt_nod, X_vec) - u_nod;
%        [ord_viol, I] = sort(violations, 'descend');
%        % violations > 0 means that the constraint isn t satisfied by X*
%        idx = find(ord_viol>1e-3);
%        I = I(idx);
%        ord_viol = ord_viol(idx);
%        % take the most violated constr and add to theta
%        n_cuts = min(1000, size(ord_viol, 1));
%        if n_cuts > 0
%            tot_cuts = tot_cuts + n_cuts;
%            cutsBt = [cutsBt Bt_nod(:, I(1:n_cuts))];
%            cuts_u = [cuts_u; u_nod(I(1:n_cuts))];
%            %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
%            [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS);
%            new_obj = new_obj(1);
%        else
%            fprintf('No violated cuts \n')
%        end
%        improve = abs(curr_obj - new_obj);
%        if improve < 1e-1 || n_cuts == 0
%            not_done = 0;
%        end
%        curr_obj = new_obj;
%        fprintf('Obj changed of %13.5e\n', improve)
%    end
%    tend = toc(start);
%    tot_impro = abs(curr_obj - init_obj);
%    fprintf('Total Improvement from theta: %13.5e\n', tot_impro);
%    fprintf('Total Cuts added: %13d\n', tot_cuts);
%    fprintf('Total time: %13.2f\n', tend);
%    
%    gr_curr_obj = curr_obj;
%    if curr_obj == inf
%        gr_curr_obj = th_plus;
%    end
%    gr_tot_cuts = tot_cuts;
%    gr_time = tend + theta_plus_time;
%    
%    %%%%%% Dukanovich-Rendl
%    model_path = strcat(models(i), '_th1+.mat');
%    load(model_path{1});
%    % Load theta model
%    At_th = At;
%    b_th = b;
%    blk = {'s', s};
%    
%    % Solve theta model (i.e. no cuts)
%    [obj_th,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th,info,runhist] = sdpnalplus(blk, {At_th}, {C}, b_th, L,[],[],[],[], OPTIONS);
%    theta_plus_time1 = info.totaltime;
%    % Find cuts from M+(FRAC(G))
%    model_path = strcat(models(i), '_dr.mat');
%    % Load nodal model
%    load(model_path{1});
%    start = tic;
%    Bt_nod = Bt;
%    u_nod = u;
%
%    curr_X = X_th;
%    curr_obj = obj_th(1);
%    new_obj = Inf;
%    n_cuts = 1;
%    cutsBt = [ ];
%    cuts_u = [ ];
%    tot_cuts = 0;
%    not_done = 1;
%    init_obj = curr_obj;
%    while not_done
%        X_vec = svec(blk, curr_X{1,1});
%        % NOD(G) constraint violation from Theta X* opt sol
%        violations = AXfun(blk, Bt_nod, X_vec) - u_nod;
%        [ord_viol, I] = sort(violations, 'descend');
%        % violations > 0 means that the constraint isn t satisfied by X*
%        idx = find(ord_viol>1e-3);
%        I = I(idx);
%        ord_viol = ord_viol(idx);
%        % take the most violated constr and add to theta
%        n_cuts = min(1000, size(ord_viol, 1));
%        if n_cuts > 0
%            tot_cuts = tot_cuts + n_cuts;
%            cutsBt = [cutsBt Bt_nod(:, I(1:n_cuts))];
%            cuts_u = [cuts_u; u_nod(I(1:n_cuts))];
%            %[new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS,X_th,s_th,y_th,S_th,Z_th,ybar_th,v_th);
%            [new_obj,curr_X,s,y,S,Z,ybar,v,info,runhist] = sdpnalplus(blk,{At_th},{C},b_th, L,[],{cutsBt},[],cuts_u, OPTIONS);
%            new_obj = new_obj(1);
%        else
%            fprintf('No violated cuts \n')
%        end
%        improve = abs(curr_obj - new_obj);
%        if improve < 1e-1 || n_cuts == 0
%            not_done = 0;
%        end
%        curr_obj = new_obj;
%        fprintf('Obj changed of %13.5e\n', improve)
%    end
%    tend = toc(start);
%    tot_impro = abs(curr_obj - init_obj);
%    fprintf('Total Improvement from theta: %13.5e\n', tot_impro);
%    fprintf('Total Cuts added: %13d\n', tot_cuts);
%    fprintf('Total time: %13.2f\n', tend);
%    
%    dr_curr_obj = curr_obj;
%    if curr_obj == inf
%        dr_curr_obj = th_plus;
%    end
%    
%    dr_tot_cuts = tot_cuts;
%    dr_time = tend + theta_plus_time1;
%
%    else
%        ls_curr_obj = 1; ls_tot_cuts = -1; ls_time = -1;
%        gr_curr_obj = 1; gr_tot_cuts = -1; gr_time = -1; 
%        dr_curr_obj = 1; dr_tot_cuts = -1; dr_time = -1;
%    end
%    
%    
%    fprintf(f, file_out_line, filename, ...
%        -th_plus, theta_plus_time,...
%        -nod_curr_obj, nod_tot_cuts, nod_time, ...
%        -ls_curr_obj, ls_tot_cuts, ls_time,...
%        -gr_curr_obj, gr_tot_cuts, gr_time, ...
%        -dr_curr_obj, dr_tot_cuts, dr_time);
    fprintf(f, file_out_line, n, d, ii, ...
    -th_plus,...
    -nod_curr_obj, nod_tot_cuts, nod_time);
end
   






