% MAIN_FILE_IRMATCHING - Implements IR-matching a la CEE specifically for our two-sector model with spillover and news.
% Stuff specific to this model are indicated throughout so that
% the structure should be general and you can use this file with modified
% functions or inputs in it to estimate different models.

% STRUCTURE
% -1.) Some notational comments
%  0.) Some playing around with the paths.
% STAGE 1. Set up a matlab function that will solve the model w/o using the
% symbolic toolbox (faster for fmincon)
% STAGE 2. Generate/import empirical IRFs from VAR
% STAGE 3. Generate bootstrap IRFs and take their variance to get the
% weighting matrix.
% STAGE 4. Compute loss function and do fmincon.

% Code by Marco Brianti and Laura Gati, Boston College, 2018
%**************************************************************

%% % --->>>>> this is specific to BriantiGati2017
% This code specifically implements a matching of responses of an IT shock
% from our justIT VAR to the responses of our two-sector model with
% spillovers to a noise shock. (It may be worth considering other shocks.)
% This sign indicates the parts that are specific to this particular model.

%%

clear
close all
disp('Running IR matching. Don''t forget to edit specific parts.')
tic
%% Add all the relevant paths
current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)

addpath(base_path)

if exist([base_path '\Code'], 'dir')
    addpath([base_path '\Code\just_IT_vecm']) %for Microsoft
else
    addpath([base_path '/Code/just_IT_vecm']) %for Mac
end

%% Choose some stuff specific to this application
nshocks = 1; % to how many shocks do we wanna match IRFs (not sure if we can actually set this... maybe it's actually = nvar_VAR?)

%% taken from Ryan directly (PS3)
%**********************************************************
% STAGE 1: Prepare a function to solve the model quickly
%          i.e. without using the symbolic toolbox for each
%          iteration
%**********************************************************

%Load Parameters
[param0,set] = paramsset;

%Compute the first-order derivative numerically
%%%%%
% How to use this model numerical approach properly:
% Write model.m (or edit it, to be correct, here: model_IRmatching_spillover_news.m)
% Then run model_func.m to generate model_prog.m (I called it model_prog_IRmatching_spillover_news.m but reverted to model_prog.m to avoid naming issues.)
% Then run model_prog.m; this will do the analytical evaluations.
%%%%%

% mod = model_IRmatching_spillover_news(param0,set);
mod = model_IRmatching_spillover_news2(param0,set); % a version that incorporates GDP and TFP

%Get dimensions of the model
nx = length(mod.X);
ny = length(mod.Y);
neq = nx+ny;

%Generate a matlab function that can compute the derivative matrices
%quickly
trans = zeros(1,neq); % this part specifies whether params should be transformed to all positive or negative
lb = -inf*trans;
ub = +inf*trans;
model_func(mod, trans, lb, ub);
rehash;
load sym_mod *idx


%Convert param and set to a vector of values: these are our starting values
%for the optimization.
param0 = struct2array(param0);
set = struct2array(set);


%Test intial values
[f, fx, fy, fxp, fyp, G, R, set]=model_prog(param0,set); % --->>>>> this is specific to BriantiGati2017
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp); 

%load checking_f % fx = fxn from the earlier solution of the model, cool

% mom_tab shows a little table of stddevs, autocorrs and corrs of the
% specified variables after a specific shock (G = eta*shock vector)
% mom_tab(gx,hx,G*G', [gamtfp_idx,gamki_idx], {'TFP','KI'});


%%
%**********************************************************
% STAGE 2. Input IRFs from VAR
%**********************************************************
% load('Workspace_Just_IT_SVAR_1LAG.mat') % --->>>>> this is specific to BriantiGati2017
load('main_just_IT_noHOURS_1LAG') % --->>>>> this is specific to BriantiGati2017

% Construct VAR IRFs to shocks of your choice
% Size of imported VAR IRFs is (nvar, T, nshocks)
nvar_VAR_orig = size(IRFs,1);
% T_VAR        = size(IRFs,2);
T_VAR        = 40;
nshocks_VAR  = size(IRFs,3);
% choose which variables responses to match from the VAR
variables_from_VAR = [1 2 4 5];
IRF_IT    = IRFs(variables_from_VAR,:,pos_IT); % response of all vars to IT shock % --->>>>> this is specific to BriantiGati2017
nvar_VAR  = size(IRF_IT,1);

% Gather the VAR IRFs to all the relevant shocks
IRFs_VAR = IRF_IT(:,1:T_VAR)'; % --->>>>> this is specific to BriantiGati2017
psi_hat  = IRFs_VAR(:);


%%
%**********************************************************
% STAGE 3. Bootstrap IRFs and weighting matrix
%**********************************************************
% Compute weighting matrix W = V^(-1), where V=var[bootstrap_IRFs] from the
% VAR.

do_we_need_V = 'no';
switch do_we_need_V
    case 'yes'
        % Get "bootstrapped IRFs" nsimul times
        % Here procedure differs from standard bootstrap because we want
        % bootstrapped IRFs, so to keep that clear I denote everything here by _s
        IRFs_s = zeros(nvar_VAR_orig, T_VAR, nshocks_VAR, nsimul);
        for i_simul=1:nsimul
            [A_s, B_s,~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
            % Get bootstrapped confidence intervals nsimul times
            disp(['Iteration ' num2str(i_simul) ' out of ' num2str(nsimul)])
            [impact_s, ~, ~]  = just_IT_ID(which_variable,which_shock,A_boot);  % --->>>>> this is specific to BriantiGati2017
            %Creating a fake matrix for the IRF due to partial ID
            fake_impact_s = zeros(nvar,nvar);
            fake_impact_s(:,which_shock) = impact_s;
            [IRFs_s(:,:,:,i_simul), ~, ~, ~, ~] = genIRFs(fake_impact_s,0,...
                B_s,0,T_VAR,sig1, sig2);
        end
        IRFs_s_IT = squeeze(IRFs_s(variables_from_VAR,:,pos_IT,:)); % --->>>>> this is specific to BriantiGati2017
        psi_boot = reshape(IRFs_s_IT,nvar_VAR*T_VAR,nsimul);
        V = diag(var(psi_boot,0,2));
        
        % alternative way to create V  -- both yield the same result, but inv(V) is
        % singular, so something is not yet quite right.
        V1 = var(IRFs_s_IT,0,3); % take variance over nsimul
        V2 = reshape(V1,size(V1,1)*size(V1,2),1);
        V_alt = diag(V2);
        if V~=V_alt
            disp('The two methods of obtaining V don''t agree.')
        end
        
        if size(V) ~= [nvar_VAR*T_VAR*nshocks,nvar_VAR*T_VAR *nshocks]
            disp('Size of bootstrap variances matrix is not correct.')
        end
        
        W = inv(V); %The weighting matrix
        W = W/norm(W);
    case 'no'
        disp 'We only have 1 single shock so we can quit V.'
        W = eye(nvar_VAR*nshocks*T_VAR);
end

disp 'Done up to weighting matrix.'


%%
%**********************************************************
% STAGE 4. Do fmincon
%**********************************************************
% What we wanna minimize here is some squared distance (IRFs_VAR - IRFs_Theory).
% I.e. we want param = argmin [IRFs_VAR - IRFs_Theory(param)].
% Thus for every new set of param values, fmincon also needs to resolve the
% model to get a new gx, hx and new IRFs_Theory. That part will be done in
% the objective function.

%Optimization Parameters
options = optimset('fmincon');
options = optimset(options, 'TolFun', 1e-9, 'display', 'iter');


%Selector matrix to select the jumps we're interested in % --->>>>> this is specific to BriantiGati2017
Sy = zeros(nvar_VAR,ny);
Sy(1,gamtfp_idx) = 1; % pretend like this is TFP for a sec
Sy(2,gamyi_idx) = 1; % IT investment
Sy(3,gamc_idx)  = 1; % consumption
Sy(4,gamp_idx)  = 1; % relative prices

%Selector matrix to select the states we're interested in
Sx = eye(size(hx,1)); % we're not selecting any states

% One-time evaluation of the objective
wlf = objective_IRmatching(param0,set,Sy,Sx,T_VAR,psi_hat,100000*W);

disp 'Evaluated loss once.'

% dbstop in objective_IRmatching if error
% dbstop in model_prog at 97 if isreal(objw(wx0)) == 0
% dbstop in objective_IRmatching at 5

%Objective with V weighting
objj = @(param) objective_IRmatching(param,set,Sy,Sx,T_VAR,psi_hat,100000*W);
%   [gam; sigma_IT; rho_IT]
LB = [.05,   .01,   0.1];
UB = [0.6,    100,   0.9];
% %   [a;        b;    gam; sigma_IT; rho_IT]
% LB = [0.01,  0.01,   .05,  .01,     0.1];
% UB = [0.5,   0.5,   0.6,    100,    0.9];
warning off
[param_opt,obj_opt] = fmincon(objj, param0,[],[],[],[],LB,UB,[],options);
% X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
warning on
toc
%% Generate IRFs with the optimal parameters
[f, fx, fy, fxp, fyp, G, R, set]=model_prog(param_opt,set); % --->>>>> this is specific to BriantiGati2017
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);

[~,ir_IT] = ir(gx,hx,G*1,H);  % IT prod level shock
ir_IT(:,gamc_idx:ny) = cumsum(ir_IT(:,gamc_idx:ny));
ir_matched = ir_IT(:,[gamtfp_idx gamyi_idx gamc_idx gamp_idx]);
ir_matched = ir_matched';

% Before printing the figs and so, re-add all the relevant paths because
% loading workspaces screws them up
current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

if exist([base_path '\Code'], 'dir')
    addpath([base_path '\Code\just_IT_vecm']) %for Microsoft
else
    addpath([base_path '/Code/just_IT_vecm']) %for Mac
end

print_figs ='no'; 
plot_matchedIRFs(IRF_IT,ir_matched,T_VAR,1,shocknames, {'TFP','IT inv', 'C', 'RP'}, print_figs, base_path, 'IRmatching_together')

param_names = {'gam', 'sigitlev', 'rhoitlev'};
% param_names = {'a', 'b', 'gam', 'sigitlev', 'rhoitlev'};
todays_date = datestr(today);

save_results = 'no';
switch save_results
    case 'yes'
        save IR_matching_results.mat param_names param_opt gx hx G ir_matched T_VAR IRF_IT base_path shocknames print_figs todays_date
end 


cd(current_dir)
disp('Done.')
return

%% Optional part: Load in saved workspaces and generate figures of SVAR and matched responses together
clear
load IR_matching_results
print_figs ='no'; 
plot_matchedIRFs(IRF_IT,ir_matched,T_VAR,1,shocknames, {'TFP','IT inv', 'C', 'RP'}, print_figs, base_path, 'IRmatching_together')
