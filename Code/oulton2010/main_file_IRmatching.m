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

% Add all the relevant paths
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

mod = model_IRmatching_spillover_news(param0,set);

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
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp); % No eq. exists?? WTF??? Problem is we're getting only 13 (instead of 27) stable eigs. 

% mom_tab shows a little table of stddevs, autocorrs and corrs of the
% specified variables after a specific shock (G = eta*shock vector)
mom_tab(gx,hx,G*G', [gamyc_idx,gamki_idx], {'YC','KI'})  %% here a dumb matrix size error 


%%
%**********************************************************
% STAGE 2. Input IRFs from VAR
%**********************************************************
load('Workspace_Just_IT_SVAR_1LAG.mat') % --->>>>> this is specific to BriantiGati2017

% Construct VAR IRFs to shocks of your choice
% Size of imported VAR IRFs is (nvar, T, nshocks)
nvar_VAR    = size(IRFs,1);
T_VAR       = size(IRFs,2);
nshocks_VAR = size(IRFs,3);
IRF_IT = IRFs(:,:,pos_IT); % response of all vars to IT shock % --->>>>> this is specific to BriantiGati2017

% Gather the VAR IRFs to all the relevant shocks
IRFs_VAR = [IRF_IT']; % --->>>>> this is specific to BriantiGati2017
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
IRFs_s = zeros(nvar_VAR, T_VAR, nshocks_VAR, nsimul);
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
IRFs_s_IT = squeeze(IRFs_s(:,:,pos_IT,:)); % --->>>>> this is specific to BriantiGati2017
psi_boot = reshape(IRFs_s_IT,nvar_VAR*T_VAR,200);
V = diag(var(psi_boot,0,2));

% alternative way to create V  -- both yield the same result, but inv(V) is
% singular, so something is not yet quite right. 
V1 = var(IRFs_s_IT,0,3); % take variance over nsimul
V2 = reshape(V1,size(V1,1)*size(V1,2),1);
V_alt = diag(V2);
if V~=V_alt
    disp('The two methods of obtaining V don''t agree.')
end

if size(V) ~= [nvar_VAR*T_VAR*nshocks,nvar_VAR*T_VAR*nshocks]
    disp('Size of bootstrap variances matrix is not correct.')
end

W = inv(V); %The weighting matrix
W = W/norm(W); 
    case 'no' 
       disp 'We only have 1 single shock so we can quit V.'
end

disp 'Done up to weighting matrix.'
return
%% Input gx hx from solved theoretical model and generate theoretical IRFs
%%%%%
% Cut out all of this part once the objective is correct.
%%%%%

load('gxhx.mat')
load('indexes.mat')
param = parameters;

% Construct theoretical IRFs
%Second Moments
nx = length(hx);
ny = size(gx,1);
eta = zeros(nx,nx);

eta(biggamitt_idx-ny,biggamitt_idx-ny) = param.sige; % --->>>>> this is specific to BriantiGati2017

shocks = zeros(nx,1);
shocks(biggamitt_idx-ny) = 1; % the noise shock  % --->>>>> this is specific to BriantiGati2017

%Selector matrix to select the variables we're interested in (only does this for jumps!!)
% which are for this project [biggamc_idx rc_idx it_idx gamyc_idx gamc_idx gamp_idx]
% Choose the variables such that they correspond to those in the VAR:
% Here I'm pretending that the variables in the model correspond to those
% in the VAR, which is not really the case and needs to be rethought.
% I wonder if this sets the order too...?
S = zeros(nvar_VAR-1,ny);
S(1,rc_idx)    = 1; 
S(2,it_idx)    = 1;
S(3,gamyc_idx) = 1;
S(4,gamc_idx)  = 1;
S(5,gamp_idx)  = 1;
g = S*gx;

%Impulse Responses Theoretical
[IRF_noise_Theory_all, IRFs_selected_jumps, IRF_states] = ir(g,hx,eta*shocks,T_VAR); %noise shock; % --->>>>> this is specific to BriantiGati2017

% I'm also selecting some states ... wonder if that's legitimate at all...?
IRFs_selected_states = IRF_states(:, [biggamc_idx-ny]); % --->>>>> this is specific to BriantiGati2017

% Gather the Theoretical IRFs to all the relevant shocks
IRFs_Theory = [IRFs_selected_states IRFs_selected_jumps]; % --->>>>> this is specific to BriantiGati2017
psi_T       = IRFs_Theory(:);


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

%%%%%
% Figure out what's going on with selecting states...?
%%%%%

%Selector matrix to select the jumps we're interested in % --->>>>> this is specific to BriantiGati2017
Sy = zeros(nvar_VAR-1,ny);
Sy(1,rc_idx)    = 1; 
Sy(2,it_idx)    = 1;
Sy(3,gamyc_idx) = 1;
Sy(4,gamc_idx)  = 1;
Sy(5,gamp_idx)  = 1;

%Selector matrix to select the states we're interested in 
% Sx

% One-time evaluation of the objective 
wlf = objective_IRmatching(param0,set,Sy,Sx,T_VAR,psi_hat,100000*W);
return
%Objective with V weighting
objj = @(param) objective_IRmatching(param,set,S,T_VAR,psi_hat,100000*W);
[param_opt,obj_opt] = fmincon(objj, param0,[],[],[],[],[.01,.01,.0,1.01,.0001,.0001],[.99,100,.99,15,1,1],[],options);



disp('Done.')