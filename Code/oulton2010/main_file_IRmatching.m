% MAIN_FILE_IRMATCHING - Implements a general version of IR-matching a la CEE. The idea is that
% the structure should be general and you can use this file with modified
% functions or inputs in it to estimate different models.

% Code by Marco Brianti and Laura Gati, Boston College, 2018
%**************************************************************

%% % --->>>>> this is specific to BriantiGati2017
% This code specifically implements a matching of responses of an IT shock
% from our justIT VAR to the responses of our two-sector model with
% spillovers to a noise shock. It may be worth considering other shocks.

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

%% Input IRFs from VAR
load('Workspace_Just_IT_SVAR_1LAG.mat') % --->>>>> this is specific to BriantiGati2017

% Construct VAR IRFs to shocks of your choice
% Size of imported VAR IRFs is (nvar, T, nshocks)
nvar_VAR    = size(IRFs,1);
T_VAR       = size(IRFs,2);
nshocks_VAR = size(IRFs,3);
IRF_IT = IRFs(:,:,pos_IT); % response of all vars to IT shock % --->>>>> this is specific to BriantiGati2017

% Gather the VAR IRFs to all the relevant shocks
IRFs_VAR = [IRF_IT]; % --->>>>> this is specific to BriantiGati2017

%% Input gx hx from solved theoretical model and generate theoretical IRFs
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

%Impulse Responses Theory
IRF_noise_Theory = ir(gx,hx,eta*shocks,T_VAR); %noise shock; % --->>>>> this is specific to BriantiGati2017

% Choose the variables such that they correspond to those in the VAR:
% Here I'm pretending that the variables in the model correspond to those
% in the VAR, which is not really the case and needs to be rethought.
IRF_noise_Theory_subset = IRF_noise_Theory(:, [biggamc_idx rc_idx it_idx gamyc_idx gamc_idx gamp_idx]); % --->>>>> this is specific to BriantiGati2017

% Gather the Theoretical IRFs to all the relevant shocks
IRFs_Theory = [IRF_noise_Theory_subset']; % --->>>>> this is specific to BriantiGati2017

%% Do GMM
% What we wanna minimize here is some squared distance (IRFs_VAR - IRFs_Theory).
% I.e. we want param = argmin [IRFs_VAR - IRFs_Theory(param)].
% Thus for every new set of param values, fmincon also needs to resolve the
% model to get a new gx, hx and new IRFs_Theory. That part will be done in
% the objective function.

% Compute weighting matrix W = V^(-1), where V=var[bootstrap_IRFs] from the
% VAR.

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
V = var(IRFs_s_IT,0,3); %The weighting matrix - not quite working yet, the wrong size.
W = inv(V);
W = W/norm(W);

disp('Done.')