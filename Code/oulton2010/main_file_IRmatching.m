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

%% Input gx hx from solved theoretical model
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

% Gather the Theoretical IRFs to all the relevant shocks
IRFs_Theory = [IRF_noise_Theory]; % --->>>>> this is specific to BriantiGati2017

%% Do GMM

disp('Done.')