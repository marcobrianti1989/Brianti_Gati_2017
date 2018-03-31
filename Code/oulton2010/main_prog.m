%**************************************************************
% MAIN_PROG - Solves the neoclassical model with a random walk
% TFP solved in class in EC860.
%
% Code by Ryan Chahrour, Boston College, 2014
%**************************************************************

clear all

current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

%Load Parameters
param = parameters;

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_exog_stat(param);

%Compute the transition and policy functions, using code by
%Stephanie Schmitt-Grohé and Martín Uribe (and available on their wedsite.)
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx))


% IRFs
% Positions of the shocks in shock vector:
nvar = size(gx,1) + size(gx,2);
nshocks = size(hx,1);
njumps  = size(gx,1);
pos_KC = 1; % shock to K stock
pos_KI = 2; % shock to IT stock
pos_BIGGAMC = 3; % a shock to TFP in final goods prod (= surprise tech shock)
pos_BIGGAMI = 4; % IT productivity shock
T = 200;
IRFs = zeros(nvar, T, nshocks);
for s=3 %1:nshocks
    x0 = zeros(nshocks,1); % impulse vector
    x0(s) = 1;
    [IR, iry, irx]=ir(gx,hx,x0,T);
    % To get levels, we need to cumsum all except RC, H, H1, H2
    IR(:,[1:6,8,12:end]) = cumsum(IR(:,[1:6,8,12:end]));
    IRFs(:,:,s) = IR'; 
end

which_shock = [pos_BIGGAMC];
% [KC KI BIGGAMC BIGGAMI CL KIL]
names = {'Hard capital shock', 'IT capital shock', 'Growth rate of Final', 'Growth rate of IT'};
% Y  = [YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P GAMC GAMKI]
varnames = {'YC', 'YI', 'C', 'IC', 'IT', 'W', 'RC', 'RI', 'H', 'H1', 'H2', 'KC1', 'KC2', 'KI1', 'KI2','P', 'Logdev C', 'Logdev KI', ...
     'KC', 'KI', 'BIGGAMC', 'BIGGAMI', 'CL', 'KIL'};
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,names, varnames, print_figs)

return
%% Oulton 2010 with news shocks (so far only 2 periods ahead)

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_news(param);

%Compute the transition and policy functions
[gx_v,hx_v]=gx_hx_alt(fyn,fxn,fypn,fxpn);

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx_v))

% IRFs
% Positions of the shocks in shock vector:
nvar = size(gx_v,1) + size(gx_v,2);
nshocks = size(hx_v,1);
pos_KC = 1; % shock to K stock
pos_KI = 2; % shock to IT stock
pos_GAMC = 3; % a shock to TFP in final goods prod (= surprise tech shock)
pos_GAMI = 4; % IT productivity shock
pos_news1 = 5; % news is confirmed
pos_news2 = 6; % news initially arrives

x0 = [0 0 0 0 0 0]'; % impulse vector
T = 60;
IRFs = zeros(nvar, T, nshocks);
for s=1:nshocks
    x0(s) = 1;
    [IR, iry, irx]=ir(gx_v,hx_v,x0,T);
    IRFs(:,:,s) = IR'; % 
end

which_shock = [pos_news2];
names = {'KC', 'KI', 'GAMC', 'GAMI', 'news confirmed', 'news arrives'};
varnames = {'YC', 'YI', 'C', 'IC', 'IT', 'W', 'RC', 'RI', 'H', 'H1', 'H2', 'KC1', 'KC2', 'KI1', 'KI2', ...
    'G', 'GI', 'GP', 'KC', 'KI', 'GAMC', 'GAMI', 'V1', 'V2'};
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,names, varnames, print_figs)
