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

%% Correctly stationarized model w/o spillover
disp('Doing Model with NO Spillover')
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
nvar = size(gx,1) + size(gx,2);
nshocks = size(hx,1);
njumps  = size(gx,1);
pos_KC = kc_idx-njumps; % shock to K stock
pos_KI = ki_idx-njumps; % shock to IT stock
pos_BIGGAMC = biggamc_idx-njumps; % a shock to TFP in final goods prod (= surprise tech shock)
pos_BIGGAMI = biggami_idx-njumps; % IT productivity shock
T = 60;
IRFs_all = zeros(nvar, T, nshocks);
for s=1:nshocks
    x0 = zeros(nshocks,1); % impulse vector
    x0(s) = 1;
    [IR, iry, irx]=ir(gx,hx,x0,T);
    % To get levels, we need to cumsum all except RC, H, H1, H2
    IR(:,gamc_idx:njumps) = cumsum(IR(:,gamc_idx:njumps));
    IRFs_all(:,:,s) = IR'; 
end
% Gather IRFs of interest:
IRFs = IRFs_all(gamc_idx:njumps,:,:);
IRFs_some = IRFs_all([c_idx, gamc_idx],:,:);

which_shock = [pos_BIGGAMC];
shocknames = {'Hard capital shock', 'IT capital shock', 'Growth rate of Final', 'Growth rate of IT'};
% GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p
varnames = {'Logdev C', 'Logdev KI', 'Logdev YC', 'Logdev YI', 'Logdev H', 'Logdev P', 'Logdev KC2', 'Logdev KI2' };
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,shocknames, varnames, print_figs)


return

%% Model with spillover

disp('Doing Model with Spillover')

%Load Parameters
param = parameters;

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_spillover(param);

%Compute the transition and policy functions, using code by
%Stephanie Schmitt-Grohé and Martín Uribe (and available on their wedsite.)
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx))


% IRFs
nvar = size(gx,1) + size(gx,2);
nshocks = size(hx,1);
njumps  = size(gx,1);
pos_KC = kc_idx-njumps; % shock to K stock
pos_KI = ki_idx-njumps; % shock to IT stock
pos_BIGGAMC = biggamc_idx-njumps; % a shock to TFP in final goods prod (= surprise tech shock)
pos_BIGGAMCL = biggamcl_idx-njumps; % a shock to TFP in final goods prod but at a different time
pos_BIGGAMI = biggami_idx-njumps; % IT productivity shock
T = 60;
IRFs_all = zeros(nvar, T, nshocks);
for s=1:nshocks
    x0 = zeros(nshocks,1); % impulse vector
    x0(s) = 1;
    [IR, iry, irx]=ir(gx,hx,x0,T);
    % To get levels, we need to cumsum all except RC, H, H1, H2
    IR(:,gamc_idx:njumps) = cumsum(IR(:,gamc_idx:njumps));
    IRFs_all(:,:,s) = IR'; 
end
% Gather IRFs of interest:
IRFs = IRFs_all([biggamc_idx gamc_idx:njumps],:,:);
IRFs_some = IRFs_all([c_idx, gamc_idx],:,:);

which_shock = [pos_BIGGAMC];
shocknames = {'KC', 'KI', '\Gamma_c', '\Gamma_i'};
% GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p
varnames = {'Biggamc', 'C', 'KI', 'YC', 'YI', 'H', 'P', 'KC2', 'KI2' };
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,shocknames, varnames, print_figs)

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
nvar = size(gx,1) + size(gx,2);
nshocks = size(hx,1);
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
