%**************************************************************
% MAIN_PROG_spillover_news - Solves the 2-sector model a la Oulton with spillover 
% and news shocks
%
% Code by Marco Brianti and Laura Gati, Boston College, 2018
%**************************************************************

clear 
close all

current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

%Load Parameters
param = parameters;

%Model with spillover

disp('Doing Model with Spillover and News')

%Load Parameters
param = parameters;

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_spillover_news(param);

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
pos_news    = v8_idx - njumps;
pos_N       = n_idx - njumps;
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
%IRFs_some = IRFs_all([c_idx, gamc_idx],:,:);

which_shock                = [pos_news];
shocknames                 = cell(1,size(hx,1));
shocknames(1,pos_KC)       = {'KC'};
shocknames(1,pos_KI)       = {'KI'};
shocknames(1,pos_BIGGAMC)  = {'\Gamma_c'};
shocknames(1,pos_BIGGAMI)  = {'\Gamma_i'};
shocknames(1,pos_news)     = {'News'};
shocknames(1, pos_N)     = {'Surprise common component'};

% GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p
varnames = {'Biggamc', 'C', 'KI', 'YC', 'YI', 'H', 'P', 'KC2', 'KI2' };
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,shocknames, varnames, print_figs)
