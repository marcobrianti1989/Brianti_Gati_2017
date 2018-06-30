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


%Model with spillover

disp('Doing Model with Spillover and News')

%Load Parameters
param = parameters;
% param = parameters_opt;

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_spillover_news(param);
% save checking_f.mat fxn  % save fxn so you can compare other fx's obtained by other codes.

%Compute the transition and policy functions, using code by
%Stephanie Schmitt-Grohé and Martín Uribe (and available on their wedsite.)
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')
save('indexes.mat', '*_idx')

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
pos_BIGGAMI = biggami_idx-njumps; % IT productivity growth shock, contemporaneous
pos_ITLEV   = itlev_idx-njumps; % IT productivity level shock, contemporaneous
pos_news    = v8_idx - njumps;
pos_N       = n_idx - njumps;
pos_noise   = biggamitt_idx - njumps; %preliminary noise shock
T = 100;
eta = eye(nshocks);
eta(pos_BIGGAMI,pos_BIGGAMI) = param.siggami;
eta(pos_ITLEV,pos_ITLEV)     = param.sigitlev;
IRFs_all = zeros(nvar, T, nshocks);
for s=1:nshocks
    x0 = zeros(nshocks,1); % impulse vector
    x0(s) = 1;
    [IR, iry, irx]=ir(gx,hx,eta*x0,T);
    % To get levels, we need to cumsum all except RC, H, H1, H2
    IR(:,gamc_idx:njumps) = cumsum(IR(:,gamc_idx:njumps));
    IRFs_all(:,:,s) = IR';
end
% Gather IRFs of interest:
IRFs = IRFs_all(gamc_idx:gamtfp_idx,:,:);
IRFs_TFP_GDP = IRFs_all([gamgdp_idx, gamtfp_idx],:,:);
IRFs_VAR     = IRFs_all([gamtfp_idx, gamit_idx, gamc_idx, gamp_idx],:,:);

which_shock                = pos_ITLEV;
shocknames                 = cell(1,size(hx,1));
shocknames(1,pos_KC)       = {'KC'};
shocknames(1,pos_KI)       = {'KI'};
shocknames(1,pos_BIGGAMC)  = {'\Gamma_c'};
shocknames(1,pos_BIGGAMI)  = {'\Gamma_i'};
shocknames(1,pos_news)     = {'News'};
shocknames(1, pos_N)       = {'Surprise common component'};
shocknames(1, pos_noise)   = {'Noise'};
shocknames(1, pos_ITLEV)   = {'IT productivity (level)'};

print_figs = 'no';

% GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p
varnames = {'C', 'KI', 'YC', 'YI', 'H', 'P', 'KC2', 'KI2', 'RI','W', 'GDP', 'TFP' };
varnames_matching = {'TFP','IT inv', 'C', 'RP'};
disp('Done.')
return
plot_single_simple_IRFs(IRFs,T,which_shock,shocknames, varnames, print_figs, base_path)

% GDP and TFP IRFs alone
plot_single_simple_IRFs(IRFs_TFP_GDP,T,which_shock,shocknames, {'GDP', 'TFP'}, print_figs, base_path)

% 'TFP','IT inv', 'C', 'RP' IRFs only
plot_single_simple_IRFs(IRFs_VAR,T,which_shock,shocknames, varnames_matching, print_figs, base_path)

save IRFs_gamma_opt.mat IRFs_VAR T which_shock shocknames varnames_matching print_figs base_path

% Generate simulated data from model
sim_shocks = zeros(1,nshocks);
sim_shocks(pos_ITLEV) = eta(pos_ITLEV,pos_ITLEV);
sim_data = simulate_model(gx,hx,sim_shocks,T);
% the question is whether we need to cumsum here (I guess so!)
sim_data(:,gamc_idx:njumps) = cumsum(sim_data(:,gamc_idx:njumps));
sim_data = sim_data';

plot(sim_data(gamtfp_idx,:))


param = parameters;
ss = oulton_spillover_stst(param);