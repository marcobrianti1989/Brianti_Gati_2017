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

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model_spillover_news(param);

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
T = 40;
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

% GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p
varnames = {'C', 'KI', 'YC', 'YI', 'H', 'P', 'KC2', 'KI2', 'RI','W', 'GDP', 'TFP' };
print_figs = 'no';
plot_single_simple_IRFs(IRFs,T,which_shock,shocknames, varnames, print_figs, base_path)

% Try to construct NIPA-consistent GDP and TFP here
IRF_YC_level  = IRFs_all(gamyc_idx,:,:);
IRF_YI_level  = IRFs_all(gamyi_idx,:,:);
IRF_P_level   = IRFs_all(gamp_idx,:,:);
IRF_GDP_level       = IRF_YC_level + IRF_YI_level.*IRF_P_level;
plot_single_simple_IRFs(IRF_GDP_level,T,which_shock,shocknames, {'GDP'}, print_figs, base_path)

a = param.a; b = param.b;
IRF_RC        = IRFs_all(rc_idx,:,:);
IRF_RI_level  = IRFs_all(gamri_idx,:,:);
IRF_W_level   = IRFs_all(gamw_idx,:,:);
IRF_TFP_level = IRF_GDP_level - 2*(1-a-b)*IRF_W_level -2*a*IRF_RC -2*b*IRF_RI_level; 
plot_single_simple_IRFs(IRF_TFP_level,T,which_shock,shocknames, {'TFP'}, print_figs, base_path)
