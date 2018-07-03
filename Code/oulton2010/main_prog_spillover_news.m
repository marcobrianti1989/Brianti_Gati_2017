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
sim_shocks(pos_ITLEV) = 1; % turn on IT level shock
sim_shocks(pos_news)  = 1; % turn on news shock
[sim_data, all_shock_series] = simulate_model(gx,hx,eta*sim_shocks',T);
shocks_series = all_shock_series(:,[pos_news pos_ITLEV])';
% the question is whether we need to cumsum here (I guess so!)
sim_data(:,gamc_idx:njumps) = cumsum(sim_data(:,gamc_idx:njumps));
sim_data = sim_data';

plot(sim_data(gamtfp_idx,:), 'k'); hold on;
plot(shocks_series(1,:),'b')
plot(shocks_series(2,:),'r')
legend('sim. tfp', 'sim. news innovations', 'sim. IT innovations')


param = parameters;
ss = oulton_spillover_stst(param);
% New st. st. values:
global yc yi h1 h2 h kc1 kc2 kc ki1 ki2 ki ic it c w 

load ss_old
kc_old = xx(1); ki_old = xx(2);
yc_old = yy(1); yi_old = yy(2); c_old = yy(3); ic_old = yy(4); it_old = yy(5); w_old = yy(6);
h_old = yy(9); h1_old = yy(10); h2_old = yy(11); kc1_old = yy(12); kc2_old = yy(13); ki1_old = yy(14); ki2_old = yy(15);
% Old st st
% xx  = [kc ki biggamc biggami ...
%     c biggamc biggami yc yi h p kc2 ki2 ki kc ri w 0 0 0 0 0 0 0 0 0 ...
%     1 ...
%     biggami biggami ...
%     1 ...
%     h1 h2 kc1 ki1 it]; 
% yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
%     expgc expgi expgc expgi 1 expgc/expgi expgc expgi expgc/expgi expgc ...
%     gamgdp gamtfp expgc ...
%     1 1 expgc expgi expgi];
diffs = [yc-yc_old, yi-yi_old, h1-h1_old, h2-h2_old, h-h_old, kc1-kc1_old, kc2-kc2_old, kc-kc_old, ...
    ki1-ki1_old, ki2-ki2_old, ki-ki_old, ic-ic_old, it-it_old, c-c_old, w-w_old]