% This file corresponds to main_file_SVAR, except it performs Ryan's
% identification strategy to identify news and IT shocks

% Marco Brianti, Laura Gáti, Oct 28 2017

% TOOK 72.1348 MIN TO RUN!!

clear 
close all

tic
%Data Reading and Transformation

filename = 'dataset_23_sept_2017';
sheet    = 'Sheet1';
range    = 'B126:L283';
[data, shocknames,varnames, which_shocks, pos_rel_prices] = ...
    read_data(filename, sheet, range); %q (pos_rel_prices) position for relative prices

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 20; %5000
nvar            = size(data,2);
sig             = 0.90; % significance level
H               = 100; %40; % horizon for generation of IRFs
h               = H; %40; % horizon for IRF plots
which_variable  = 1; % select TFP as the variable whose FEV we wanna max

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
    nlags = 1; % 1
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run VAR imposing Cholesky
[A,B,res,sigma] = sr_var(data, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% TO DO: check info sufficiency


% dbstop in Ryan_two_stepsID at 74

% Implement Ryan's ID strategy
% [impact, FEV_opt, IRFs, gamma_opt, FEV_news, FEV_IT] = ryansID(which_variable,which_shocks,H,B,A,q);
[impact, FEV_opt, ~, gam3_opt, FEV_news, FEV_IT] ...
    = Ryan_two_stepsID(which_variable,which_shocks,H,B,A,pos_rel_prices);
% impact is the nvar x 2 impact matrix of news and IT.

% Bootstrap
which_ID = 'Ryan_two_stepsID';
which_correction = 'blocks'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks

[beta_tilde, data_boot2, beta_tilde_star, nonstationarities] = bootstrap_with_kilian( ...
    B, nburn, res, nsimul, which_correction, blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
    [A_boot, ~,~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
    % Get bootstrapped confidence intervals nsimul times
    disp(['Iteration ' num2str(i_simul) ' out of ' num2str(nsimul)])
    [impact_boot(:,:,i_simul),~,~,~,~,~] = Ryan_two_stepsID(which_variable,which_shocks,H, ...
        beta_tilde_star(:,:,i_simul),A_boot, pos_rel_prices);
end

%Creating a fake matrix for the IRF of the point estimation
fake_impact = zeros(nvar,nvar);
fake_impact(:,which_shocks) = impact;
% Create a fake matrix for the IRFs of the bootstrap
fake_impact_boot = zeros(nvar,nvar,nsimul);
for i_simul = 1:nsimul
    fake_impact_boot(:,which_shocks,i_simul) = impact_boot(:,:,i_simul);
end

%Creating and Printing figures
print_figs = 'no';
[IRFs, ub, lb] = genIRFs(fake_impact,fake_impact_boot,B,beta_tilde_star,H,sig);
plotIRFs(IRFs,ub,lb,40,which_shocks,shocknames,varnames, which_ID,print_figs)


% %%%%%%% Saved work from before really implementing Kilian, uncomment all
% %%%%%%% of this if the real Kilian stuff are superbad - although the two
% seem to be working nearly identically.
% [A_boot, B_boot] = ...
%     bootstrappedCIs(B, nburn, res, nsimul, which_correction, blocksize, nvar, ...
%     nlags,pos_rel_prices, which_shocks,which_variable,H,which_ID,impact); % automatically does Kilian correction
% 
% %Creating a fake matrix for the IRF
% fake_impact = zeros(nvar,nvar);
% fake_impact(:,which_shocks) = impact;
% 
% %Creating and Printing figures
% print_figs = 'yes';
% [IRFs, ub, lb] = genIRFs(fake_impact,A_boot,B,B_boot,H,sig);
% plotIRFs(IRFs,ub,lb,h,which_shocks,shocknames,varnames, which_ID,print_figs)
% %%%%%%
toc

