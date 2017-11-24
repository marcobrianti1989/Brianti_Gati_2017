% This file corresponds to main_file_SVAR, except it performs Ryan's
% identification strategy to identify news and IT shocks

% Marco Brianti, Laura Gáti, Oct 28 2017

% 500 bootstrap: Mac takes ? min
% 500 bootstrap: Server takes 30 min

clear
close all

tic

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:M286';
[data, varnames] = read_data2(filename, sheet, range);
shocknames = {'News Shock','IT Shock'};
which_shocks = [3 4];
pos_rel_prices = 6;

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 500; %5000
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

% dbstop in Ryan_two_stepsID at 74

% Implement Ryan's ID strategy
LR_hor = 8; % at what horizon to impose the LR restriction
% [impact, FEV_opt, IRFs, gamma_opt, FEV_news, FEV_IT] = ryansID(which_variable,which_shocks,H,B,A,q);
[impact, FEV_opt, ~, gam_opt, FEV_news, FEV_IT] ...
    = Ryan_two_stepsID(which_variable,which_shocks,H,LR_hor,B,A,pos_rel_prices);
% impact is the nvar x 2 impact matrix of news and IT.

% With Barsky & Sims-type ID, since you do a abs max, there are two
% solutions: gam and -gam. So choose the one that makes sense. I'm doing
% that so that news and IT have + effects on TFP.
if gam_opt(1,1) < 0 % If news have - effects on TFP
   impact(:,1) = -impact(:,1); 
end
if gam_opt(1,2) < 0 % if IT has - effects on TFP
    impact(:,2) = -impact(:,2);
end



% % TO DO: check info sufficiency (Forni's orthogonality test)
% [s, obj_opt] = get_structral_shocks_alternative(A,gam_opt,res);
% s3 = s(3,:)';
% s4 = s(4,:)';
% % [s_old, obj_opt] = get_structural_shocks(A,gam_opt,res); % this old way
% % didn't work well, cancel it completely after a while
% % s1 = s_old(1,:)';
% % s2 = s_old(2,:)';
%
% pc = get_principal_components(data);
% pc = pc(1+nlags:end,:); % adjust because due to the lags, the structural shocks are shorter
%
% [B, yhat, res] = quick_ols(s3,pc);
% [B, yhat, res] = quick_ols(s4,pc);
%
% signi = f_test(s4,res, nvar,156,nvar)
%
% return





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
    [impact_boot(:,:,i_simul),~,~,~,~,~] = Ryan_two_stepsID(which_variable,which_shocks,H,LR_hor, ...
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
comment = [which_ID '_' char(varnames(6)) '_LR_hor_' num2str(LR_hor)];

print_figs = 'no';
[IRFs, ub, lb] = genIRFs(fake_impact,fake_impact_boot,B,beta_tilde_star,H,sig);
plotIRFs(IRFs,ub,lb,40,which_shocks,shocknames,varnames, which_ID,print_figs)

fev_matrix = {'News', 'IT', 'Total'};
fev_matrix(2,:) = {num2str(FEV_news), num2str(FEV_IT), num2str(FEV_opt)};

disp('% of FEV of TFP explained:')
fev_matrix

export_FEV_matrix = 'yes';
if strcmp(export_FEV_matrix,'yes') ==1
    fev_matrix_out = [FEV_news, FEV_IT, FEV_opt];
    rowLabels = {'Share of TFP FEV explained'};
    columnLabels = {'News', 'IT', 'Total'};
    matrixname   = 'FEVs';
    invoke_matrix_outputting(fev_matrix_out,matrixname,rowLabels,columnLabels,comment);
end




%save('capital_price_workspace')

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
disp(varnames)
disp(datestr(now))

