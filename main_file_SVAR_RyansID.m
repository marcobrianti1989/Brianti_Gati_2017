% This file corresponds to main_file_SVAR, except it performs Ryan's
% identification strategy to identify news and IT shocks

% Marco Brianti, Laura Gáti, Oct 28 2017

clear
close all

tic

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:P286';
[data, varnames] = read_data2(filename, sheet, range);
shocknames = {'News Shock','IT Shock'};
% which_shocks = [3 4];
which_shocks   = [2 3];
pos_rel_prices = 6;
pos_news       = which_shocks(1);
pos_IT         = which_shocks(2);
if pos_rel_prices == 6
else
      warning('Position of Relative Price is not anymore 6.')
end
if strcmp(varnames(pos_news),'Real SP') == 1
    % good
else
    warning('Position of News is not anymore 2.')
end
if strcmp(varnames(pos_IT),'Real IT Investment') == 1
    % good
else
    warning('Position of IT is not anymore 3.')
end
%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 20; %5000
nvar            = size(data,2);
sig1            = 0.80; % significance level
sig2            = 0.90; % significance level
H               = 100; %40; % horizon for generation of IRFs
h               = H; %40; % horizon for IRF plots
which_variable  = 1; % select TFP as the variable whose FEV we wanna max

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
    nlags = 1; % 2;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run VAR imposing Cholesky
[A,B,res,sigma] = sr_var(data, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% Implement Ryan's ID strategy
LR_hor = 8; % at what horizon to impose the LR restriction
% [impact, FEV_opt, IRFs, gamma_opt, FEV_news, FEV_IT] = ryansID(which_variable,which_shocks,H,B,A,q);
[impact, FEV_opt, ~, gam_opt, FEV_news, FEV_IT] ...
    = Ryan_two_stepsID(which_variable,which_shocks,H,...
    LR_hor,B,A,pos_rel_prices);

% Bootstrap
which_ID = 'Ryan_two_stepsID';
which_correction = 'blocks'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks
[beta_tilde, data_boot2, beta_tilde_star, nonstationarities] = ...
      bootstrap_with_kilian(B, nburn, res, ...
      nsimul, which_correction, blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
    [A_boot, ~,~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
    % Get bootstrapped confidence intervals nsimul times
    disp(['Iteration ' num2str(i_simul) ' out of ' num2str(nsimul)])
    [impact_boot(:,:,i_simul),~,~,~,~,~] = Ryan_two_stepsID(which_variable,...
          which_shocks,H,LR_hor,beta_tilde_star(:,:,i_simul),A_boot, pos_rel_prices);
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
comment = [which_ID '_' char(varnames(5)) '_LR_hor_' num2str(LR_hor)];
print_figs = 'no';

[IRFs, ub1, lb1, ub2, lb2] = genIRFs(fake_impact,fake_impact_boot,...
      B,beta_tilde_star,H,sig1,sig2);

% % With Barsky & Sims-type ID, since you do a abs max, there are two
% % solutions: gam and -gam. So choose the one that makes sense. I'm doing
% % that so that news and IT have + effects on TFP.
% if sum(IRFs(1,:,pos_news)) < 0 % if the majority of TFP response is negative (news)
%     IRFs(:,:,pos_news)   = -IRFs(:,:,pos_news);
%     ub1(:,:,pos_news)    = - ub1(:,:,pos_news);
%     lb1(:,:,pos_news)    = - lb1(:,:,pos_news);
%     ub2(:,:,pos_news)    = - ub2(:,:,pos_news);
%     lb2(:,:,pos_news)    = - lb2(:,:,pos_news);
% end
% if sum(IRFs(1,:,pos_IT)) < 0 % if the majority of TFP response is negative (IT)
%     IRFs(:,:,pos_IT)   = -IRFs(:,:,pos_IT);
%     ub1(:,:,pos_IT)    = - ub1(:,:,pos_IT);
%     lb1(:,:,pos_IT)    = - lb1(:,:,pos_IT);
%     ub2(:,:,pos_IT)    = - ub2(:,:,pos_IT);
%     lb2(:,:,pos_IT)    = - lb2(:,:,pos_IT);
% end

% plotIRFs(IRFs,ub1,lb1,40,which_shocks,shocknames,varnames, which_ID,print_figs)

%Printing/Showing IRFs
h = 100;
% plotIRFs(IRFs,ub1,lb1,40,which_shocks,shocknames,varnames, ...
%       which_ID,print_figs)
plot_single_IRFs(IRFs,ub1,lb1,ub2,lb2,h,which_shocks,...
    shocknames, varnames, which_ID, print_figs)

%Forni&Gambetti Orthogonality Test
do_FG_test = 'no';
switch do_FG_test
    case 'yes'
filename_PC       = 'Dataset_test_PC';
sheet_PC          = 'Quarterly';
range_PC          = 'B2:DC287';
first_n_PCs       = 10;
mlags             = 1;
[pvalue_news_shock, pvalue_IT_shock] = ...
      Forni_Gambetti_orthogonality_test(filename_PC,...
      sheet_PC,range_PC,first_n_PCs,A,gam_opt,res,which_shocks,mlags)
end

%Saving in Tex format the Variance Decomposition Matrix
fev_matrix = {'News', 'IT', 'Total'};
fev_matrix(2,:) = {num2str(FEV_news), num2str(FEV_IT), num2str(FEV_opt)};
disp('% of FEV of TFP explained:')
fev_matrix

export_FEV_matrix = 'no';
if strcmp(export_FEV_matrix,'yes') ==1
    fev_matrix_out = [FEV_news, FEV_IT, FEV_opt];
    rowLabels = {'Share of TFP FEV explained'};
    columnLabels = {'News', 'IT', 'Total'};
    matrixname   = 'FEVs';
    invoke_matrix_outputting(fev_matrix_out,matrixname,rowLabels,...
          columnLabels,comment);
end

toc
disp(varnames)
disp(datestr(now))

