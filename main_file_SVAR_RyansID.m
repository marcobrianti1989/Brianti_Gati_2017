% This file corresponds to main_file_SVAR, except it performs Ryan's
% identification strategy to identify news and IT shocks

% Marco Brianti, Laura Gáti, Oct 28 2017

clear
close all

tic

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:W287';
[data, varnames] = read_data2(filename, sheet, range);
shocknames = {'News Shock','IT Shock'};
varnames

% Find positions of shocks and of relative prices in a correct way
if sum(strcmp('Real SP', varnames))==1 % i.e. 'Real SP' exists as a variable
    pos_news = find(strcmp('Real SP', varnames));
elseif sum(strcmp('Michigan Index ', varnames))==1 % i.e. we use Mich for news
    pos_news = find(strcmp('Michigan Index ', varnames));
elseif sum(strcmp('SP deflated', varnames))==1 % i.e. we use SP deflated by GDPDEF for news    
    pos_news = find(strcmp('SP deflated', varnames));
elseif sum(strcmp('SP deflated per capita', varnames))==1 % i.e. we use SP deflated by GDPDEF for news    
    pos_news = find(strcmp('SP deflated per capita', varnames));
end
pos_IT = find(strcmp('Real IT Investment', varnames));
if find(strcmp('Relative Price', varnames)) > 0
    pos_rel_prices = find(strcmp('Relative Price', varnames));
elseif find(strcmp('Relative price PCE', varnames))
    pos_rel_prices = find(strcmp('Relative price PCE', varnames));
end
which_shocks = [pos_news pos_IT];

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 30; %5000
nvar            = size(data,2);
sig1            = 0.9; % significance level
sig2            = 0.95; % a 2nd sig. level
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
H_max = 60;
[impact, FEV_opt, ~, gam_opt, FEV_news, FEV_IT] ...
    = Ryan_two_stepsID(which_variable,which_shocks,H_max,...
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
          which_shocks,H_max,LR_hor,beta_tilde_star(:,:,i_simul),A_boot, pos_rel_prices);
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
      B,beta_tilde_star,H,sig1, sig2);

% % With Barsky & Sims-type ID, since you do a abs max, there are two
% % solutions: gam and -gam. So choose the one that makes sense. I'm doing
% % that so that news and IT have + effects on TFP.
if sum(IRFs(1,:,pos_IT)) < 0 % if the majority of TFP response is negative
    IRFs(:,:,pos_IT)   = -IRFs(:,:,pos_IT);
    ub1(:,:,pos_IT)     = - ub1(:,:,pos_IT);
    lb1(:,:,pos_IT)     = - lb1(:,:,pos_IT);
    ub2(:,:,pos_IT)     = - ub2(:,:,pos_IT);
    lb2(:,:,pos_IT)     = - lb2(:,:,pos_IT);
end
if sum(IRFs(1,:,pos_news)) < 0 % if the majority of TFP response is negative
    IRFs(:,:,pos_news)   = -IRFs(:,:,pos_news);
    ub1(:,:,pos_news)     = - ub1(:,:,pos_news);
    lb1(:,:,pos_news)     = - lb1(:,:,pos_news);
    ub2(:,:,pos_news)     = - ub2(:,:,pos_news);
    lb2(:,:,pos_news)     = - lb2(:,:,pos_news);
end

%Printing/Showing IRFs
h = 40;
% plotIRFs(IRFs,ub,lb,40,which_shocks,shocknames,varnames,which_ID,print_figs)
% plot_single_IRFs(IRFs,ub1,lb1,h,which_shocks,shocknames, varnames, which_ID, print_figs)
use_current_time = 0; % don't save the time
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,h,which_shocks,shocknames, varnames, '_', print_figs, use_current_time)

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

% Counterfactual
[s, IR] = get_structural_shocks_Forni(A,gam_opt,res);

which_shock_zero = 3; % set negative realizations of the IT shock to 0
y0 = data(1,:); % initial values of data
counterfactuals = counterfactual(s,IR,B,y0,which_shock_zero);

time = datetime(1989,7,1) + calquarters(0:size(s,1)); %ALWAYS CHECK IF THE INITIAL PERIOD IS THE SAME

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
plot(time, data(:,1), 'b', 'linewidth', 2)
hold on
plot(time, counterfactuals(:,1), 'r', 'linewidth', 2)
legend('Original TFP series', 'Counterfactual TFP','Location','Southeast')
title('Shutting off negative IT productivity shocks')
set(gca, 'FontSize', 28)
grid on

save_counterfactual_fig = 'yes';
switch save_counterfactual_fig
    case 'yes'
        invoke_export_fig('counterfactual','',0)
end

% % Historical decompositions of TFP from the 2 structral shocks
% % IT:
% hd_IT = historical_decomposition(s(:,pos_IT),fake_impact,B, pos_IT, 1);
% % News:
% hd_news = historical_decomposition(s(:,pos_news),fake_impact,B, pos_news, 1);
% % sum of the two:
% hd = hd_IT + hd_news;
% 
% % Get counterfactual when shutting off IT shocks 2000-Q3 onward
% IT_bubble = find(time=='01-Jul-2000');
% s_alternative = s;
% s_alternative(IT_bubble:end,3) = 0;
% % Redo counterfactuals for alternative scenario
% % IT:
% hd_IT_alt = historical_decomposition(s_alternative(:,pos_IT),fake_impact,B, pos_IT, 1);
% % News:
% hd_news_alt = historical_decomposition(s_alternative(:,pos_news),fake_impact,B, pos_news, 1);
% % sum of the two:
% hd_alt = hd_IT_alt + hd_news_alt;





